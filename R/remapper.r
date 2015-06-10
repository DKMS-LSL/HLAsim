#' @rdname make_mapper
#' @export
make_remapper <- function(lookup) {
  ## global funcs
  detect  <- stringi::stri_detect_regex
  `%do%`  <- foreach::`%do%`
  memoise <- memoise::memoise

  expand_code <- function(x) {
    if (is_nmdp(x)) {
      remap_nmdp_code(x, ntbl4)
    } else if (is_g_code(x)) {
      remap_g_code(x, gtbl)
    } else x
  }

  ptn <- function(x) paste0("^(", paste0(x, collapse = "|"), ")(:.*)?$")

  get_partials <- function(ex3, prt) {
    anum1 <- ex3[, .(allele_num = allele_num[1]), by = eag_num]
    setkeyv(anum1, "allele_num")
    prt2 <- na.omit(prt[anum1, .(eag_allele, eag_num, allele_num), allow.cartesian = TRUE])
    setkeyv(prt2, "eag_num")
    prt2
  }

  crosscheck_nmdp_subtypes <- function(rep) {
    prt3.2 <- copy(prt3[detect(eag_allele, rep)])
    prt3.2[, allele_num := NULL]
    prt3.2[, eag_allele := vapply(hla_burst(eag_allele), function(x) hla_collapse(x[1:2]), "")]
    prt3.2 <- unique(prt3.2)
    Map(function(nums) enp[eag_num == nums][, prob := prob/sum(prob)],
        nums = split(prt3.2$eag_num, prt3.2$eag_allele))
  }

  crosscheck_eag_nums <- function(rep) {
    nums <- unique(prt3[detect(eag_allele, rep), eag_num])
    enp[eag_num == nums][, prob := prob/sum(prob)]
  }

  ## global objects
  ex2 <- exon(lookup, 2)[, .(eag_allele, eag_num, allele_num)]
  ex2 <- na.omit(ex2[!duplicated(ex2)])
  setkeyv(ex2, "eag_allele")

  ex3 <- exon(lookup, 3)[, .(eag_allele, eag_num, allele_num)]
  ex3 <- na.omit(ex3[!duplicated(ex3)])
  setkeyv(ex3, "eag_allele")

  jkr   <- lookup$jkr_tbl$eag_allele
  prt3  <- rbind(ex3, get_partials(ex3, prt = lookup$prt_tbl))
  gtbl  <- lookup$gtbl
  ntbl4 <- lookup$ntbl4

  ## DT of unique exon 3 EAGs and their approximate probability
  ## of accurance simply based on how many alleles they stand
  ## for: eag_num_prob => enp
  enp <- ex3[, .(prob = .N/nrow(ex3)), by = eag_num]

  structure(memoise(function(a) {
    if (is.na(x = a))
      return(eag_numbers(NA, NA))

    a1   <- allele(a, 1)
    exa1 <- expand_code(a1)
    rep1 <- ptn(exa1)

    a2   <- allele(a, 2)
    exa2 <- expand_code(a2)
    rep2 <- ptn(exa2)

    nums2.1 <- unique(ex2[detect(eag_allele, rep1), eag_num])
    nums2.2 <- unique(ex2[detect(eag_allele, rep2), eag_num])

    ## Exon 3 gets special treatment due to the fact that some
    ## alleles may not be found in the lookup table but only
    ## exist as jokers.
    nums3.1 <- ex3[detect(eag_allele, rep1), eag_num]
    nums3.2 <- ex3[detect(eag_allele, rep2), eag_num]

    # If both alleles on exon 3 have a single EAG number
    # each, no further checks are required.
    if (not_single(nums3.1) || not_single(nums3.2)) {
      # Allele 1 and allele 2 have one or more EAG nums.
      if (none_empty(nums3.1, nums3.2)) {
        if (a1 == "05:RGPW") {
          nums3.1 <- unique(nums3.1[duplicated(nums3.1)])
        }
        if (a2 == "05:RGPW") {
          nums3.2 <- unique(nums3.2[duplicated(nums3.2)])
        }
        # Both alleles have ambiguity codes and an intersection
        # resulting in two distinct EAG numbers can be formed.
        if (maybe_exon_shuffling(a) && length(nums3 <- intersect(nums3.1, nums3.2)) == 2) {
          nums3.1 <- nums3[1]
          nums3.2 <- nums3[2]
        }
        else {
          if (length(nums3.1 <- unique(nums3.1)) > 1) {
            #message("EAG nums ", slash(nums3.1), " in allele 1 of ", slash(a))
            enp2 <- crosscheck_eag_nums(rep = rep1)
            nums3.1 <- sample(enp2$eag_num, 1, prob = enp2$prob)
          }
          if (length(nums3.2 <- unique(nums3.2)) > 1) {
            #message("EAG nums ", slash(nums3.2), " in allele 2 of ", slash(a))
            enp2 <- crosscheck_eag_nums(rep = rep2)
            nums3.2 <- sample(enp2$eag_num, 1, prob = enp2$prob)
          }
        }
      }
      # Only allele 1 has one or more EAG numbers, that is,
      # allele 2 is part of the joker or partials group.
      #
      # Since EAGs are defined with respect to exon 2, the subtypes of
      # an ambiguous allele on exon 3 can be spread across 2 EAGs. If
      # this is the case we pick these two EAGs
      #
      # If this is not the case we randomly pick a second EAG in
      # proportion to the numbers of alleles associated with an EAG.
      #
      # We have to make sure, however, that there is no intersection
      # (outside of the jokers/partials) between this EAG and the alleles at
      # exon 2
      else if (is_empty(nums3.2) && !is_empty(nums3.1)) {
        if (a1 == "05:RGPW") {
          nums3.1 <- nums3.2 <- unique(nums3.1[duplicated(nums3.1)])
        }
        else {
          if (length(nums3 <- unique(nums3.1)) == 2) {
            nums3.1 <- nums3[1]
            nums3.2 <- nums3[2]
          }
          else if (length(nums3) > 2) {
            #message("EAG nums ", slash(nums3), " in allele 1 of ", slash(a))
            enp2 <- crosscheck_nmdp_subtypes(rep = rep1)
            if (length(enp2) == 2) {
              nums3.1 <- sample(enp2[[1]]$eag_num, 1, prob = enp2[[1]]$prob)
              nums3.2 <- sample(enp2[[2]]$eag_num, 1, prob = enp2[[2]]$prob)
            }
            else {
              enp2 <- crosscheck_eag_nums(rep = rep1)
              nums3 <- sample(enp2$eag_num, 2, replace = TRUE, prob = enp2$prob)
              nums3.1 <- nums3[1]
              nums3.2 <- nums3[2]
            }
          }
          else {
            nums3.1 <- nums3
            enp2 <- if (!all(a2 %in% jkr)) {
              crosscheck_eag_nums(rep = rep2)
            } else if (length(aa <- ex2[eag_num == nums2.1][!detect(eag_allele, rep1), eag_allele]) > 0) {
              ## check for additional diagonal allele matches
              rs <- enp[eag_num == setdiff(enp$eag_num, crosscheck_eag_nums(rep = ptn(aa))$eag_num)]
              rs[, prob := prob/sum(prob)]
            } else enp
            nums3.2 <- sample(enp2$eag_num, 1, prob = enp2$prob)
          }
        }
      }
      # Only allele 2 has one or more EAG numbers, that is,
      # allele 1 is part of the joker/partials group.
      #
      # Since EAGs are defined with respect to exon 2, the subtypes of
      # an ambiguous allele on exon 3 can be spread across 2 EAG. If
      # this is the case we pick these two EAGs.
      #
      # If this is not the case we randomly pick a second EAG in
      # proportion to the numbers of alleles associated with an EAG.
      #
      # We have to make sure, however, that there is no intersection
      # (outside of the jokers/partials) between this EAG and the alleles at
      # exon 2
      else if (is_empty(nums3.1) && !is_empty(nums3.2)) {
        if (a2 == "05:RGPW") {
          nums3.1 <- nums3.2 <- unique(nums3.2[duplicated(nums3.2)])
        }
        else {
          if (length(nums3 <- unique(nums3.2)) == 2) {
            nums3.1 <- nums3[1]
            nums3.2 <- nums3[2]
          }
          else if (length(nums3) > 2) {
            #message("EAG nums ", slash(nums3), " in allele 2 of ", slash(a))
            enp2 <- crosscheck_nmdp_subtypes(rep = rep2)
            if (length(enp2) == 2) {
              nums3.1 <- sample(enp2[[1]]$eag_num, 1, prob = enp2[[1]]$prob)
              nums3.2 <- sample(enp2[[2]]$eag_num, 1, prob = enp2[[2]]$prob)
            }
            else {
              enp2 <- crosscheck_eag_nums(rep = rep2)
              nums3 <- sample(enp2$eag_num, 2, replace = TRUE, prob = enp2$prob)
              nums3.1 <- nums3[1]
              nums3.2 <- nums3[2]
            }
          }
          else {
            nums3.2 <- nums3
            enp2 <- if (!all(a1 %in% jkr)) {
              crosscheck_eag_nums(rep = rep1)
            } else if (length(aa <- ex2[eag_num == nums2.2][!detect(eag_allele, rep2), eag_allele]) > 0) {
              ## check for additional diagonal allele matches
              rs <- enp[eag_num == setdiff(enp$eag_num, crosscheck_eag_nums(rep = ptn(aa))$eag_num)]
              rs[, prob := prob/sum(prob)]
            } else enp

            nums3.1 <- sample(enp2$eag_num, 1, prob = enp2$prob)
          }
        }
      }
      # Neither allele has an EAG number because both alleles
      # are part of the joker group
      else {
        ## Randomly pick EAG numbers in proportion to the
        ## numbers of alleles associated with a number. Because
        ## both alleles are part of the joker group this should
        ## be safe.
        enp2 <- if (!all(a1 %in% jkr)) {
          crosscheck_eag_nums(rep = rep1)
        } else enp
        nums3.1 <- sample(enp2$eag_num, 1, prob = enp2$prob)
        enp2 <- if (!all(a2 %in% jkr)) {
          crosscheck_eag_nums(rep = rep2)
        } else enp
        nums3.2 <- sample(enp2$eag_num, 1, prob = enp2$prob)
      }
    }

    nums <- eag_numbers(
      e2 = c(nums2.1, nums2.2),
      e3 = c(nums3.1, nums3.2)
    )
    return(nums)
  }), class = c("remapper", "function"))
}

#' @export
print.remapper <- function(x, ...) {
  cat("function <remapper> (a)", sep = "")
}

