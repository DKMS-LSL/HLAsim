#' Generate a mapping function
#'
#' Make mappers for mapping \code{EAG numbers -> Report Alleles} (\code{make_mapper}) and
#' for mapping \code{Report Alleles -> EAG numbers} (\code{make_remapper}).
#'
#' @param lookup A \code{\link{lookup_list}} object.
#' @return A function
#'
#' @export
#' @examples
#'
#' ## load the prepackaged date for DPB1 report alleles and
#' ## the DPB1 EAG table
#' data(rep_dpb1, package = "HLAdata"); data(eag, package = "HLAdata")
#'
#' lookup <- lookup_list(rep_dpb1, eag)
#' map <- make_mapper(lookup)
#'
#' nums <- eag_numbers(e2 = c("6680148", "6673401"), e3 = "6670791")
#' map(nums)
#'
#' remap <- make_remapper(lookup)
#' new_nums <- remap(map(nums))
#'
#' nums == new_nums
make_mapper <- function(lookup) {
  stopifnot(is(lookup, "lookup_list"))
  memoise <- memoise::memoise

  ## global variables
  ex2   <- exon(lookup, 2)
  ex3   <- exon(lookup, 3)
  ntbl2 <- lookup$ntbl2
  ntbl4 <- lookup$ntbl4
  gtbl  <- lookup$gtbl
  jkr   <- lookup$jkr_tbl$eag_allele
  prt   <- lookup$prt_tbl

  structure(memoise(function(nums) {
    stopifnot(is(nums, "eag_numbers"))
    a <- nums[1]; b <- nums[2]
    c <- nums[3]; d <- nums[4]

    ##  |  ex2  |  ex3  |
    ##  |-------|-------|
    ##  |   a   |   c   |
    ##  |-------|-------|
    ##  |   b   |   d   |
    ##  |-------|-------|

    ## Fail if there are no EAG nums for exon 2
    rs <- if (is.na(a) && is.na(b)) {
      eag_allele()
    }
    ## There are no EAG nums for exon 3
    else if (is.na(c) && is.na(d)) {
      a_ <- ex2[eag_num == a, eag_allele]
      b_ <- ex2[eag_num == b, eag_allele]
      eag_allele(a_, b_)
    }
    ## One EAG on Exon 2; one EAG on Exon 3
    else if (a == b && c == d) {
      a_ <- ex2[eag_num == a, eag_allele]
      c_ <- ex3[eag_num == c, eag_allele]
      c_partials <- if (NROW(prt) > 0) {
        prt[ex3[eag_num == c, allele_num][1], eag_allele]
      } else NULL
      c_ <- c(c_, c_partials, jkr)

      ac <- unique(intersect(a_, c_))
      if (none_empty(ac)) {
        eag_allele(ac, ac)
      } else {
        eag_allele()
      }
    }
    ## One EAG on Exon 2; two EAGs on Exon 3
    else if (a == b) {
      a_ <- ex2[eag_num == a, eag_allele]
      c_ <- ex3[eag_num == c, eag_allele]
      d_ <- ex3[eag_num == d, eag_allele]
      if (NROW(prt) > 0) {
        c_partials <- prt[ex3[eag_num == c, allele_num][1], eag_allele]
        d_partials <- prt[ex3[eag_num == d, allele_num][1], eag_allele]
      } else {
        c_partials <- d_partials <- NULL
      }
      c_ <- c(c_, c_partials, jkr)
      d_ <- c(d_, d_partials, jkr)
      ## match horizontally
      ac <- unique(intersect(a_, c_))
      ## match diagonally
      ad <- unique(intersect(a_, d_))
      if (none_empty(ac, ad)) {
        eag_allele(ac, ad)
      } else {
        eag_allele()
      }
    }
    ## Two EAGs on Exon 2; one EAG on Exon 3
    else if (c == d) {
      a_ <- ex2[eag_num == a, eag_allele]
      b_ <- ex2[eag_num == b, eag_allele]
      c_ <- ex3[eag_num == c, eag_allele]
      c_partials <- if (NROW(prt) > 0) {
        prt[ex3[eag_num == c, allele_num][1], eag_allele]
      } else NULL
      c_ <- c(c_, c_partials, jkr)
      ## match horizontally
      ac <- unique(intersect(a_, c_))
      ## match diagonally
      bc <- unique(intersect(b_, c_))
      if (none_empty(ac, bc)) {
        eag_allele(ac, bc)
      } else {
        eag_allele()
      }
    }
    ## Two EAGs on Exon 2; two EAGs on Exon 3
    else {
      a_ <- ex2[eag_num == a, eag_allele]
      b_ <- ex2[eag_num == b, eag_allele]
      c_ <- ex3[eag_num == c, eag_allele]
      d_ <- ex3[eag_num == d, eag_allele]
      if (NROW(prt) > 0) {
        c_partials <- prt[ex3[eag_num == c, allele_num][1], eag_allele]
        d_partials <- prt[ex3[eag_num == d, allele_num][1], eag_allele]
      } else {
        c_partials <- d_partials <- NULL
      }
      c_ <- c(c_, c_partials, jkr)
      d_ <- c(d_, d_partials, jkr)
      ## match horizontally
      ac <- unique(intersect(a_, c_))
      bd <- unique(intersect(b_, d_))
      ## match diagonally
      ad <- unique(intersect(a_, d_))
      bc <- unique(intersect(b_, c_))
      ## Full reciprocal match → exon shuffling
      ## a ∩ c ≠ {} ∧ a ∩ d ≠ {} ∧ b ∩ d ≠ {} ∧ b ∩ c ≠ {} ⇒
      ## allele 1 = ac ∪ ad, allele 2 = bd ∪ bc
      if (none_empty(ac, ad, bd, bc)) {
        eag_allele(union(ac, ad), union(bd, bc))
      }
      ## Full horizontal match
      ## a ∩ c ≠ {} ∧ b ∩ d ≠ {} ⇒ allele 1 = ac, allele 2 = bd
      else if (none_empty(ac, bd)) {
        eag_allele(ac, bd)
      }
      ## Full diagonal match
      ## a ∩ d ≠ {} ∧ b ∩ c ≠ {} ⇒ allele 1 = ad, bllele 2 = bc
      else if (none_empty(ad, bc)) {
        eag_allele(ad, bc)
      }
      ## Semi-reciprocal match ⇒ no valid report alleles possible
      else {
        eag_allele()
      }
    }
    if (is.na(rs)) {
      return(rep_allele())
    }
    a1_rep <- encode_ambiguities(x = allele(rs, 1), gtbl, ntbl2)
    a2_rep <- encode_ambiguities(x = allele(rs, 2), gtbl, ntbl2)
    if (a1_rep < a2_rep) {
      rep_allele(a1_rep, a2_rep)
    } else {
      rep_allele(a2_rep, a1_rep)
    }
  }), class = c("mapper", "function"))
}

#' @export
print.mapper <- function(x, ...) {
  cat("function <mapper> (nums)", sep = "")
}
