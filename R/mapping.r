## Helpers for mapping NMDP and G codes ####

merge_subtypes <- function(x, collapse = "/") {
  paste0(hla_sort(x), collapse = collapse)
}

merge_nmdp_subtypes <- function(x) {
  ## c("04:01:01:01", "04:01:01:02", "350:01")
  ## --> c("04:01", "350:01")
  f_ <- function(x) hla_collapse(x[1:2])
  merge_subtypes(unique(vapply(hla_burst(hla_sort(x)), f_, "")))
}

burst_subtypes <- function(x) {
  strsplit(x, split = "/", fixed = TRUE)[[1]]
}

#' @keywords internal
expand_nmdp_lookup <- function(exon, ntbl) {
  detect  <- stringi::stri_detect_regex
  `%do%`  <- foreach::`%do%`
  subtypes <- strsplit(ntbl$subtype, split = "/", fixed = TRUE)
  iPTN <- iterators::iter(vapply(subtypes, function(x) {
    paste0("^(", paste0(x, collapse = "|"), ")")
  }, FUN.VALUE = ""))
  iNM <- iterators::iter(ntbl$code)
  # ptn <- nextElem(iPTN)
  # nm <- nextElem(iNM)
  rs <- foreach::foreach(ptn = iPTN, nm = iNM, .combine = "rbind") %do% {
    tbl <- exon[detect(eag_allele, ptn)]
    styp <- if ((idx <- anyDuplicated(tbl$eag_num)) != 0) {
      tbl[eag_num == tbl[idx, eag_num], merge_subtypes(eag_allele)]
    } else {
      tbl[, merge_subtypes(eag_allele)]
    }
    data.table(code = nm, subtype = styp)
  }
  structure(rs, class = c("nmdp_table", "data.table", "data.frame"))
}

#' Encode subtypes to G or NMDP codes
#'
#' @details
#' Expect gtbl and ntbl2 (with 2-field subtypes)
#' in the parent environment.
#'
#' @param x A vector of alleles
#' @keywords internal
encode_ambiguities <- function(x, gtbl, ntbl2) {
  if (all(is.na(x))) {
    return(NA_character_)
  }
  if (length(ans <- gtbl[subtype == merge_subtypes(x), code]) == 0) {
    if (length(ans <- strip_field_four(x)) > 1) {
      styp <- merge_nmdp_subtypes(x)
      if (length(code <- ntbl2[subtype == styp, code]) == 0) {
        ## if this also fails, return the first allele
        return(x[1L])
      }
      ans <- paste0(field1(styp), ":", code)
    }
  }
  ans
}

#' Disambiguate NMDP and G codes
#'
#' @details
#' Expect gtbl and ntbl4 (4-field subtypes)
#' in the parent environment
#'
#' @param x An allele
#' @keywords internal
disambiguate_codes <- function(x, gtbl, ntbl4) {
  if (is_g_code(x)) {
    styp <- gtbl[code == x, subtype] %||% NA_character_
    burst_subtypes(styp)
  } else if (is_nmdp(x)) {
    styp <- ntbl4[code == field2(x), subtype] %||% NA_character_
    burst_subtypes(styp)
  } else x
}

remap_g_code <- function(x, gtbl) {
  styp <- gtbl[code == x, subtype] %||% NA_character_
  burst_subtypes(styp)
}

remap_nmdp_code <- function(x, ntbl4) {
  ## x <- "04:ADCGE"
  styp <- ntbl4[code == field2(x), subtype] %||% NA_character_
  burst_subtypes(styp)
}

#' Expand ambiguity codes
#'
#' Input: a vector of alleles with ambiguity codes (NMDP and G)
#' Output: a vector of alleles with ambiguity codes expanded.
#'
#' @param alleles
#' @param lookup
#' @keywords internal
expand_ambiguity_codes <- function(alleles, lookup) {
  `%do%`  <- foreach::`%do%`
  rs <- if (any(ambi <- is_ambiguous(alleles))) {
    iAMBI <- iterators::iter(alleles[ambi])
    foreach::foreach(a = iAMBI, .combine = "union") %do% {
      disambiguate_codes(a, gtbl = lookup$gtbl, ntbl4 = lookup$ntbl4)
    }
  } else character(0)
  hla_sort(union(alleles[!ambi], rs))
}


#' Make a function that disambiguates report alleles
#'
#' @param lookup
#' @return A function that takes a <rep_allele> as input
#' and returns an <eag_allele>
#' @keywords internal
make_disambiguator <- function(lookup) {
  ntbl <- lookup$nmdp_lookup
  gtbl <- lookup$g_lookup
  function(a) {
    stopifnot(is(a, "eag_allele"))
    a1 <- disambiguate_allele(allele(a, 1))
    a2 <- disambiguate_allele(allele(a, 2))
    eag_allele(a1, a2)
  }
}

#' Make a function that encodes eag alleles
#'
#' @param lookup
#' @return A function that takes an <eag_allele> as input
#' and returns a <rep_allele>
#' @keywords internal
make_encoder <- function(lookup) {
  ntbl <- lookup$nmdp_lookup
  gtbl <- lookup$g_lookup
  function(a) {
    stopifnot(is(a, "rep_allele"))
    a1 <- encode_ambiguities(allele(a, 1))
    a2 <- encode_ambiguities(allele(a, 2))
    rep_allele(a1, a2)
  }
}

## To be remmoved soon ####

#' Map eag numbers to an eag allele
#'
#' @param nums EAG numbers for exons 2 and 3 as
#' an <\code{\link{eag_numbers}}> object.
#' @param lookup A <\code{\link{lookup_list}}> object.
#' @return An <\code{\link{eag_allele}}>.
#' @keywords internal
#' @examples
#' nums <- eag_numbers(e2 = c(6673401, 6673404), e3 = c(6672138, 6670791))
#' lookup <- lookup_list(rep_dpb1, eag)
#'
#' map_eag_allele(nums, lookup)
#'
map_eag_allele <- function(nums, lookup) {
  stopifnot(
    is(nums, "eag_numbers"),
    is(lookup, "lookup_list")
  )
  a <- nums[1]; b <- nums[2]
  c <- nums[3]; d <- nums[4]
  ##  |  ex2  |  ex3  |
  ##  |-------|-------|
  ##  |   a   |   c   |
  ##  |-------|-------|
  ##  |   b   |   d   |
  ##  |-------|-------|
  ex2 <- exon(lookup, 2)
  ex3 <- exon(lookup, 3)
  jkr <- lookup[["joker_lookup"]][["eag_allele"]]

  ## No EAG nums for exon 2
  if (is.na(a) && is.na(b)) {
    ## no valid report possible
    return(eag_allele())
  }

  ## No EAG nums for exon 3
  if (is.na(c) && is.na(d)) {
    a_ <- unique(ex2[eag_num == a, eag_allele])
    b_ <- unique(ex2[eag_num == b, eag_allele])
    return(eag_allele(a_, b_))
  }

  ## One EAG on Exon 2; one EAG on Exon 3
  if (a == b && c == d) {
    a_ <- unique(ex2[eag_num == a, eag_allele])
    c_ <- c(unique(ex3[eag_num == c, eag_allele]), jkr)
    ac <- unique(intersect(a_, c_))
    if (none_empty(ac)) {
      return(eag_allele(ac, ac))
    }

    ## no valid report possible
    return(eag_allele())
  }

  ## One EAG on Exon 2; two EAGs on Exon 3
  if (a == b) {
    a_ <- unique(ex2[eag_num == a, eag_allele])
    c_ <- c(unique(ex3[eag_num == c, eag_allele]), jkr)
    d_ <- c(unique(ex3[eag_num == d, eag_allele]), jkr)
    ## match horizontally
    ac <- unique(intersect(a_, c_))
    ## match diagonally
    ad <- unique(intersect(a_, d_))
    if (none_empty(ac, ad)) {
      return(eag_allele(ac, ad))
    }

    ## no valid report possible
    return(eag_allele())
  }

  ## Two EAGs on Exon 2; one EAG on Exon 3
  if (c == d) {
    a_ <- unique(ex2[eag_num == a, eag_allele])
    b_ <- unique(ex2[eag_num == b, eag_allele])
    c_ <- c(unique(ex3[eag_num == c, eag_allele]), jkr)
    ## match horizontally
    ac <- unique(intersect(a_, c_))
    ## match diagonally
    bc <- unique(intersect(b_, c_))
    if (none_empty(ac, bc)) {
      return(eag_allele(ac, bc))
    }

    ## no valid report possible
    return(eag_allele())
  }

  ## Two EAGs on Exon 2; two EAGs on Exon 3
  a_ <- unique(ex2[eag_num == a, eag_allele])
  b_ <- unique(ex2[eag_num == b, eag_allele])
  c_ <- c(unique(ex3[eag_num == c, eag_allele]), jkr)
  d_ <- c(unique(ex3[eag_num == d, eag_allele]), jkr)
  ## match horizontally
  ac <- unique(intersect(a_, c_))
  bd <- unique(intersect(b_, d_))
  ## match diagonally
  ad <- unique(intersect(a_, d_))
  bc <- unique(intersect(b_, c_))

  ## If a ∩ c != {} && a ∩ d != {} && b ∩ d != {} && b ∩ c != {}
  ## Allele2 1 = union(ac, ad); Allele 2 = union(bd, bc)
  if (none_empty(ac, ad, bd, bc)) {
    return(eag_allele(union(ac, ad), union(bd, bc)))
  }
  ## If a ∩ c != {} && b ∩ d != {}
  ## Allele 1 = ac, Allele 2 = bd
  if (none_empty(ac, bd)) {
    return(eag_allele(ac, bd))
  }
  ## If a ∩ d != {} && b ∩ c != {}
  ## Allele 1 = ad, Allele 2 = bc
  if (none_empty(ad, bc)) {
    return(eag_allele(ad, bc))
  }

  ## no valid report possible
  return(eag_allele())
}

#' Map an eag allele to a report allele
#'
#' @param a An <\code{\link{eag_allele}}>.
#' @param lookup A <\code{\link{lookup_list}}>
#' @return A <\code{\link{rep_allele}}>.
#' @keywords internal
#' @examples
#' nums <- eag_numbers(e2 = c(6673401, 6673404), e3 = c(6672138, 6670791))
#' lookup <- lookup_list(rep_dpb1, eag)
#'
#' a <- map_eag_allele(nums, lookup)
#' rp <- map_rep_allele(a, lookup)
#'
#' rp1 <- rep_allele(a1 = "04:ADCGE", a2 = "13:01:01G")
#' rp == rp1
#'
map_rep_allele <- function(a, lookup) {
  stopifnot(is(a, "eag_allele"))
  a1_rep <- encode_ambiguities(allele(a, 1), lookup$gtbl, lookup$ntbl2)
  a2_rep <- encode_ambiguities(allele(a, 2), lookup$gtbl, lookup$ntbl2)
  rep_allele(a1_rep, a2_rep)
}
