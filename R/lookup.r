#' Class: lookup_list
#'
#' Constructor for a <lookup_list> object.
#'
#' Inherits from <list>.
#'
#' @param alleles A character vector of alleles.
#' @param eag EAG alleles <\code{\link{eag_tbl}}>.
#'
#' @return
#' A list with the elements:
#' \describe{
#'  \item{e2_tbl}{A mapping of rep_allele <-> eag_alllele <-> eag_num for exon 2}
#'  \item{e3_tbl}{A mapping of rep_allele <-> eag_alllele <-> eag_num for exon 2}
#'  \item{jkr_tbl}{A table of exon 3 joker alleles}
#'  \item{prt_tbl}{A mapping table with partials}
#'  \item{ntbl2}{2-field NMDP-code lookup}
#'  \item{ntbl4}{4-field NMDP-code lookup}
#'  \item{gtbl}{G-code lookup}
#' }
#' and the attributes:
#' \itemize{
#'  \item{nextype_basis_id}
#'  \item{gene}
#' }
#' @export
#' @examples
#' \dontrun{
#' library("HLA")
#'
#' ## get a vector of reported DPB1 alleles
#' dpb1 <- HLA::HLA("DPB1")
#' alleles <- HLA::get_table(dpb1)[, .(allele1, allele2)]
#' alleles <- sort(unique(c(alleles[, allele1], alleles[, allele2])))
#' alleles <- alleles[-grep("(:XXX|NEW)", alleles)]
#'
#' ## get an EAG table for DPB1 and a NeXtype Basis ID
#' eag <- eag_table("DPB1", "1267")
#'
#' lookup <- lookup_list(alleles, eag)
#' }
lookup_list <- function(alleles, eag) {
  assertthat::assert_that(
    all(is.hla(alleles)),
    is(eag, "eag_tbl")
  )
  stwf <- stringi::stri_startswith_fixed
  `%do%` <- foreach::`%do%`

  ## load prepackaged NMDP codes or get them from the internet
  if (!requireNamespace("HLAdata", quietly = TRUE) ||
      (rs <- data(nmdp_codes, package = "HLAdata", envir = environment())) != "nmdp_codes") {
    message("Fetching NMDP codes ...")
    nmdp_codes <- nmdp_table()
  }
  ## load prepackaged G codes or get them from the internet
  if (!requireNamespace("HLAdata", quietly = TRUE) ||
      (rs <- data(g_codes, package = "HLAdata", envir = environment())) != "g_codes") {
    message("Fetching G-Codes ...")
    g_codes <- g_table()
  }
  ## generate lookup table for relevant nmdp codes
  nmdp_lookup <- generate_nmdp_lookup(alleles, nmdp_codes)
  ## generate lookup table for relevant expanded nmdp codes
  exp_nmdp_lookup <- expand_nmdp_lookup(exon = exon(eag, 2), ntbl = nmdp_lookup)
  ## generate lookup table for relevant g codes
  g_lookup <- generate_g_lookup(alleles, eag, g_codes)

  ## The NMDP and G codes are kicked out of the allele list before
  ## constructing the exon-specific lookup table.
  alleles <- alleles[!is_ambiguous(alleles)]

  ## Alleles reported as G codes are expanded into their subtypes.
  ## We check for the occurrence of these alleles in the EAG table
  ## and if alleles do not occur under the present NeXtype Basis
  ## we update the g_lookup table accordingly.
  g_styp <- strip_field_four(unlist(lapply(g_lookup$subtype, burst_subtypes)))
  missing_g_alleles <- Filter(function(x) !any(stwf(eag$eag_allele, x)), g_styp)
  g_lookup <- update_g_lookup(g_lookup, missing_g_alleles)
  g_styp <- g_styp[g_styp %ni% missing_g_alleles]

  ## Alleles reported as NMDP codes are expanded into their subtypes
  all_eag_alleles <- eag$eag_allele
  nmdp_styp <- strip_field_four(unlist(lapply(exp_nmdp_lookup$subtype, burst_subtypes)))

  ## Reported alleles and subtypes are merged.
  alleles2 <- hla_sort(unique(c(alleles, nmdp_styp, g_styp)))
  allele_tbl <- data.table(rep_allele = alleles2 , rep_allele_f3 = alleles2)
  setkeyv(allele_tbl, "rep_allele_f3")

  structure(list(
    e2_tbl  = lookup_table(allele_tbl, eagn = exon(eag, 2)),
    e3_tbl  = lookup_table(allele_tbl, eagn = exon(eag, 3)),
    jkr_tbl = joker_table(eag),
    prt_tbl = partials_table(eag),
    ntbl2   = nmdp_lookup,
    ntbl4   = exp_nmdp_lookup,
    gtbl    = g_lookup
  ),
  nextype_basis_id = nextype_basis_id(eag),
  gene             = gene(eag),
  class            = c("lookup_list", "list"))
}

generate_nmdp_lookup <- function(alleles, nmdp_codes) {
  x <- data.table(f1 = field1(alleles), f2 = field2(alleles))
  x <- x[grepl("^[A-Z]", f2)]
  if (is.null(key(nmdp_codes)) || key(nmdp_codes) != "code") {
    setkeyv(nmdp_codes, "code")
  }
  setkeyv(x, "f2")
  x <- nmdp_codes[x]
  x[, subtype := unlist(purrr::map2(f1, subtype, function(a, b) {
    if (!grepl(":", b, fixed = TRUE)) {
      slash(paste0(a, ":", strsplit(b, split = "/", fixed = TRUE)[[1]]))
    } else b
  }))]
  x[, f1 := NULL]
  x
}

generate_g_lookup <- function(alleles, eag, g_codes) {
  g_lookup <- g_codes[gene == gene(eag), .(code, subtype)]
  setkeyv(g_lookup, "code")
  g_lookup
}

update_g_lookup <- function(g_lookup, missing_g_alleles) {
  lookup <- copy(g_lookup)
  styp <- strsplit(lookup$subtype, split = "/", fixed = TRUE)
  lookup[, subtype := unlist(Map(function(x) merge_subtypes(x[x %ni% missing_g_alleles]), styp))]
  lookup
}

#' @describeIn lookup_list
#' @keywords internal
#' @export
print.lookup_list <- function(x, ..., n = 2) {
  fmt1 <- "Lookup list [Gene: %s, BasisID: %s]\n"
  fmt2 <- "Lookup table [Exon: %s] [%s x %s]\n"
  cat(sprintf(fmt1, gene(x), nextype_basis_id(x)), sep = "")
  cat("\n")
  cat(sprintf(fmt2, "2", nrow(exon(x, 2)), ncol(exon(x, 2))), sep = "")
  print(head(exon(x, 2), n = n), ...)
  cat("\n")
  cat(sprintf(fmt2, "3", nrow(exon(x, 3)), ncol(exon(x, 3))), sep = "")
  print(head(exon(x, 3), n = n), ...)
  cat("\n")
  cat("Joker lookup: [", nrow(x$jkr_tbl), " x ", ncol(x$jkr_tbl) ,"]\n", sep = "")
  print(head(x$jkr_tbl, n = n), ...)
  cat("\n")
  cat("Partials lookup: [", nrow(x$prt_tbl), " x ", ncol(x$prt_tbl) ,"]\n", sep = "")
  print(head(x$prt_tbl, n = n), ...)
  cat("\n")
  cat("NMDP code lookup: [", nrow(x$ntbl2), " x ", ncol(x$ntbl2) ,"]\n", sep = "")
  print(head(x$ntbl2, n = n), ...)
  cat("\n")
  cat("Expanded NMDP code lookup: [", nrow(x$ntbl4), " x ", ncol(x$ntbl4) ,"]\n", sep = "")
  print(head(x$ntbl4, n = n), ...)
  cat("\n")
  cat("G code lookup: [", nrow(x$gtbl), " x ", ncol(x$gtbl) ,"]\n", sep = "")
  print(head(x$gtbl, n = n), ...)
}

#' @describeIn lookup_list
#' @keywords internal
#' @export
gene.lookup_list <- function(x) {
  attr(x, "gene")
}

#' @describeIn lookup_list
#' @keywords internal
#' @export
nextype_basis_id.lookup_list <- function(x) {
  attr(x, "nextype_basis_id")
}

#' @describeIn lookup_list
#' @keywords internal
#' @export
exon.lookup_list <- function(x, n) {
  n <- match.arg(as.character(n), c("2", "3"))
  if (n == "2") x$e2_tbl else x$e3_tbl
}

#' Class: lookup_table
#'
#' Constructor for a <\code{lookup_table}>.
#'
#' Inherits from <\code{data.table}>
#'
#' @param allele_tbl
#' @param eagn An \code{\link{eag_tbl}} for either Exon 2 or 3 .
#'
#' @return
#' A table with the fields:
#' \describe{
#'  \item{rep_allele}{The reported allele codes (including NMDP and G codes)}
#'  \item{eag_allele}{The fully resolved HLA allele codes}
#'  \item{eag_num}{Unique numeric code for an Exon Allele Group.}
#' }
#' and the attributes:
#' \itemize{
#'  \item{nextype_basis_id}
#'  \item{gene}
#'  \item{exon}
#' }
#'
#' @seealso \code{\link{lookup_list}}
#' @export
lookup_table <- function(allele_tbl, eagn) {
  stwf <- stringi::stri_startswith_fixed
  setkeyv(eagn, "eag_allele")
  ## cut off the fourth field (intronic differences)
  ## from the eag_alleles to generate three-field alleles
  eag_allele_full <- hla_sort(eagn$eag_allele)
  eag_allele_f3   <- strip_field_four(eag_allele_full, unique = FALSE)
  rp1 <- unique(data.table(eag_allele_f3, eag_allele_full, key = "eag_allele_f3,eag_allele_full"))
  setkeyv(rp1, "eag_allele_f3")
  rp2 <- rp1[allele_tbl, allow.cartesian = TRUE][, .(rep_allele, eag_allele_full)]
  setkeyv(rp2, "eag_allele_full")
  # eagn[rp2][, .(rep_allele, eag_allele, eag_num)]
  ## generate lookup table
  tbl <- na.omit(eagn[rp2][, .(rep_allele, eag_allele, eag_num)])
  eag_nums <- unique(tbl$eag_num)
  setkeyv(eagn, "eag_num")
  eagn2 <- eagn[eag_nums, .(eag_allele, eag_num, allele_num)]
  setkeyv(eagn2, c("eag_num", "eag_allele"))
  setkeyv(tbl, c("eag_num", "eag_allele"))
  tbl <- tbl[eagn2]
  tbl[, rep_allele := ifelse(is.na(rep_allele), strip_field_four(eag_allele), rep_allele)]
  lookup <- structure(
    tbl,
    nextype_basis_id = nextype_basis_id(eagn),
    gene             = gene(eagn),
    exon             = eagn[, unique(exon)],
    class            = c("lookup_table", "data.table", "data.frame")
  )
  setkeyv(lookup, c("rep_allele", "eag_allele"))
  lookup
}

#' @describeIn lookup_table
#' @keywords internal
#' @export
print.lookup_table <- function(x, ..., n = 3) {
  fmt <- "Lookup table [Exon: %s, Gene: %s, BasisID: %s] [%s x %s]\n"
  cat(sprintf(fmt, exon(x), gene(x), nextype_basis_id(x), nrow(x), ncol(x)), sep = "")
  NextMethod(topn = n)
}

#' @describeIn lookup_table
#' @keywords internal
#' @export
exon.lookup_table <- function(x) {
  attr(x, "exon")
}

#' @describeIn lookup_table
#' @keywords internal
#' @export
gene.lookup_table <- function(x) {
  attr(x, "gene")
}

#' @describeIn lookup_table
#' @keywords internal
#' @export
nextype_basis_id.lookup_table <- function(x) {
  attr(x, "nextype_basis_id")
}


