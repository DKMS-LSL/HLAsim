#' Class: nmdp_tbl
#'
#' Constructor for an <\code{nmdp_tbl}> object.
#'
#' @source \url{https://bioinformatics.bethematchclinical.org/HLA-Resources/Allele-Codes/Allele-Code-Lists/Allele-Code-List-in-Numerical-Order/}.
#'
#' @return
#' A table of mappings between NMDP codes and allelic subtypes with the fields:
#' \itemize{
#'  \item "code": An NMDP code.
#'  \item "subtype": The '/'-separated subtypes into which an NMDP code expands.
#' }
#' @seealso \code{\link{g_table}}
#' @export
#' @examples
#' \dontrun{
#' nmdp_codes <- nmdp_table()
#' }
nmdp_table <- function() {
  u <- "https://bioinformatics.bethematchclinical.org/HLA/numeric.v3.zip"
  tmp <- tempfile(fileext = ".zip")
  on.exit(unlink(tmp))
  curl::curl_download(u, tmp)
  con <- unz(tmp, filename = "numer.v3.txt")
  rs <- dtplyr::tbl_dt(
    data.table::setDT(
      scan(con, what = list("character", "character"),
           nlines = -1, sep = "\t", skip = 3, quiet = TRUE)
    )
  )
  close(con)
  data.table::setnames(rs, names(rs), c("code", "subtype"))
  data.table::setkeyv(rs, "code")
  structure(rs, class = c("nmdp_tbl", class(rs)))
}

#' @export
print.nmdp_tbl <- function(x, ..., n = 5) {
  cat("NMDP table: ", sep = "")
  NextMethod(n = n)
  cat("...\n", sep = "")
}

#' Class: g_tbl
#'
#' Constructor for a <\code{g_tbl}> object.
#'
#' @source \url{http://hla.alleles.org/nomenclature/g_groups.html}.
#'
#' @return
#' A table with the fields:
#' \itemize{
#'  \item "gene": A HLA gene.
#'  \item "code": The G Code.
#'  \item "subtype": The '/'-separated subtypes into which a G Code expands.
#' }
#' @seealso \code{\link{nmdp_table}}
#' @export
#' @examples
#' \dontrun{
#' g_codes <- g_table()
#' }
g_table <- function() {
  con <- curl::curl("http://hla.alleles.org/wmda/hla_nom_g.txt")
  on.exit(close(con))
  tryCatch(open(con), error = function(e) {
    stop("Trying to access http://hla.alleles.org/wmda/: ", e$message, call. = FALSE)
  })

  if (readLines(con, n = 1) != "# file: hla_nom_g.txt") {
    warning("Possibly malformed file \"hla_nom_g.txt\" ",
            "downloaded from http://hla.alleles.org/wmda/.", immediate. = TRUE)
  }
  rs <- read.csv(con, header = FALSE, colClasses = "character", sep = ";", comment.char = "#")
  rs <- dtplyr::tbl_dt(data.table::setDT(rs))
  data.table::setnames(rs, names(rs), c("gene", "subtype", "code"))
  rs <- rs[, gene := paste0("HLA-", sub("*", "", gene, fixed = TRUE))]
  rs <- rs[nzchar(rs[, code])][, .(gene, code, subtype)]
  data.table::setkeyv(rs, "gene")
  structure(rs, class = c("g_tbl", class(rs)))
}

#' @export
print.g_tbl <- function(x, ..., n = 5) {
  cat("G-Codes: ", sep = "")
  NextMethod(n = n)
  cat("...\n", sep = "")
}

#' Class: eag_tbl
#'
#' @description
#' Constructor for an <\code{eag_tbl}> object.
#'
#' THIS FUNCTION REQURIES ACCESS TO INTERNAL DATABASES AT DKMS LSL.
#'
#' Use the prepackaged datasets (\code{\link{eag_tables}}) instead.
#' @param gene An HLA gene. One of \sQuote{A}, \sQuote{B}, \sQuote{C}, \sQuote{DRB1},
#' \sQuote{DQB1}, \sQuote{DPB1}.
#' @param nextype_basis_id "NeXtype Basis ID".
#'
#' @return
#' A table with the fields:
#' \itemize{
#'  \item "eag_num": Numeric code describing an Exon Allele Group.
#'  \item "eag_allele": HLA allele code.
#'  \item "exon": Exon '2' or '3'.
#'  \item "allele_id": Allele identifier.
#' }
#' and the attributes:
#' \itemize{
#'  \item "nextype_basis_id"
#'  \item "gene"
#' }
#' @export
#' @examples
#' \dontrun{
#' dpb1_eag1412 <- eag_table("DPB1", "1412")
#' }
eag_table <- function(gene = "DPB1", nextype_basis_id = "1412") {
  if (!requireNamespace("orcl", quietly = TRUE)) {
    stop("Need access to internal databases at DKMS LSL. Please use the packaged tables.",
         call. = FALSE)
  }
  assertive::assert_is_scalar(nextype_basis_id, "length")
  gene <- match_hla_gene(gene)
  fmt <- "
  SELECT
    a.eag_num AS eag_num, a.nmdp_new AS eag_allele,
    TO_NUMBER(b.exon) AS exon, a.allele_id, a.dna_version_id,
    a.allele_num
  FROM ngstest.nextype_alleles_per_eag a
  INNER JOIN ngstest.nextype_eags b
  ON a.eag_num           = b.eag_num
  WHERE a.gene           = '%s'
  AND a.nextype_basis_id = %s"
  rs <- dtplyr::tbl_dt(orcl::ora_query(sprintf(fmt, gene, nextype_basis_id), user = "ngstest"))
  rs[, allele_id := strsplitN(allele_id, ";", 1)]
  dna_version_id <- rs[, unique(dna_version_id)]
  rs[, dna_version_id := NULL]
  structure(
    rs,
    dna_version_id = dna_version_id,
    nextype_basis_id = nextype_basis_id,
    gene = gene,
    class = c("eag_tbl", class(rs))
  )
}

#' @export
print.eag_tbl <- function(x, ..., n = 5) {
  fmt <- "EAG table [Gene: %s, BasisID: %s, DNAversion: %s]\n"
  cat(sprintf(fmt, gene(x), nextype_basis_id(x), dna_version_id(x)), sep = "")
  NextMethod(n = n)
}

#' @export
exon.eag_tbl <- function(x, n) {
  n <- match.arg(as.character(n), c("2", "3"))
  structure(
    x[exon == n],
    nextype_basis_id = nextype_basis_id(x),
    gene = gene(x),
    class = class(x)
  )
}

#' @export
nextype_basis_id.eag_tbl <- function(x) {
  attr(x, "nextype_basis_id")
}

#' @export
dna_version_id.eag_tbl <- function(x) {
  attr(x, "dna_version_id")
}

#' @export
gene.eag_tbl <- function(x) {
  attr(x, "gene")
}

#' Class: gtf_tbl
#'
#' Constructor for <\code{gtf_tbl}> objects. Calculates observed
#' and expected genotype frequencies.
#'
#' @param x A <\code{\link{HLA}}> object.
#' @return A table with the fields:
#' \itemize{
#'  \item "genotype": \strong{Key}. The genotype in the format "A1/A2".
#'  \item "pexp": Expected genotype frequencies.
#'  \item "pobs": Observed genotype frequencies.
#'  \item "log_pexp": Log of expected genotype frequenciies.
#'  \item "log_pobs": Log of expected genotype frequenciies.
#'  \item "nexp": Expected number of samples for a genotype.
#'  \item "nobs": Observed number of samples for a genotype.
#'  \item "chisq": \eqn{\chi^2}-value describing the difference
#'    between expected and observed genotype frequency.
#'  \item "log_fold_diff": Log-fold diffence between expected
#'    and observed genotype frequency.
#' }
#'
#' @keywords internal
#' @export
#' @examples
#' \dontrun{
#' ## Extract HLA-DPB1 genotype frequencies
#' dpb1 <- HLA("DPB1", "01/01/2014", "23/03/2015")
#'
#' ## Restrict the data to the German sample
#' dpb1.de <- dpb1[provenance == "DE"]
#'
#' ## Remove low-resolution four-digit codes
#' dpb1.de <- dpb1.de[allele1 %ni% dpb1_four_digit_codes | allele2 %ni% dpb1_four_digit_codes]
#'
#' ## Remove cases of unknown alleles
#' dpb1.de <- dpb1.de[field2(allele1) != "XXX" | field2(allele2) != "XXX"]
#' dpb1.de <- dpb1.de[toupper(allele1) != "NEW" | toupper(allele2) != "NEW"]
#'
#' ## Calculate observed and expected genotype frequencies
#' gtf_table(dpb1.de)
#' }
gtf_table.HLA <- function(x) {
  ## purge existing genotype/allele frequencies
  x$purge_frequencies()
  ## calculate observed genotype frequencies
  gtf_obs <- dplyr::select(x$genotype_frequency(), genotype, pobs = frequency)
  ## calculate expected genotype frequencies
  af <- x$allele_frequency()
  gtf_exp <- expected_genotype_frequencies(alleles = af$allele, frequencies = af$frequency)
  ## Which theoretically possible genotypes are missing from the sample
  missing <- gtf_exp[!as.character(gtf_exp$genotype) %chin% as.character(gtf_obs$genotype)]
  ## Setting P_obs zero and merge with the table of observed frequencies.
  missing <- missing[, `:=`(pexp = NULL, pobs = 0)]
  gtf_obs <- rbind(gtf_obs, missing)
  gtf <- gtf_exp[gtf_obs][order(-pexp)]
  gtf <- gtf[, `:=`(log_pexp = base::log(pexp), log_pobs = base::log(pobs))]
  ## replace -Inf with -30 as lower bound at the log scale
  gtf <- gtf[, log_pobs := ifelse(log_pobs == -Inf, -30, log_pobs)]
  ## transform frequencies to counts.
  N <- af[, sum(count)] / 2 # sample size
  gtf <- gtf[, `:=`(nexp = pexp*N, nobs = pobs*N )]
  # calculate Chi^2 and log fold difference as measures of deviation between
  # expected and observed genotype frequencies
  gtf <- gtf[, chisq := (nobs - nexp)^2 / nexp]
  gtf <- gtf[, log_fold_diff := log_pobs - log_pexp]
  setkeyv(gtf, "genotype")
  if ("eag_status" %in% names(tbl <- x$get_table())) {
    tbl <- tbl[, .(eag_status = unique(eag_status)), by = genotype]
    setkeyv(tbl, "genotype")
    gtf <- tbl[gtf]
  }
  gtf <- dtplyr::tbl_dt(gtf)
  structure(gtf, class = c("gtf_tbl", class(gtf)))
}

#' @export
print.gtf_tbl <- function(x, n = 4, ...) {
  cat("GTF table [Nonevents: ", sum(x$Class == "N"), "; Events: ", sum(x$Class == "E"), "]\n", sep = "")
  NextMethod(n = n)
}

joker_table <- function(eag) {
  gene           <- gene(eag)
  dna_version_id <- dna_version_id(eag)
  exon           <- 3
  if (gene %in% c("HLA-DQB1", "HLA-DRB1", "HLA-DPB1") &&
      dna_version_id == "52" &&
      !requireNamespace("orcl", quietly = TRUE)) {
    switch(gene,
           `HLA-DQB1` = dqb1_jokers1412,
           `HLA-DRB1` = drb1_jokers1412,
           `HLA-DPB1` = dpb1_jokers1412)
  } else {
    ## TODO: Revert table names once migration is finished
    fmt <- "
    SELECT joker_num, allele_num, nmdp_new AS eag_allele
    FROM ngsrep.nextype_jokers
    WHERE gene         = '%s'
    AND dna_version_id = %s
    AND exon           = %s"
    rs <- orcl::ora_query(sprintf(fmt, gene, dna_version_id, exon))
    setkeyv(rs, "eag_allele")
  }
}

partials_table <- function(eag) {
  gene             <- gene(eag)
  dna_version_id   <- dna_version_id(eag)
  nextype_basis_id <- nextype_basis_id(eag)
  exon             <- 3
  if (gene %in% c("HLA-DQB1", "HLA-DRB1", "HLA-DPB1") &&
      dna_version_id == "52" &&
      !requireNamespace("orcl", quietly = TRUE)) {
    switch(gene,
           `HLA-DQB1` = dqb1_partials1412,
           `HLA-DRB1` = drb1_partials1412,
           `HLA-DPB1` = dpb1_partials1412)
  } else {
    ## TODO: Revert table names once migration is finished
    fmt <- "
    SELECT a.allele_num,
    b.nmdp_new AS eag_allele
    FROM ngsrep.nextype_partials a
    INNER JOIN ngsrep.nextype_alleles_per_eag b
    ON a.allele_num_partial = b.allele_num
    AND a.dna_version_id    = b.dna_version_id
    AND a.gene              = b.gene
    WHERE a.gene            = '%s'
    AND b.dna_version_id    = '%s'
    AND a.exon              = %s
    AND b.nextype_basis_id  = %s"
    rs <- orcl::ora_query(sprintf(fmt, gene, dna_version_id, exon, nextype_basis_id))
    setkeyv(rs, "allele_num")
  }
}


