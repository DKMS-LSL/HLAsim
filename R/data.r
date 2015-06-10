#' EAG tables
#'
#' Datasets containing the EAG numbers and the corresponding alleles
#' for DRB1, DQB1, and DPB1 based on NeXtype Basis ID 1412 for Exons
#' 2 and 3
#'
#' @format Data tables with 4 columns:
#' \describe{
#'  \item{eag_num}{Unique numeric code for an Exon Allele Group.}
#'  \item{eag_allele}{HLA allele code}
#'  \item{exon}{Exon "2" or "3"}
#'  \item{allele_id}{Allele identifier}
#' }
#' @name eag_tables
#' @examples
#'   drb1_eag1412
"drb1_eag1412"

#' @rdname eag_tables
"dqb1_eag1412"

#' @rdname eag_tables
"dpb1_eag1412"

#' cdfA and cdfB
#'
#' Empirical cumulative distribution functions \code{\link{ecdf}} describing
#' DNA-concentration-dependent error rates based on data for \emph{HLA-A} and
#' \emph{HLA-B} extracted for the period between 01/09/2014 and 23/05/2015.
#'
#' @examples
#' data(cdfA)
#' plot(cdfA)
#'
#' ## extract the number of observations underlying an ecdf
#' get("nobs", environment(cdfA))
"cdfA"

#' @rdname cdfA
"cdfB"

#' Genotyping results for class II HLA genes.
#'
#' Genotyping results for the class II HLA genes, DRB1, DQB1, and DPB1,
#' processed between 01/01/2014 and 23/03/2015 at DKMS Life Science Lab.
#'
#' @format
#' \itemize{
#'   \item \code{drb1}: A <\code{\link{HLA}}> instance with 926,228 rows.
#'   \item \code{dqb1}: A <\code{\link{HLA}}> instance with 922,495 rows.
#'   \item \code{dpb1}: A <\code{\link{HLA}}> instance with 900,913 rows.
#' }
#'
#' @name geno_data
#' @examples
#' drb1
#' drb1$get_table()
#' drb1.de <- drb1[provenance == "DE"]
#' drb1.de$allele_frequency()
#' drb1.de$genotype_frequency()
"drb1"

#' @rdname geno_data
"dqb1"

#' @rdname geno_data
"dqb1"

#' Template DNA concentrations
#'
#' DNA concentration measures for 24150, 24304, and 25049 random samples
#' from Germany, Poland, and the UK, respectively.
#'
#' @format Data tables with 4 columns:
#' \describe{
#'  \item{lims_donor_id}{Sample identifier.}
#'  \item{ymd}{Date of measurement}
#'  \item{conc}{DNA concentration in ng/µl}
#'  \item{provenance}{Origin of sample. One of "DE", "PL", or "UK"}
#' }
#'
#' @examples
#'   concentration
"concentration"

