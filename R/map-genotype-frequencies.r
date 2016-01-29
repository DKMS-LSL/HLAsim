#' Map observed/expected genotype frequencies to sampled genotypes
#'
#' @param x Sample returned by \code{\link{inject_errors}}.
#'
#' @return A \code{gtf_tbl} object that can be used as a training/test
#' data set with the slots:
#' \itemize{
#'   \item "Class": Has an error been injected in the EAGs of a sample (\dQuote{event},
#'          \emph{E}), or not, (\dQuote{nonevent}, \emph{N}).
#'   \item "genotype": The genotype in format: \dQuote{03:FYKD/04:ADCGE}.
#'   \item "zygosity": One of \sQuote{heterozygous} or \sQuote{homozygous}.
#'   \item "conc": Randomly attributed DNA concentrations.
#'   \item "pexp": The expected genotype frequency.
#'   \item "pobs": The observed genotype frequency.
#'   \item "log_pexp": The log of the expected genotype frequency.
#'   \item "log_pobs": The log of the observed genotype frequency.
#'   \item "log_fold_diff": The log-fold difference between observed
#'          and expected genotype frequency.
#' }
#' and the attributes:
#' \itemize{
#'  \item "penv": The enclosing environment of the <geno_table> constructor containing
#'  memoised mapping and remapping functions created by \code{\link{make_mapper}}, as
#'  well as genotype frequencies calculated by \code{\link{gtf}}.
#'  \item "breaks": The cut points used to slice up the DNA concentration.
#' }
#' @family simulation functions
#' @export
#' @examples
#' \dontrun{
#' ## Extract HLA-DPB1 genotype frequencies
#' dpb1 <- HLA("DPB1", "01/01/2014", "23/03/2015")
#'
#' ## Restrict the data to the German sample
#' dpb1.de <- dpb1[provenance == "DE"]
#'
#' ## Generate an EAG table
#' dpb1_eag1412 <- eag_table(gene = "DPB1", nextype_basis_id = "1412")
#'
#' ## Generate a distribution of DNA concentrations
#' conc <- sample_dna_concentration(dpb1.de, n = 1000, ncores = 8)
#'
#' ## Generate a sampling function
#' sample_dpb1_de <- make_genotype_sampler(dpb1.de, dpb1_eag1412)
#'
#' ## Load precomputed DNA-concentration-dependent error distribution for HLA-A
#' data(cdfA)
#'
#' ## Sample genotypes
#' n <- 10000
#' bin_size <- 3
#' perr <- 0.01
#' odds <- 0.25
#' ans <- sample_dpb1_de(conc, n, bin_size) %>%
#'   inject_errors(cdfA, perr, odds) %>%
#'   map_genotype_frequencies()
#' ans
#' }
map_genotype_frequencies <- function(x) {
  assertive::assert_is_all_of(x, "geno_table")
  ## get genotype frequencies based on remapped
  ## allele distribution
  g <- get("gtfrq", envir = attr(x, "penv"))
  ## merge the sample with the injected errors
  ## remove samples that did not remap properly
  s <- data.table::copy(samples(x))
  s[, Class := "N"]
  s[, bins := NULL]
  e <- data.table::copy(errors(x))
  e[, conc := s[e$idx, conc]]
  e[, Class := "E"]
  s <- s[!e$idx]
  se <- rbind(s, e)
  setkeyv(se, "idx")
  se[, idx := NULL]
  setkeyv(se, "genotype")
  rs <- g[se][, .(Class, genotype, eag_status, zygosity, conc, pexp, pobs, log_pexp, log_pobs, log_fold_diff)]
  rs <- na.omit(rs)
  rs[, genotype := factor(genotype)]
  rs[, eag_status := factor(eag_status)]
  rs[, zygosity := factor(zygosity)]
  rs[, Class := factor(Class)]
  setkeyv(rs, "genotype")
  rs <- dplyr::tbl_dt(rs)
  structure(rs, class = c("gtf_tbl", class(rs)))
}

