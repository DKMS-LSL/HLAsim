ichunk <- function(x, chunksize) {
  l <- length(x)
  i <- seq.int(from = 1L, to = l, by = chunksize)
  j <- c(i[-1] - 1, l)
  iterators::iter(
    lapply(.mapply(seq.int, list(i, j), NULL), function(k) x[k])
  )
}

comma <- function(...) paste0(..., collapse = ",")

wrap <- function(x, wrap = "'") sprintf('%s%s%s', wrap, x, wrap)

#' Extract DNA concentration for LIMS_DONOR_ID from NGSREP
#'
#' @param lims_donor_id LIMS_DONOR_ID
#' @param gene One of \sQuote{A}, \sQuote{B}, \sQuote{C},
#' \sQuote{DRB1}, \sQuote{DQB1}, or \sQuote{DPB1}.
#' @param chunksize Submit queries in chunks of LIMS_DONOR_IDs.
#' @param ncores Parallelise queries.
#' @return A table with fields:
#' \itemize{
#'  \item{\code{lims_donor_id}}
#'  \item{\code{ymd}}
#'  \item{\code{conc}}
#'  \item{\code{provenance}}
#' }
#'
#' @keywords internal
#' @export
dkms_dna_concentration <- function(lims_donor_id, gene, chunksize = 1000, ncores = 1) {
  if (ncores > 1 && !requireNamespace("doParallel", quietly = TRUE)) {
    stop("doParallel package required if ncores > 1", call. = FALSE)
  }
  if (missing(lims_donor_id)) {
    stop("Argument 'lims_donor_id' required", call. = FALSE)
  }
  if (missing(gene)) {
    stop("Argument 'gene' required", call. = FALSE)
  }
  gene <- match_hla_gene(gene)
  if (chunksize > 1000) {
    warning("The maximum chunksize allowed is 1000", call. = FALSE, immediate. = TRUE)
    chunksize <- 1000
  }
  if (ncores > 1) {
    doParallel::registerDoParallel(ncores)
  } else {
    foreach::registerDoSEQ()
  }
  fmt <- "
  SELECT
    N.LIMS_DONOR_ID, TRUNC(M.DATETIME) AS YMD,
    AVG(M.DNA_CONC) AS CONC, A.AUFTRAGGEBER AS PROVENANCE
  FROM LIMSREP.NEXTYPE_DONOR_BEFU N
  INNER JOIN NGSREP.MV_DNA_MESSUNG M
    ON N.PLATE_POS = M.PLATE_POS
  INNER JOIN LIMSREP.LIMS_PROBENTRACKING T
    ON N.SRC_PANEL          = T.SOURCE_PANEL
    AND T.TARGET_PANEL_NAME = M.DNA_PANEL
  INNER JOIN LIMSREP.LIMS_AUFTRAEGE A
    ON A.LIMS_DONOR_ID = N.LIMS_DONOR_ID
  WHERE N.GENE              = '%s'
  AND N.LIMS_DONOR_ID       IN (%s)
  AND N.FREIGABE            IS NOT NULL
  GROUP BY N.LIMS_DONOR_ID, TRUNC(M.DATETIME), A.AUFTRAGGEBER
  ORDER BY N.LIMS_DONOR_ID
  "
  rs <- dplyr::tbl_dt(foreach(ldid = ichunk(lims_donor_id, chunksize), .combine = "rbind") %dopar% {
    rs <- orcl::ora_query(sprintf(fmt, gene, comma(wrap(ldid))))
    rs[, ymd := lubridate::floor_date(ymd, unit = "day")]
    data.table::setkeyv(rs, "lims_donor_id")
    rs
  })
  rs
}

#' Extract a random sample of DNA concentrations.
#'
#' @description
#' THIS FUNCTION REQURIES ACCESS TO INTERNAL DATABASES AT DKMS LSL.
#'
#' Use the prepackaged dataset (\code{\link{concentration}}) instead.
#'
#' @param x A <\code{\link{HLA}}> object.
#' @param n Number of samples.
#' @param ncores How many cores to we use for querying the database.
#' @return
#' A \code{data.table} with the fields:
#' \itemize{
#'  \item "lims_donor_id": LIMS_DONOR_ID
#'  \item "ymd": A timestamp.
#'  \item "conc": The DNA concentration in ng/µl.
#' }
#' @export
#' @examples
#' \dontrun{
#' ## Generate a random distribution of DNA concentrations
#' data("dpb1")
#' conc <- sample_dna_concentration(dpb1, n = 12000, ncores = 12)
#' }
sample_dna_concentration <- function(x, n = 1e4, ncores = 4) {
  assertive::assert_is_any_of(x, "HLA")
  gene <- strsplitN(x$gene, '-', 2)
  tbl <- x$get_table()
  ldid <- unique(tbl$lims_donor_id) %>%
    sample(size = n, replace = FALSE)
  rs <- dkms_dna_concentration(ldid, gene, chunksize = 1000, ncores = ncores)
  rs[conc > 0]
}

