#' Genotyping data generated at DKMS Life Science Lab.
#'
#' @description
#' An R6 reference class holding genotyping results for any of the
#' six HLA genes routinely typed at DKMS LSL.
#'
#' THIS FUNCTION REQURIES ACCESS TO INTERNAL DATABASES AT DKMS LSL.
#'
#' Use the prepackaged datasets (\code{\link{geno_data}}) instead.
#'
#' @param gene Fetch alleles for one of "A", "C", "B", "DRB1", "DQB1", "DPB1".
#' @param from Start query at DD/MM/YYYY.
#' @param to End query at DD/MM/YYYY.
#' @return A <HLA> object containing a table with the following fields:
#' \describe{
#'    \item{lims_donor_id}{Unique sample tracking ID.}
#'    \item{provenance}{One of \emph{DE}, \emph{PL}, or \emph{UK}.}
#'    \item{allele1}{First allele.}
#'    \item{allele2}{Other allele.}
#'    \item{genotype}{The genotype in the format "allele1/allele2", where
#'    the alleles are alphabetically ordered.}
#'    \item{zygosity}{One of \emph{homozygous} or \emph{heterozygous}}
#' }
#' @keywords internal
#' @export
#' @examples
#' \dontrun{
#' dpb1 <- HLA("DPB1", "01/01/2014", "23/04/2015")
#'
#' # restrict the sample to a country of origin
#' dpb1.de <- dpb1[provenance == "DE"]
#'
#' # access the data table
#' dpb1.de.tbl <- dpb1.de$get_table()
#'
#' # calculate allele frequencies
#' dpb1.de.af <- dpb1.de$allele_frequency()
#'
#' # calculate genotype frequencies
#' dpb1.de.gf <- dpb1.de$genotype_frequency()
#' }
HLA <- function(gene = NULL, from = NULL, to = NULL) {
  .HLA$new(gene, from, to)
}

#' @export
"[.HLA" <- function(x, i, ...) {
  out <- HLA()
  out$gene <- x$gene
  out$time_from <- x$time_from
  out$time_to <- x$time_to
  tbl <- x$get_table()
  out$set_table(tbl[eval(substitute(i)), ...])
  invisible(out)
}

#' @export
print.HLA <- function(x, ...) {
  x$print(...)
}

# HLA-class --------------------------------------------------------------

.HLA <- R6::R6Class(
  classname = "HLA",
  public = list(
    gene      = NULL,
    time_from = NULL,
    time_to   = NULL,
    initialize = function(gene = NULL, from = NULL, to = NULL) {
      if ( !is.null(gene) &&
          (!is.null(from) && valid_date(from)) &&
          (!is.null(to)   && valid_date(to))) {
        self$gene <- match_hla_gene(gene)
        self$time_from <- from
        self$time_to <- to
        private$db <- orcl::Ngsread()
        private$fetch()
      }
    },
    get_table = function() {
      private$allele_table
    },
    set_table = function(tbl) {
      private$allele_table <- tbl
    },
    refresh = function(gene = NULL, from = NULL, to = NULL) {
      self$gene <- if (!is.null(gene)) match_hla_gene(gene) else self$gene
      self$time_from <- if (!is.null(from) && valid_date(from)) from else self$time_from
      self$time_to <- if (!is.null(to) && valid_date(to)) to else self$time_to
      private$db$connect()
      private$fetch()
      private$allele_freq <- data.table()
      private$genotype_freq <- data.table()
    },
    print = function(...) {
      cat(private$to_string(), sep = "")
      print(private$allele_table)
      invisible(self)
    }
  ),
  private = list(
    ## Fields
    db            = NULL,
    allele_table  = data.table(),
    allele_freq   = data.table(),
    genotype_freq = data.table(),
    ## Functions
    to_string = function(...) {
      sprintf("HLA genotyping data - Locus <%s>\n", self$gene)
    },
    has_allele_freq = function() {
      nrow(private$allele_freq) > 0L
    },
    has_genotype_freq = function() {
      nrow(private$genotype_freq) > 0L
    },
    length = function() {
      nrow(private$allele_table)
    }
  )
)

.HLA$set("private", "fetch", function() {
  gene <- self$gene
  from <- self$time_from
  to   <- self$time_to
  if (is.null(self$gene)) {
    return(NULL)
  }

  cat("Fetching data from <lims_befunde> for locus <", gene, ">",
      if (!is.null(from)) paste0(" from ", from),
      if (!is.null(to)) paste0(" to ", to),
      sep ="")

  fmt <- "
  SELECT
  b.LIMS_DONOR_ID AS lims_donor_id,
  (CASE a.AUFTRAGGEBER
    WHEN 'DKMS'   THEN 'DE'
    WHEN 'NKR'    THEN 'DE'
    WHEN 'DKMSUK' THEN 'UK'
    WHEN 'DKMSPL' THEN 'PL'
    WHEN 'DKMSUS' THEN 'US'
    WHEN 'EUROD'  THEN 'NL'
  END) AS provenance,
  b.ALLELE1 AS allele1,
  b.ALLELE2 AS allele2
  FROM LIMSREP.LIMS_BEFUNDE b
  INNER JOIN LIMSREP.LIMS_AUFTRAEGE a
  ON b.LIMS_DONOR_ID = a.LIMS_DONOR_ID
  WHERE b.GENE         = '%s'%s
  AND a.AUFTRAGGEBER IN ('DKMS', 'DKMSUK', 'DKMSPL', 'DKMSUS', 'EUROD', 'NKR')
  AND (b.AUFLOESUNG  = 'H0' OR b.AUFLOESUNG = 'H1')
  AND (UPPER(b.ALLELE1) != 'NEW' OR UPPER(b.ALLELE2) != 'NEW')
  AND a.AUFTRAGSART  = 'NEUSPENDER'
  "

  time_constraint <- if (!is.null(from) && !is.null(to)) {
    sprintf("\nAND b.UEBERMITTLUNGSDATUM BETWEEN TO_DATE('%s', 'DD/MM/YYYY') AND TO_DATE('%s', 'DD/MM/YYYY')",
            from, to)
  } else if (!is.null(from)) {
    sprintf("\nAND b.UEBERMITTLUNGSDATUM >= TO_DATE('%s', 'DD/MM/YYYY')", from)
  } else if (!is.null(to)) {
    sprintf("\nAND b.UEBERMITTLUNGSDATUM <= TO_DATE('%s', 'DD/MM/YYYY')", to)
  } else ""

  q <- sprintf(fmt, gene, time_constraint)
  rs <- private$db$query(q)
  on.exit(private$db$disconnect())

  ## remove any incomplete samples
  rs <- rs[!is.na(allele1) & !is.na(allele2)]
  ## Determine genotype taking care to sort Allele1 and Allele2 alphabetically
  ## before joining them together.
  rs[, genotype := hla_allele_to_genotype(allele1, allele2)]
  ## Deterimine zygosity
  rs[, zygosity := ifelse(allele1 == allele2, 'homozygous', 'heterozygous')]
  data.table::setkeyv(rs, 'lims_donor_id') # sorting by key
  private$allele_table <- dplyr::tbl_dt(rs)
})

.HLA$set("public", "allele_frequency", function() {
  if (!private$has_allele_freq()) {
    rs <- data.table::copy(self$get_table())
    rs[, `:=`(provenance = NULL, genotype = NULL, eag_status = NULL)]
    rs <- reshape2::melt(rs, id.vars = c("lims_donor_id", "zygosity"), value.name = "allele")
    data.table::setkeyv(rs, 'lims_donor_id')
    n <- nrow(rs)
    rs <- rs[, list(
      count     = .N,
      frequency = .N/n
    ), by = "allele"][order(-count)]
    rs[, cumsum := cumsum(frequency)]
    rs[, allele := factor(allele, levels = allele, ordered = TRUE)]
    private$allele_freq <- dplyr::tbl_dt(rs)
  }
  private$allele_freq
})

.HLA$set("public", "genotype_frequency", function() {
  if (!private$has_genotype_freq()) {
    rs <- data.table::copy(self$get_table())
    rs[, `:=`(provenance = NULL)]
    ## sort by genotype
    data.table::setkeyv(rs, 'genotype')
    n <- nrow(rs)
    rs <- rs[, list(
      count     = .N,
      frequency = .N/n
    ), by = c("genotype")][order(-count)]
    rs[, cumsum := cumsum(frequency)]
    alleles <- strsplit(rs[, genotype], "/", fixed = TRUE)
    rs[, `:=`(
      allele1 = vapply(alleles, "[", 1L, FUN.VALUE = character(1)),
      allele2 = vapply(alleles, "[", 2L, FUN.VALUE = character(1))
    )]
    rs[, genotype := factor(genotype, levels = genotype, ordered = TRUE)]
    private$genotype_freq <- dplyr::tbl_dt(rs)
  }
  private$genotype_freq
})

.HLA$set("public", "purge_frequencies", function() {
  if (private$has_genotype_freq()) {
    private$genotype_freq <- data.table()
  }
  if (private$has_allele_freq()) {
    private$allele_freq <- data.table()
  }
})
