#' Generate a genotype resampling function.
#'
#' @note
#' The next step is to inject errors in the sampled genotypes:
#' \code{\link{inject_errors}}.
#' @param x A <\code{\link{HLA}}> object.
#' @param eag An <\code{\link{eag_table}}> object.
#' @return A constructor function for <\code{\link{geno_table}}> objects.
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
#' ## Sample a distribution of DNA concentrations
#' concentration <- sample_dna_concentration(dpb1.de, n = 24000, ncores = 8)
#'
#' ## Generate a sampling function
#' sample_dpb1_de <- make_genotype_sampler(dpb1.de, dpb1_eag1412)
#'
#' ## Sample genotypes
#' n <- 10000
#' bin_size <- 3
#' rs <- sample_dpb1_de(concentration, n, bin_size)
#' rs
#' }
make_genotype_sampler <- function(x, eag) {
  assertthat::assert_that(
    is(x, "HLA"),
    is(eag, "eag_tbl")
  )

  `%do%`   <- foreach::`%do%`
  icount   <- iterators::icount
  iter     <- iterators::iter
  foreach  <- foreach::foreach

  ## check if we run in knitr
  in_knitr <- isTRUE(getOption('knitr.in.progress'))

  ## Remove low resolution four-digit codes and cases of unknown alleles
  x <- clean_hla_data(x)
  xtbl <- x$get_table()

  ## create lookup tables mapping rep_alleles (before assignment of ambigutiy
  ## codes), eag alleles (after assignment of ambiguity codes), eag_num, and allele_num
  if (!in_knitr) {
    message("Creating lookup table and mappers ...\n", sep = "")
  }
  alleles <- xtbl[, .(allele1, allele2)]
  alleles <- hla_sort(unique(c(alleles$allele1, alleles$allele2)))

  ## a list of lookup tables matching alleles and EAG numbers
  ## for exons 2 and 3, jokers, and partials
  lookup  <- lookup_list(alleles, eag)
  ## a function that maps eag_nums to alleles
  map     <- make_mapper(lookup)
  ## a function that maps alleles to eag_nums
  remap   <- make_remapper(lookup)

 # xtbl[genotype == "04:02:01G/04:HJMR"]
 # allele2string(list(map(remap(string2allele("04:02:01G/04:HJMR")[[1]]))))
 # dt_remap["04:02:01G/04:HJMR"]

  ## check the genotypes in x for mappability
  unique_genotypes <- xtbl[, unique(genotype)]
  if (!in_knitr) {
    message("Ensuring mappability of ", length(unique_genotypes),
            " genotypes. This step will take a while ...\n", sep = "")
  }
  pbar <- utils::txtProgressBar(min = 0, max = length(unique_genotypes), style = 3)
  on.exit(close(pbar))
  iREP <- iter(string2allele(unique_genotypes))
  unique_remapped_genotypes <- foreach(rep = iREP, cnt = icount()) %do% {
    if (!in_knitr)
      utils::setTxtProgressBar(pb = pbar, value = cnt)
    tryCatch(map(remap(rep)), error = function(e) {
      message("Error '", e$message, "' in ", cnt)
    })
  }
  dt_remap <- data.table(
    genotype   = unique_genotypes,
    genotype2  = allele2string(unique_remapped_genotypes),
    eag_status = unlist(purrr::map(unique_remapped_genotypes, ~eag_status(remap(.))))
  )
  setkeyv(dt_remap, "genotype")

  ## check for and remove genotypes that could not
  ## be mapped reliably: NA/NA or something like 03:/06:02:02
  invalid_gtps <- dt_remap[!is.genotype(genotype2), genotype]
  setkeyv(xtbl, "genotype")
  if (length(invalid_gtps) > 0) {
    xtbl <- xtbl[!invalid_gtps]
    dt_remap <- dt_remap[!invalid_gtps]
  }

  ## Calculate genotype frequencies based on
  ## the remapped genotypes
  tbl2 <- xtbl[dt_remap]
  #tbl2[genotype == "04:02:01G/04:HJMR"]

  tbl2 <- tbl2[, genotype := NULL]
  setnames(tbl2, "genotype2", "genotype")
  tbl2 <- tbl2[, `:=`(
    allele1 = allele1(genotype),
    allele2 = allele2(genotype),
    zygosity = ifelse(allele1 == allele2, "homozygous", "heterozygous")
  )]
  setkeyv(tbl2, "lims_donor_id")
  x$set_table(tbl2)
  gtfrq <- gtf_table(x)

  ## Globals:
  ##   gtfrq
  ## Passed on in enclosing environment:
  ##   map
  ##   remap

  ## cleanup
  rm(
    `%do%`, cnt, dt_remap, eag, foreach, icount, invalid_gtps,
    iREP, iter, lookup, rep, unique_remapped_genotypes, tbl2,
    unique_genotypes, x, xtbl
  )

  structure(function(conc, n, bin_size) {
    ## Generate a random sample of genotypes based on the observed genotype frequency
    gtp <- as.character(gtfrq[, genotype])
    frq <- gtfrq[, pobs]
    rgeno <- sample(gtp, size = n, replace = TRUE, prob = frq)
    ## EAG status of the random sample
    eag_status <- gtfrq[rgeno][, eag_status]
    ## Is the genotype homozygous?
    rzygosity <- ifelse(is_homozygous(rgeno), "homozygous", "heterozygous")
    ## Sample the emperically observed DNA concentration
    ## and distribute them randomly to the generated genotypes.
    rconc <- sample(conc[, conc], size = n, replace = TRUE)
    ## Divide the range of DNA concentrations into slices.
    breaks <- slice_bins(conc, bin_size)
    bins <- cut(rconc, breaks, include.lowest = TRUE, dig.lab = 4, ordered_result = TRUE)
    structure(
      list(
        samples = dplyr::tbl_dt(data.table(
          genotype   = rgeno,
          eag_status = eag_status,
          zygosity   = rzygosity,
          conc       = rconc,
          bins       = bins,
          idx        = seq_len(n),
          key        = "idx"
        )),
        errors = dplyr::tbl_dt(data.table(
          genotype   = character(),
          eag_status = character(),
          zygosity   = character(),
          idx        = integer(),
          key        = "idx"
        ))
      ),
      ## pass along the closure containing the mapper, remapper, and gtf_table
      penv        = parent.env(environment(NULL)),
      breaks      = breaks,
      class = c("geno_table", "list")
    )
  }, class = c("genotype_sampler", "function"))
}

#' Class: geno_table
#'
#' Constructor function for a <\code{geno_table}> object. Created
#' by the factory function \code{\link{make_genotype_sampler}}.
#'
#' A function that takes a vector of DNA concentrations
#' to draw \emph{n} random genotypes based on observed genotype
#' frequencies and randomly attributes DNA concentrations to these
#' genotypes.
#'
#' @usage geno_table(conc, n, bin_size)
#' @return
#' A <\code{geno_table}> object with the slots:
#' \itemize{
#'   \item "samples": table with the fields:
#'     \itemize{
#'       \item "genotype": The original genotype in format: \dQuote{03:FYKD/04:ADCGE}.
#'       \item "zygosity": One of \sQuote{heterozygous} or \sQuote{homozygous}.
#'       \item "conc": Randomly attributed DNA concentrations.
#'       \item "bins": Intervals of DNA concentrations used to attribute
#'         concentration-dependent errors.
#'       \item "idx": \strong{key}; A running index.
#'     }
#'   \item "errors": table with the fields:
#'     \itemize{
#'       \item "genotype": The genotype after error injection.
#'       \item "zygosity": The zygosity after error injection.
#'       \item "idx": \strong{key}; A running index.
#'     }
#' }
#' and the attributes:
#' \itemize{
#'  \item "penv": The enclosing environment of the <geno_table> constructor
#'  containing memoised mapping and remapping functions created by
#'  \code{\link{make_mapper}} and a table of observed and expected genotype
#'  frequencies generated by \code{\link{gtf_table}}.
#'  \item "breaks": The cut points used to slice up the DNA concentration.
#' }
#' @name geno_table
#' @examples
#' \dontrun{
#' data(dpb1)
#' dpb1.pl <- dpb1[provenance == "PL"]
#'
#' ## Generate a distribution of DNA concentrations
#' conc <- sample_dna_concentration(dpb1.pl, n = 1000, ncores = 8)
#'
#' ## Generate a sampling function
#' data(dpb1_eag1412)
#' geno_table <- make_genotype_sampler(dpb1.pl, dpb1_eag1412)
#'
#' ## Access mapping and remapping functions as well as genotype frequencies
#' mapper(geno_table)
#' remapper(geno_table)
#' gtf_table(genotable)
#'
#' ## Sample genotypes
#' n <- 10000
#' bin_size <- 3
#' ans <- geno_table(conc, n, bin_size)
#' ans
#'
#' summary(ans)
#' samples(ans)
#' }
NULL

#' @export
print.geno_table <- function(x, n = 5, ...) {
  cat("<geno_table> [samples: ", nrow(x$samples), "; errors: ", nrow(x$errors), "]\n\n", sep = "")
  cat(" Samples: ", sep = "")
  print(x$samples, n = n)
  cat("\n Errors: ", sep = "")
  print(x$errors, n = n)
}

#' @export
summary.geno_table <- function(object, ...) {
  j <- samples(object)[errors(object)]
  j[zygosity == "heterozygous", .(
    pHom = sum(i.zygosity == "homozygous", na.rm = TRUE)/nrow(.SD),
    pHet = sum(i.zygosity == "heterozygous", na.rm = TRUE)/nrow(.SD),
    pHetIdent = sum(i.zygosity == "heterozygous" & genotype == i.genotype, na.rm = TRUE)/nrow(.SD),
    pHetDiff  = sum(i.zygosity == "heterozygous" & genotype != i.genotype, na.rm = TRUE)/nrow(.SD),
    pNA  = sum(is.na(i.zygosity))/nrow(.SD),
    nerr = nrow(.SD)
  )]
}

#' @export
merge.geno_table <- function(x, y = NULL) {
  samples(x)[errors(x)] %>%
    dplyr::select(genotype, i.genotype, zygosity, i.zygosity, eag_status, i.eag_status, conc, bins) %>%
    rename(genotype2 = i.genotype, zygosity2 = i.zygosity, eag_status2 = i.eag_status)
}

#' @export
samples.geno_table <- function(x) {
  x$samples
}

#' @export
errors.geno_table <- function(x) {
  x$errors
}

#' @export
gtf_table.geno_table <- function(x) {
  get("gtfrq", envir = attr(x, "penv"))
}

#' @export
gtf_table.genotype_sampler <- function(x) {
  get("gtfrq", envir = environment(x))
}

#' @export
mapper.geno_table <- function(x) {
  get("map", envir = attr(x, "penv"))
}

#' @export
mapper.genotype_sampler <- function(x) {
  get("map", envir = environment(x))
}

#' @export
remapper.geno_table <- function(x) {
  get("remap", envir = attr(x, "penv"))
}

#' @export
remapper.genotype_sampler <- function(x) {
  get("remap", envir = environment(x))
}

slice_bins <- function(conc, bin_size) {
  if (is(conc, "data.table")) {
    assertthat::assert_that('conc' %in% names(conc))
    conc <- conc[, conc]
  }
  conc_max <- ceiling(max(conc))
  unique(c(seq(0, conc_max, bin_size), conc_max))
}

