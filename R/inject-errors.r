#' Inject errors into a sample of genotypes.
#'
#' Takes a \code{\link{geno_table}} object and injects a total of
#' \emph{perr} errors with a concentration-dependent rate described
#' by a cumulative distribution functions, \emph{cdf}, using the
#' probability calculation described in \code{\link{loose_eag}}..
#'
#' @note The next step is to map genotype frequencies onto the sampled
#'   genotypes: \code{\link{map_genotype_frequencies}}.
#' @param x A \code{\link{geno_table}} object with attribute \code{has_errors = FALSE}.
#' @param cdf A cumulative distribution function describing a DNA-concentration-
#'   dependent error rate. See \code{\link{cdfA}} and \code{\link{cdfB}} for two
#'   precomputed error functions based on HLA-A and HLA-B.
#' @param perr The overall error rate.
#' @param odds The odds of loosing an EAG on exon 2 or exon 3, respectively.
#' E.g. \code{odds = 2} means a 2:1 chance that the loss occurs on Exon 2;
#' \code{odds = 0.25} means a 1:4 chance that the loss occurs on Exon 3.
#'
#' @return A \code{\link{geno_table}} object with attribute \code{has_errors = TRUE}.
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
#' sample_dpb1_de <- make_genotype_sampler(dpb1.de, eag_1267)
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
#'   inject_errors(cdfA, perr, odds)
#' ans
#'
#' summary(ans)
#' samples(ans)
#' errors(ans)
#' merge(ans)
#' }
inject_errors <- function(x, cdf, perr, odds = 1) {
  assertthat::assert_that(
    is(x, "geno_table")
  )

  if (nrow(errors(x)) > 0) {
    warning("This <geno_table> has already errors injected", immediate. = TRUE)
    return(x)
  }

  `%do%`  <- foreach::`%do%`
  iter    <- iterators::iter
  foreach <- foreach::foreach

  ## Fetch mappers from the enclosing environment of the
  ## <geno_table> constructor. Make sure they get re-assigned
  ## after use to take full advantage of the memoisation.
  mapargs <- list(
    map   = get("map", envir = attr(x, "penv")),
    remap = get("remap", envir = attr(x, "penv")),
    odds  = odds
  )
  on.exit(assign("remap", mapargs$remap, envir = attr(x, "penv")))
  on.exit(assign("map", mapargs$map, envir = attr(x, "penv")), add = TRUE)

  ## Number of samples
  n <- nrow(samples(x))
  ## Total number of mistyped samples
  nerr <- ceiling(n*perr)
  ## Conditional probability of an error at some slice of DNA
  ## concentrations given an error occurred
  cond_perr <- diff(cdf(attr(x, "breaks")))
  cond_perr <- setNames(object = cond_perr, nm = levels(samples(x)[, bins]))
  cond_perr <- cond_perr[cond_perr > 0]
  ## Compute the number of errors per slice of DNA concentration
  sampling_size <- round(nerr*cond_perr, digits = 0)

  ## Sample heterozygotes that will be homogenised, merge them with the true
  ## homozygotes and join the resulting table with the genotype frequencies.
  iSS <- iter(sampling_size)
  iSN <- iter(names(sampling_size))
  samples <- samples(x)
  #ss <- nextElem(iSS)
  #sn <- nextElem(iSN)
  rs <- foreach(ss = iSS, sn = iSN, .combine = "rbind") %do% {
    dt <- samples[bins == sn]
    if ((nr <- nrow(dt)) < ss) {
      fmt <- paste0(
        "\n      The number of errors (%s) in bin %s is larger",
        " than the number of available samples (%s).\n",
        "Try decreasing the overall error rate (perr = %s)")
      msg <- sprintf(gsub(" +", " ", trimws(fmt)), ss, sn, nrow(dt), perr)
      warning(msg, call. = TRUE)
      ss <- nrow(dt)
    }
    dt[sample.int(nr, size = ss, replace = FALSE)]
  }
  setkeyv(rs, "bins")
  rs2 <- data.table(pE = cond_perr, bins = names(cond_perr), key = "bins")
  dt <- rs[rs2, allow.cartesian = TRUE]
  dt <- dt[!is.na(genotype)]
  dt[, c("genotype", "eag_status") := do_inject(dt, mapargs)]
  dt[, zygosity := ifelse(is_homozygous(genotype), "homozygous", "heterozygous")]
  dt[, genotype := ifelse(genotype == "NA/NA", NA_character_, genotype)]
  dt[, zygosity := ifelse(is.na(genotype), NA_character_,  zygosity)]
  dt[, conc := NULL]
  dt[, bins := NULL]
  dt[, pE := NULL]
  setkeyv(dt, "idx")
  x$errors <- dplyr::tbl_dt(dt)
  x
}

EA <- function(E, O) {
  0.5*(1 + O - O*sqrt((1 + 1/O)^2 - (4*E)/O))
}

#' Loose an EAG
#'
#' Randomly loose an EAG number.
#'
#' @details
#' We generally deal with the following scenario on any HLA locus currently:
#'
#' \tabular{rcc}{
#'         \tab Exon A \tab Exon B \cr
#'   EAG 1 \tab a      \tab c      \cr
#'   EAG 2 \tab b      \tab d
#' }
#'
#' We know (have a guesstimate of) the total error rate due to loss of EAGs
#' at a given DNA concentration, \eqn{P(E)}, and we know the odds, \eqn{\omega},
#' that we'll loose an EAG at one or the other exon.
#'
#' E.g. if we perform two independent PCR reactions for exon 2 and only one
#' PCR reaction for exon 3, the odds of loosing an EAG are 4 in 1 for exon
#' 3. I. e. in order for an error to manifest on exon 2 we must loose the
#' same EAG in both reactions, the probability of which is \out{&frac14;}
#' of loosing an EAG in a single reaction.
#'
#' What we want to know is:
#' \itemize{
#'   \item \out{P(A &cap; &not; B | E)}
#'   The probability that we loose an EAG on exon A but not on exon B,
#'   given that we make an error at all.
#'   \item \out{P(&not; A &cap; B | E)}
#'   The probability that we loose an EAG on exon B but not on exon A,
#'   given that we make an error at all.
#'   \item \out{P(A &cap; B | E)}
#'   The probability that we loose an EAG on exon A and on exon B,
#'   given that we make an error at all.
#' }
#'
#' To get there we need to know the error rate for exon A, \out{&epsilon;<sub>A</sub>},
#' and the error rate for exon B, \out{&epsilon;<sub>B</sub>}, since:
#' \itemize{
#'  \item \out{P(A &cap; &not; B; | E) = &epsilon;<sub>A</sub>(1 - &epsilon;<sub>B</sub>)
#'                                       &frasl;
#'                                       (&epsilon;<sub>A</sub> + &epsilon;<sub>B</sub> -
#'                                       &epsilon;<sub>A</sub>&epsilon;<sub>B</sub>)}
#'  \item \out{P(&not; A &cap; B; | E) = &epsilon;<sub>B</sub>(1 - &epsilon;<sub>A</sub>)
#'                                       &frasl;
#'                                       (&epsilon;<sub>A</sub> + &epsilon;<sub>B</sub> -
#'                                       &epsilon;<sub>A</sub>&epsilon;<sub>B</sub>)}
#'  \item \out{P(&not; A &cap; B; | E) = &epsilon;<sub>B</sub>&epsilon;<sub>A</sub>
#'                                       &frasl;
#'                                       (&epsilon;<sub>A</sub> + &epsilon;<sub>B</sub> -
#'                                       &epsilon;<sub>A</sub>&epsilon;<sub>B</sub>)}
#' }
#'
#' Considering that \out{&epsilon;<sub>A</sub> + &epsilon;<sub>B</sub> -
#' &epsilon;<sub>A</sub>&epsilon;<sub>B</sub> = P(E)} and
#' \out{&epsilon;<sub>B</sub> = <sup>1</sup>&frasl;<sub>&omega;</sub>
#' &epsilon;<sub>A</sub>} we can arrive after some rearrangements at:
#'
#' \itemize{
#'  \item \out{&epsilon;<sub>A</sub> = &frac12;(1 + &omega; - &omega;
#'  &radic;<span style="text-decoration: overline">
#'  (1 + 1&frasl;&omega;)<sup>2</sup> -
#'  4<em>P(E)</em>&frasl;&omega;</span>)}
#'  \item \out{&epsilon;<sub>B</sub> = <sup>1</sup>&frasl;<sub>&omega;</sub>
#' &epsilon;<sub>A</sub>}
#'  }
#'
#'  These error rates allow now to randomly loose EAGs with the
#'  correct probabilities.
#'
#' @param x A <\code{\link{eag_numbers}}> instance.
#' @param pE The absolute probability that an error occurs.
#' @param odds The odds in favour of loosing an EAG on Exon 2 rather
#' than Exon 3. E.g. \code{odds = 2} means a 2:1 chance that the loss
#' occurs on Exon 2; \code{odds = 0.25} means a 1:4 chance that
#' the loss occurs on Exon 3.
#' @return An <eag_numbers> instance
#' @keywords internal
loose_eag <- function(x, pE, odds = 1) {
  ea   <- EA(pE, odds)
  eb   <- (1/odds)*ea
  P_a  <- ea*(1 - eb)/pE
  P_b  <- eb*(1 - ea)/pE
  P_ab <- (ea*eb)/pE
  ## Cast the die to determine if we are going to sample both exons
  exons <- sample.int(2L, 1L, replace = FALSE, prob = c(1 - P_ab, P_ab))
  if (exons == 2L) {
    ## Sample one EAG from each exon to loose and keep the other one
    loose <- sample.int(2L, 2L, replace = TRUE, prob = c(0.5, 0.5))
    keep <- ifelse(loose == 1, 2, 1)
    loose[2] <- loose[2] + 2
    keep[2] <- keep[2] + 2
  }
  else {
    ## Sample one of the four EAGs to loose and keep the other one
    loose <- sample.int(4L, 1L, replace = FALSE, prob = c(P_a/2, P_a/2, P_b/2, P_b/2))
    keep <- switch(loose, 2, 1, 4, 3)
  }
  x[loose] <- x[keep]
  x
}

#' How many eags do we find on exons 2 and 3?
#'
#' 11 <-> hom:hom
#' 12 <-> hom:het
#' 21 <-> het:hom
#' 22 <-> het:het
#'
#' @param nums Eag number
#' @keywords internal
eag_status <- function(nums) {
  paste0(nunique(nums[1:2]), nunique(nums[3:4]))
}

#' Inject errors
#'
#' @param dt A data.table with genotypes and error probabilities
#' @param args A list containing a mapper function <map>, a
#'   remapper function <remap>, and <odds>.
#' @keywords internal
do_inject <- function(dt, args) {
  dots <- list(a = string2allele(dt$genotype), p = dt$pE)
  data.table::rbindlist(
    .mapply(inject, dots = dots, MoreArgs = args)
  )
}

inject <- function(a, p, odds, map, remap) {
  nums <- loose_eag(remap(a), pE = p, odds = odds)
  list(slash(map(nums)), eag_status(nums))
}



