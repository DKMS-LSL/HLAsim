multinomial_coefficient <- function(x) {
  n <- nrow(x)
  coef <- numeric(n)
  for (i in seq_len(n)) coef[i] <- 2/prod(factorial(tabulate(x[i, ])))
  coef
}

expand_genotype <- function(alleles, as.int = FALSE, matrix = TRUE) {
  n <- length(alleles)
  first.allele <- rep.int(1:n, times = n:1)
  seq.n <- purrr::partial(seq.int, to = n)
  second.allele <- unlist(lapply(1:n, seq.n))
  rs <- cbind(first.allele, second.allele)
  dimnames(rs) <- NULL
  if (!as.int)
    rs <- matrix(alleles[rs], ncol = 2)
  if (!matrix) {
    rs <- apply(rs, 1, paste0, collapse = "/")
  }
  rs
}

#' Calculate expected genotype frequencies.
#'
#' @param alleles A character vector of alleles.
#' @param frequencies The corresponding allele frequencies.
#' @return A data table with two columns:
#' \itemize{
#'  \item{genotype}
#'  \item{pexp}
#' }
#' @keywords internal
expected_genotype_frequencies <- function(alleles, frequencies) {
  ## Make sure that alleles and frequencies are properly sorted so
  ## that they yield correct genotypes
  dt <- data.table(
    a = as.character(alleles),
    f = as.numeric(frequencies),
    key = "a"
  )[hla_sort(as.character(alleles))]
  genos <- expand_genotype(dt$a, as.int = TRUE, matrix = TRUE)
  p <- multinomial_coefficient(genos) * apply(genos, 1, function(i) prod(dt$f[i]))
  np <- apply(matrix(dt$a[genos], ncol = 2), 1L, paste, collapse = "/")
  dplyr::tbl_dt(data.table(genotype = np, pexp = p, key = "genotype"))
}
