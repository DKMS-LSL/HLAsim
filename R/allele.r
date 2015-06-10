#' Class: eag_numbers
#'
#' Constructor for <eag_numbers>.
#'
#' @param e2, e3 EAG numbers for exon 2 and 3.
#' @seealso \code{\link{eag_allele}},  \code{\link{rep_allele}}
#' @keywords internal
#' @export
eag_numbers <- function(e2, e3) {
  stopifnot(length(e2) <= 2, length(e3) <= 2)
  e2 <- at_least_two(e2)
  e3 <- at_least_two(e3)
  structure(c(e2, e3), class = "eag_numbers")
}

#' @describeIn eag_numbers
#' @keywords internal
#' @export
print.eag_numbers <- function(x) {
  fmt <- "a:%-9sc:%-9s\nb:%-9sd:%-9s"
  cat(sprintf(fmt, x[1], x[3], x[2], x[4]), sep = "\n")
  invisible(x)
}

#' @describeIn eag_numbers
#' @keywords internal
#' @export
exon.eag_numbers <- function(x, n) {
  n <- as.character(n)
  if (n == "2") {
    x[1:2]
  } else if (n == "3") {
    x[3:4]
  }
}

#' @describeIn eag_numbers
#' @keywords internal
#' @export
`==.eag_numbers` <- function(e1, e2) {
  e1[1] == e2[1] && e1[2] == e2[2] && e1[3] == e2[3] && e1[4] == e2[4]
}

#' @describeIn eag_numbers
#' @keywords internal
#' @export
`!=.eag_numbers` <- function(e1, e2) !`==`(e1, e2)

#' @describeIn eag_numbers
#' @keywords internal
#' @export
`is.na.eag_numbers` <- function(x) {
  any(NextMethod())
}

#' Class: eag_allele
#'
#' Constructor for an <eag_allele>.
#'
#' Inherits from <list>.
#'
#' @param a1, a2 Alleles
#' @seealso \code{\link{eag_numbers}},  \code{\link{rep_allele}}
#' @keywords internal
#' @export
eag_allele <- function(a1, a2) {
  if (missing(a1) && missing(a2)) {
    return(structure(
      list(a1 = NA_character_, a2 = NA_character_),
      class = c("eag_allele", "list")
    ))
  }
  if (missing(a2) && is.genotype(a1)) {
    aa <- strsplit(a1, split = "/", fixed = TRUE)[[1]]
    return(structure(
      list(a1 = aa[1], a2 = aa[2]),
      class = c("eag_allele", "list")
    ))
  }
  a1 <- if (length(a1) == 1) a1 else hla_sort(a1)
  a2 <- if (length(a2) == 1) a2 else hla_sort(a2)
  structure(
    list(a1 = a1, a2 = a2),
    class = c("eag_allele", "list")
  )
}

#' @describeIn eag_allele
#' @keywords internal
#' @export
is.na.eag_allele <- function(x) {
  is.na(x$a1) && is.na(x$a2) || grepl(":$", x$a1) || grepl(":$", x$a2)
}

#' @describeIn eag_allele
#' @keywords internal
#' @export
print.eag_allele <- function(x, ...) {
  fmt <- "a: %s\nb: %s"
  cat(sprintf(fmt, paste0(allele(x, 1), collapse = "/"), paste0(allele(x, 2), collapse = "/")),
      sep = "\n")
  invisible(x)
}

#' @describeIn eag_allele
#' @keywords internal
#' @export
`==.eag_allele` <- function(e1, e2) {
  (e1$a1 == e2$a1 && e1$a2 == e2$a2) || (e1$a1 == e2$a2 && e1$a2 == e2$a1)
}

#' @describeIn eag_allele
#' @keywords internal
#' @export
`!=.eag_allele` <- function(e1, e2) !`==`(e1, e2)

#' @describeIn eag_allele
#' @keywords internal
#' @export
allele.eag_allele <- function(x, n) {
  n <- as.character(n)
  if (n == "1") {
    x$a1
  } else if (n == "2") {
    x$a2
  }
}

#' Class: rep_allele
#'
#' Constructor for a <rep_allele>.
#'
#' Inherits from <list>.
#'
#' @param a1, a2 Alleles
#' @seealso \code{\link{eag_numbers}},  \code{\link{eag_allele}}
#' @keywords internal
#' @export
rep_allele <- function(a1 = NA_character_, a2 = NA_character_) {
  structure(
    list(a1 = a1, a2 = a2),
    class = c("rep_allele", "list")
  )
}

#' @describeIn rep_allele
#' @keywords internal
#' @export
print.rep_allele <- function(x, ...) {
  fmt <- "a: %s\nb: %s"
  cat(sprintf(fmt, allele(x, 1), allele(x, 2)), sep = "\n")
  invisible(x)
}

#' @describeIn rep_allele
#' @keywords internal
#' @export
is.na.rep_allele <- function(x) {
  is.na(x$a1) && is.na(x$a2) || grepl(":$", x$a1) || grepl(":$", x$a2)
}

#' @describeIn rep_allele
#' @keywords internal
#' @export
`==.rep_allele` <- function(e1, e2) {
  (e1$a1 == e2$a1 && e1$a2 == e2$a2) || (e1$a1 == e2$a2 && e1$a2 == e2$a1)
}

#' @describeIn rep_allele
#' @keywords internal
#' @export
`!=.rep_allele` <- function(e1, e2) !`==`(e1, e2)

#' @describeIn rep_allele
#' @keywords internal
#' @export
allele.rep_allele <- function(x, n) {
  n <- as.character(n)
  if (n == "1") {
    x$a1
  } else if (n == "2") {
    x$a2
  }
}

