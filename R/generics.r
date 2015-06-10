## Generics ####

#' Extract the gene name
#'
#' @param x An object.
#' @param ... Further arguments passed to methods
#' @keywords internal
#' @export
gene <- function(x, ...) UseMethod("gene")

#' Extract an exon
#' @inheritParams gene
#' @keywords internal
#' @export
exon <- function(x, ...) UseMethod("exon")

#' Extract an allele
#' @inheritParams gene
#' @keywords internal
#' @export
allele <- function(x, ...) UseMethod("allele")

#' Extract the NeXtype Basis ID
#' @inheritParams gene
#' @keywords internal
#' @export
nextype_basis_id <- function(x, ...) UseMethod("nextype_basis_id")

#' Extract the DNA Version ID
#' @inheritParams gene
#' @keywords internal
#' @export
dna_version_id <- function(x, ...) UseMethod("dna_version_id")

#' Extract samples
#'
#' @param x An object.
#' @param ... Further arguments passed to methods
#' @keywords internal
#' @export
samples <- function(x, ...) UseMethod("samples")

#' Extract errors
#'
#' @param x An object.
#' @param ... Further arguments passed to methods
#' @keywords internal
#' @export
errors <- function(x, ...) UseMethod("errors")

#' Extract table of genotype frequencies
#'
#' @param x An object.
#' @param ... Further arguments passed to methods
#' @keywords internal
#' @export
gtf_table <- function(x, ...) UseMethod("gtf_table")

#' Extract mapper function
#'
#' @param x An object.
#' @param ... Further arguments passed to methods
#' @keywords internal
#' @export
mapper <- function(x, ...) UseMethod("mapper")

#' Extract remapper function
#'
#' @param x An object.
#' @param ... Further arguments passed to methods
#' @keywords internal
#' @export
remapper <- function(x, ...) UseMethod("remapper")





