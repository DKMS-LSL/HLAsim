# This file was generated by Rcpp::compileAttributes
# Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

hla_sort <- function(alleles) {
    .Call('HLAsim_hla_sort', PACKAGE = 'HLAsim', alleles)
}

hla_allele_to_genotype <- function(a1, a2) {
    .Call('HLAsim_hla_allele_to_genotype', PACKAGE = 'HLAsim', a1, a2)
}

field1 <- function(a) {
    .Call('HLAsim_field1', PACKAGE = 'HLAsim', a)
}

field2 <- function(a) {
    .Call('HLAsim_field2', PACKAGE = 'HLAsim', a)
}

