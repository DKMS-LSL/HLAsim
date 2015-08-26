
basis2dna_ <- map_nextype_basis_to_dna_version()
devtools::use_data(basis2dna_, overwrite = TRUE)

dna2basis_ <- map_dna_version_to_nextype_basis()
devtools::use_data(dna2basis_, overwrite = TRUE)
