library(HLAsim)
library(lubridate)

# Genotyping data ---------------------------------------------------------

## Genotyping results for DRB1, DQB1, and DPB1 between 01/01/2014 and 23/03/2015
DRB1 <- HLA("DRB1", "01/01/2014", "23/03/2015")
DQB1 <- HLA("DQB1", "01/01/2014", "23/03/2015")
DPB1 <- HLA("DPB1", "01/01/2014", "23/03/2015")

## remove low-resolution four-digit codes and unknown alleles
drb1 <- HLAsim:::clean_hla_data(DRB1)
dqb1 <- HLAsim:::clean_hla_data(DQB1)
dpb1 <- HLAsim:::clean_hla_data(DPB1)

devtools::use_data(drb1, overwrite = TRUE)
devtools::use_data(dqb1, overwrite = TRUE)
devtools::use_data(dpb1, overwrite = TRUE)

# EAG tables --------------------------------------------------------------

drb1_eag1412 <- eag_table(gene = "DRB1", nextype_basis_id = "1412")
devtools::use_data(drb1_eag1412, overwrite = TRUE)

dqb1_eag1412 <- eag_table(gene = "DQB1", nextype_basis_id = "1412")
devtools::use_data(dqb1_eag1412, overwrite = TRUE)

dpb1_eag1412 <- eag_table(gene = "DPB1", nextype_basis_id = "1412")
devtools::use_data(dpb1_eag1412, overwrite = TRUE)

# Error distribution for HLA-A and HLA-B ----------------------------------

A <- HLA::lost_heterozygotes("A", "01/09/2014", "23/03/2015")
cdfA <- HLA::concentration_dependent_error_distribution(A, chunksize = 500, ncores = 8)
devtools::use_data(cdfA, overwrite = TRUE)

B <- HLA::lost_heterozygotes("B", "01/09/2014", "23/03/2015")
cdfB <- HLA::concentration_dependent_error_distribution(B, chunksize = 500, ncores = 8)
devtools::use_data(cdfB, overwrite = TRUE)

# DNA concentrations ------------------------------------------------------

## Sample a distribution of DNA concentrations
conc.de <- sample_dna_concentration(dqb1[provenance == "DE"], n = 24000, ncores = 12)
conc.de <- conc.de[ymd >= ymd("2014-01-01")]
conc.de[, provenance := "DE"]

conc.pl <- sample_dna_concentration(dqb1[provenance == "PL"], n = 24000, ncores = 12)
conc.pl <- conc.pl[ymd >= ymd("2014-01-01")]
conc.pl[, provenance := "PL"]

conc.uk <- sample_dna_concentration(dqb1[provenance == "UK"], n = 24000, ncores = 12)
conc.uk <- conc.uk[ymd >= ymd("2014-01-01")]
conc.uk[, provenance := "UK"]

concentration <- tbl_dt(rbind(conc.de, conc.pl, conc.uk))
devtools::use_data(concentration, overwrite = TRUE)

# Jokers and Partials -----------------------------------------------------

drb1_jokers1412 <- joker_table(drb1_eag1412)
dqb1_jokers1412 <- joker_table(dqb1_eag1412)
dpb1_jokers1412 <- joker_table(dpb1_eag1412)
devtools::use_data(drb1_jokers1412, overwrite = TRUE)
devtools::use_data(dqb1_jokers1412, overwrite = TRUE)
devtools::use_data(dpb1_jokers1412, overwrite = TRUE)

drb1_partials1412 <- partials_table(drb1_eag1412)
dqb1_partials1412 <- partials_table(dqb1_eag1412)
dpb1_partials1412 <- partials_table(dpb1_eag1412)
devtools::use_data(drb1_partials1412, overwrite = TRUE)
devtools::use_data(dqb1_partials1412, overwrite = TRUE)
devtools::use_data(dpb1_partials1412, overwrite = TRUE)


