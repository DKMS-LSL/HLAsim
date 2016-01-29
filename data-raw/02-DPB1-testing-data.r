
# DPB1 --------------------------------------------------------------------

library(devtools)
library(HLAdata)
library(HLAtest)
library(HLAsim)

get_testdata <- HLAsim:::get_testdata
clean_testing_data <- HLAsim:::clean_testing_data

## DNA version ID 52  ##
HLA_DPB1_TEST_52 <- get_testdata(gene = "DPB1", nextype_basis_id = HLAdata::dna2basis("52"))
tdnew <- HLA_DPB1_TEST_52
tdnew2 <- clean_testing_data(td = tdnew)
HLA_DPB1_TEST_52 <- tdnew2
devtools::use_data(HLA_DPB1_TEST_52, pkg = "~/Devel/HLAtest", overwrite = TRUE)

## DNA version ID 53  ##
HLA_DPB1_TEST_53 <- get_testdata(gene = "DPB1", nextype_basis_id = HLAdata::dna2basis("53"))
tdnew <- HLA_DPB1_TEST_53
tdnew2 <- clean_testing_data(td = tdnew)
HLA_DPB1_TEST_53 <- tdnew2
devtools::use_data(HLA_DPB1_TEST_53, pkg = "~/Devel/HLAtest", overwrite = TRUE)

## DNA version ID 56  ##
HLA_DPB1_TEST_56 <- get_testdata(gene = "DPB1", nextype_basis_id = HLAdata::dna2basis("56"))
tdnew <- HLA_DPB1_TEST_56
tdnew2 <- clean_testing_data(td = tdnew)
HLA_DPB1_TEST_56 <- tdnew2
devtools::use_data(HLA_DPB1_TEST_56, pkg = "~/Devel/HLAtest", overwrite = TRUE)

## DNA version ID 58  ##
HLA_DPB1_TEST_58 <- get_testdata(gene = "DPB1", nextype_basis_id = HLAdata::dna2basis("58"))
tdnew <- HLA_DPB1_TEST_58
tdnew2 <- clean_testing_data(td = tdnew)
HLA_DPB1_TEST_58 <- tdnew2
devtools::use_data(HLA_DPB1_TEST_58, pkg = "~/Devel/HLAtest", overwrite = TRUE)

