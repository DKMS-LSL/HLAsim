# HLAsim - Simulating HLA locus allelic dropout 

Simulation code and data associated with the _Human Immumology_ article:

__Schöfl G, Schmidt AH, Lange V (2016) Prediction of spurious HLA class II typing results using probabilistic classification.__ _Human Immunology, accepted for publication_.

The article is available online, [doi:10.1016/j.humimm.2016.01.012](http://dx.doi.org/10.1016/j.humimm.2016.01.012)

## Abstract

While modern high-throughput sequence-based HLA genotyping methods generally provide highly accurate typing results, artefacts may nonetheless arise for numerous reasons, such as sample contamination, sequencing errors, read misalignments, or PCR amplification biases. To help detecting spurious typing results, we tested the performance of two probabilistic classifiers (binary logistic regression and random forest models) based on population-specific genotype frequencies. We trained the model using high-resolution typing results for HLA-DRB1, DQB1, and DPB1 from large samples of German, Polish and UK-based donors. The high predictive capacity of the best models replicated both in 10-fold cross-validation for each gene and in using independent evaluation data (AUC 0.820–0.893). While genotype frequencies alone provide enough predictive power to render the model generally useful for highlighting potentially spurious typing results, the inclusion of workflow-specific predictors substantially increases prediction specificity. Low initial DNA concentrations in combination with low-volume PCR reactions form a major source of stochastic error specific to the Fluidigm chip-based workflow at DKMS Life Science Lab. The addition of DNA concentrations as a predictor variable thus substantially increased AUC (0.947–0.959) over purely frequency-based models. 

## Installation

HLAsim is not available from CRAN, but you can get it from Github with:

```R
# install.packages("devtools")
devtools::install_github("dkms-lsl/HLAsim")
```

Note that building the packages requires compiler support for C++11 features as well as the BOOST libraries.

## Example

An extended example, partially reproducing Figure 2 of the paper is included as a vignette. 
