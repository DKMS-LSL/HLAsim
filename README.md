# HLAsim - Allelic Dropout Simulations

Simulation code and data behind the article:

Prediction of spurious HLA class II typing results using probabilistic classification.
Human Immunology (under review).

## Installation

HLAsim is not available from CRAN, but you can get it from Github with:

```R
# install.packages("devtools")
devtools::install_github("hadley/purrr")
devtools::install_github("gschofl/HLAsim")
```

Note that building the packages requires compiler support for C++11 features as well as the BOOST libraries.

## Example

An extended example, partially reproducing Figure 2 of the paper is included as a vignette. 
