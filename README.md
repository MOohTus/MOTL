
<!-- README.md is generated from README.Rmd. Please edit that file -->

# MOTL

<!-- badges: start -->

<!-- badges: end -->

MOTL infers latent factor values for a multi-omics target dataset,
consisting of a small number of samples, by incorporating latent factor
values already inferred with a MOFA factorization of a large,
heterogeneous, learning dataset.

## Installation

You can install the development version of MOTL from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
# devtools::install_github("MOohTus/MOTL-pkg")
```

## How to use MOTL

### Preprocessing of learning dataset factorization

``` r
LrnDir <- "LrnData"
```

### Preproessing of target dataset factorization

``` r
YTrg_list <- list(
  mRNA = "expdat_mRNA",
  miRNA = "expdat_miRNA",
  DNAme = "expdat_DNAme",
  SNV = "expdat_SNV"
)
```

### Transfer learning factorization with MOTL

``` r
TL_data <- ""
```
