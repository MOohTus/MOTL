## CREATE ENVIRONMENT
## conda create --name tcga_data --file 00_Ressources/tcga_data_env.txt

## https://github.com/ArnaudDroitLab/DemoPaquetR/blob/master/paquets_R.md
## In Rstudio
library("usethis")
library("devtools")
library("BiocCheck")

## setwd("/home/morgane/Documents/01_Projects/03_OtherProjects/05_David_Thesis/test_package_aout25/")

## CREATE PACKAGE NAMED MOTL
## usethis::create_package("MOTL")

## ACTIVATE THE PROJECT (not always necessary ...)
setwd("/home/morgane/Documents/01_Projects/03_OtherProjects/05_David_Thesis/test_package_aout25/MOTL/")
usethis::proj_set(".", force = TRUE)

## UNIT TEST
usethis::use_testthat()

## VIGNETTES
usethis::use_vignette("Vignette")

## ROXYGEN TO DOCUMENT
## LOAD COMPLETE PACKAGE
devtools::document(); devtools::load_all()

## BUILD and CHECK PACKAGE (devtools::check("../"))
devtools::build(); devtools::check_built("../MOTL_0.99.0.tar.gz")

## Check for bioconductor
BiocCheck("../MOTL_0.99.0.tar.gz")
