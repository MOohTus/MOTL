## CREATE ENVIRONMENT
## conda create --name tcga_data --file 00_Ressources/tcga_data_env.txt

## https://github.com/ArnaudDroitLab/DemoPaquetR/blob/master/paquets_R.md
## In Rstudio
library("usethis")
library("devtools")
library("BiocCheck")

.libPaths(c(.libPaths(), "/home/morgane/Documents/00_Tools/miniconda3/envs/tcga_data/"))
.libPaths(c(.libPaths(), "/home/morgane/R/x86_64-pc-linux-gnu-library/4.2/"))

## setwd("/home/morgane/Documents/01_Projects/03_OtherProjects/05_David_Thesis/test_package_aout25/")

## CREATE PACKAGE NAMED MOTL
## usethis::create_package("MOTL")

## ACTIVATE THE PROJECT (not always necessary ...)
setwd("/home/morgane/Documents/01_Projects/03_OtherProjects/05_David_Thesis/test_package_aout25/MOTL/")
usethis::proj_set(".", force = TRUE)
## SITUATION REPORT OF PROJECT PATHS
proj_sitrep()

## IMPORT SOME PACKAGES BY HAND ...
## usethis::use_import_from("stats", c("dpois", "dbinom"))

## TO IGNORE FILES WHEN CREATE THE BUNDLE
usethis::use_build_ignore("file")

## UNIT TEST
## https://r-pkgs.org/testing-basics.html
## https://r-pkgs.org/testing-design.html
## https://r-pkgs.org/testing-advanced.html
usethis::use_testthat(3)
## If the current page is a R script, this command will create a corresponding test file
usethis::use_test()
## Perfom test
devtools::test()
##  Add withr in SUGGEST, because we ued it in tests
usethis::use_package("withr", type = "Suggests")





## VIGNETTES
## https://r-pkgs.org/vignettes.html
## Install package before building vignette
devtools::install(build_vignettes = TRUE)
##
devtools::install(quick = TRUE)

## CREATE VIGNETTE FILE
usethis::use_vignette("MOTL")
usethis::use_vignette("UseCase_TCGA")
usethis::use_vignette("UseCase_Glioblastoma")





usethis::use_package_doc()

## CREATE README TEMPLATE FOR THE PACKAGE
usethis::use_readme_rmd()
rmarkdown::render(input = "README.Rmd")
devtools::build_readme()
usethis::use_news_md()

## CREATE WEBSITE WITH EVERY DOCUMENTATION/VIGNETTES
## INIT WEBSITE
usethis::use_pkgdown()

##
pkgdown::clean_site(force = TRUE)
pkgdown::clean_cache()
devtools::document(); devtools::load_all()

## RENDER ALL VIGNETTE FROM VIGNETTES INTO ARTICLES IN DOCS
pkgdown::build_articles(lazy = TRUE)

## IMPORT DOC FROM FUNCTIONS/DATA INTO REFERENCES
pkgdown::build_reference()

## PREVIEW SITE LOCALLY
pkgdown::build_site()


## INIT GITHUB
gitcreds::gitcreds_get()
gitcreds::gitcreds_set()


## CREATE WEBPAGE ON GITHUB
## PROJECT NEED TO BE PUBLIC
usethis::gh_token_help()
usethis::use_pkgdown_github_pages()

## TIDY DESCRIPTION FILE
usethis::use_tidy_description()

## DATA
usethis::use_data_raw()

## ROXYGEN TO DOCUMENT
## LOAD COMPLETE PACKAGE
devtools::document(); devtools::load_all()
## BUILD and CHECK PACKAGE (devtools::check("../"))
devtools::check()
devtools::build(); devtools::check_built("../MOTL_0.99.0.tar.gz")
devtools::build(vignettes = FALSE); devtools::check_built("../MOTL_0.99.0.tar.gz")





## PERFORM TESTS
devtools::document(); devtools::load_all()
devtools::test()

## VIGNETTES
devtools::document(); devtools::load_all()
devtools::install(build_vignettes = TRUE)
#vignette(package = "MOTL")
??MOTL

## RENDER VIGNETTE HTML
devtools::document(); devtools::load_all()
devtools::build_rmd("vignettes/MOTL.Rmd")
devtools::build_rmd("vignettes/UseCase_Glioblastoma.Rmd")

## RENDER README
devtools::document(); devtools::load_all()
devtools::build_readme()

## ROXYGEN DOCUMENT UPDATE AND LOAD COMPLETE PACKAGE
## BUILD and CHECK PACKAGE
## BIOCONDUCTOR CHECKING
devtools::document(); devtools::load_all()
devtools::build(); devtools::check_built("../MOTL_0.99.0.tar.gz")
BiocCheck::BiocCheck("../MOTL_0.99.0.tar.gz")
devtools::build(vignettes = TRUE, manual = TRUE)
devtools::install(build_vignettes = TRUE)
