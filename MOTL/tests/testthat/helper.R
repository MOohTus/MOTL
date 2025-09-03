##
## INIT PARAMTER
InputModel <- base::system.file("inst/extdata", "Lrn_5000D_Fctrzn_model.hdf5", package = "MOTL")
LrnFctrnDir <- base::system.file("inst/extdata", "", package = "MOTL")
CenterTrg <- FALSE

## IMPORT LEARNING DATA
expdat_meta_Lrn <- readRDS(test_path("fixtures","expdat_meta_Lrn.rds"))
Fctrzn <- readRDS(test_path("fixtures","Lrn_5000D_Fctrzn_model.rds"))

## INIT VALUES FROM FACTORIZATION
viewsLrn <- Fctrzn@data_options$views
likelihoodsLrn <- Fctrzn@model_options$likelihoods
MLrn <- Fctrzn@dimensions$M

## SELECT RANDOM FEATURE NAMES
## REMOVE FEATURE NAMES VERSION
## FOR mRNA_addVersion TESTING
ftrs <- sample(x = expdat_meta_Lrn$ftrs_mRNA, size = 100)
ftrs_wov <- unlist(lapply(strsplit(ftrs, "[.]"), function(x){return(x[[1]])}))
df <- data.frame(row.names = ftrs_wov, "sample1" = seq(1:100))
Lrndat = data.frame(row.names = expdat_meta_Lrn$ftrs_mRNA, "view" = rep("mRNA", length(expdat_meta_Lrn$ftrs_mRNA)))
