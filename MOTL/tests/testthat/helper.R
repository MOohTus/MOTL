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








#
# Fctrzn@expectations[["Tau"]] = Tau_init(viewsLrn, Fctrzn, InputModel)
# Fctrzn@expectations[["TauLn"]] = sapply(viewsLrn, TauLn_calculation, likelihoodsLrn, Fctrzn, LrnFctrnDir)
# Fctrzn@expectations[["WSq"]] = sapply(viewsLrn, WSq_calculation, Fctrzn, LrnFctrnDir)
# Fctrzn@expectations[["W0"]] = sapply(viewsLrn, W0_calculation, CenterTrg, Fctrzn, LrnFctrnDir)
