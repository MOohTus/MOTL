## MT - 16/09/2025
##

## LEARNING DATA SET
Lrn_FctrnDir <- base::system.file("tests/testthat/fixtures/Lrn_4test_5000D/", "", package = "MOTL")
Lrn_ModelFile <- base::system.file("tests/testthat/fixtures/Lrn_4test_5000D/Fctrzn_50K_01TH", "Model.hdf5", package = "MOTL")

Lrn_meta <- readRDS(test_path("fixtures/Lrn_4test_5000D", "Lrn_meta.rds"))
Lrn_Fctrzn = load_model(file = Lrn_ModelFile)

## TARGET DATA SET
Trg <- readRDS(test_path("fixtures", "Trg_4test_4omics.rds"))
YTrg_list <- Trg$YTrg_list
YTrg_meta <- Trg$Trg_meta
Trg <- readRDS(test_path("fixtures", "Trg_4test_4omics.rds"))

## PARAMETERS
views <- c("mRNA", "miRNA", "DNAme", "SNV")
smpls <- colnames(YTrg_list$mRNA)

## SELECT RANDOM FEATURE NAMES
## REMOVE FEATURE NAMES VERSION
## FOR mRNA_addVersion TESTING
ftrs <- sample(x = Lrn_meta$ftrs_mRNA, size = 100)
ftrs_wov <- unlist(lapply(strsplit(ftrs, "[.]"), function(x){return(x[[1]])}))
df <- data.frame(row.names = ftrs_wov, "sample1" = seq(1:100))
Lrndat <- data.frame(row.names = Lrn_meta$ftrs_mRNA, "view" = rep("mRNA", length(Lrn_meta$ftrs_mRNA)))




# Lrn_Fctrn <- readRDS(test_path("fixtures/Lrn_4test_5000D", "Lrn_5000D_Fctrzn_model.rds"))
# Fctrzn <- readRDS(test_path("fixtures","Lrn_5000D_Fctrzn_model.rds"))
#
#
#
# InputModel <- base::system.file("inst/extdata", "Lrn_5000D_Fctrzn_model.hdf5", package = "MOTL")
# LrnFctrnDir <- base::system.file("inst/extdata", "", package = "MOTL")
#
# ## IMPORT LEARNING DATA
# expdat_meta_Lrn <- readRDS(test_path("fixtures","expdat_meta_Lrn.rds"))
# Fctrzn <- readRDS(test_path("fixtures","Lrn_5000D_Fctrzn_model.rds"))
#
# ## INIT VALUES FROM FACTORIZATION
# viewsLrn <- Fctrzn@data_options$views
# likelihoodsLrn <- Fctrzn@model_options$likelihoods
# MLrn <- Fctrzn@dimensions$M
#
# ## SELECT RANDOM FEATURE NAMES
# ## REMOVE FEATURE NAMES VERSION
# ## FOR mRNA_addVersion TESTING
# ftrs <- sample(x = expdat_meta_Lrn$ftrs_mRNA, size = 100)
# ftrs_wov <- unlist(lapply(strsplit(ftrs, "[.]"), function(x){return(x[[1]])}))
# df <- data.frame(row.names = ftrs_wov, "sample1" = seq(1:100))
# Lrndat <- data.frame(row.names = expdat_meta_Lrn$ftrs_mRNA, "view" = rep("mRNA", length(expdat_meta_Lrn$ftrs_mRNA)))
#
# ## IMPORT TARGET DATA
# Trg <- readRDS(test_path("fixtures", "Trg_4test_4omics.rds"))
# YTrg_list <- Trg$YTrg_list
#
# ## INIT PARAMETER FOR TRANSFER LEARNING
# CenterTrg <- FALSE
# smpls <- colnames(YTrg_list[[1]])
# viewsTrg <- names(YTrg_list)
# views <- viewsLrn[is.element(viewsLrn, viewsTrg)]
# likelihoods = likelihoodsLrn[views]
#
#
# ##
# TL_param <- readRDS(test_path("fixtures", "TL_param.rds"))
# YTrg <- TL_param$YTrg
# Fctrzn_Lrn_W0 <- TL_param$Fctrzn_Lrn_W0
# Fctrzn_Lrn_W <- TL_param$Fctrzn_Lrn_W
# Fctrzn_Lrn_WSq <- TL_param$Fctrzn_Lrn_WSq
# Tau <- TL_param$Tau
# TauLn <- TL_param$TauLn
# ZMu <- TL_param$ZMu
# ZVar <- TL_param$ZVar
# ZMu_0 <- TL_param$ZMu_0
# ZMuSq <- TL_param$ZMuSq
#
# PoisRateCstnt = 0.0001
