## MT - 16/09/2025
##

## DIRECTORIES FOR TEST FILES
Lrn_ModelDir <- base::system.file("tests/testthat/fixtures/", package = "MOTL")

## LEARNING DATASET METADATA
## LEARNING DATASET FACTORIZATION MODEL
## LEARNING DATASET FACTORIZATION INITIALIZED BEFORE TRANSFER LEARNING
#data(Lrn, package = "MOTL")
Lrn <- MOTL::Lrn
Lrn_meta <- Lrn$Lrn_meta
Lrn_Fctrzn <- Lrn$Fctrzn
Lrn_Fctrzn_init <- Lrn$Fctrzn_init

## TARGET DATASET METADATA
## TARGET DATASET LIST
## TARGET DATASET LIST PREPARED FOR TRANSFER LEARNING
#data(Trg, package = "MOTL")
Trg <- MOTL::Trg
Trg_meta <- Trg$Trg_meta
YTrg_list <- Trg$YTrg_list
YTrg_prep <- Trg$YTrg_prep

## PARAMETERS FOR TRANSFER LEARNING
CenterTrg <- FALSE
views <- c("mRNA" = "mRNA", "miRNA" = "miRNA", "DNAme" = "DNAme", "SNV" = "SNV")
smpls <- colnames(YTrg_list$mRNA)
likelihoods <- Lrn_Fctrzn@model_options$likelihoods
PoisRateCstnt <- 0.0001

## TRANSFER LEARNING
##
#data(TL_param, package = "MOTL")
TL_param <- MOTL::TL_param
YTrg <- TL_param$YTrg
Fctrzn_Lrn_W0 <- TL_param$Fctrzn_Lrn_W0
Fctrzn_Lrn_W <- TL_param$Fctrzn_Lrn_W
Fctrzn_Lrn_WSq <- TL_param$Fctrzn_Lrn_WSq
ZMu_0 <- TL_param$ZMu_0
ZMu <- TL_param$ZMu
ZMuSq <- TL_param$ZMuSq
Tau <- TL_param$Tau
TauLn <- TL_param$TauLn

## EXPECTATION CALCULATION FOR EACH OMICS
E_ZE_W <- sapply(views, E_ZE_W_update, ZMu_0, ZMu, Fctrzn_Lrn_W0, Fctrzn_Lrn_W)
E_Z_SqE_W_Sq <- sapply(views, E_Z_SqE_W_Sq_update, ZMu_0, ZMu, Fctrzn_Lrn_W0, Fctrzn_Lrn_W)
E_ZSqE_WSq <- sapply(views, E_ZSqE_WSq_update, ZMu_0, ZMuSq, Fctrzn_Lrn_W0, Fctrzn_Lrn_WSq)
E_ZWSq <- sapply(views, E_ZWSq_update, E_ZE_W, ZMuSq, E_Z_SqE_W_Sq, E_ZSqE_WSq)

## ZETA CALCULATION FOR EACH OMICS
Zeta <- sapply(views, Zeta_calculation, likelihoods, E_ZWSq, E_ZE_W)

## TAU CALCULATION FOR EACH OMICS
Tau <- sapply(views, Tau_calculation, likelihoods, Zeta, Tau)

## SELECT RANDOM FEATURE NAMES
## REMOVE FEATURE NAMES VERSION
## FOR mRNA_addVersion TESTING
ftrs <- sample(x = Lrn_meta$ftrs_mRNA, size = 100)
ftrs_wov <- unlist(lapply(strsplit(ftrs, "[.]"), function(x) {
    return(x[[1]])
}))
df <- data.frame(row.names = ftrs_wov, "sample1" = seq(1:100))
Lrndat <- data.frame(row.names = Lrn_meta$ftrs_mRNA, "view" = rep("mRNA", length(Lrn_meta$ftrs_mRNA)))

