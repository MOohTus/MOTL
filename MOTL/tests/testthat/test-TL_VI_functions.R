skip_issue_2_solve <- function(message) {
    skip(paste0("ISSUE 2 SOLVE - ASK_DAVID - ", message))
}

test_that("TCGATargetDataPrefiltering", {
  print("to do")
})

test_that("TCGATargetDataPreparation", {
  print("to do")
})

test_that("mRNA_addVersion", {
    expdat_mRNA_wov_up = mRNA_addVersion(expdat = df, Lrndat)
    expect_in(rownames(expdat_mRNA_wov_up), rownames(Lrndat))
})

test_that("GeoMeanFun", {
  x <- YTrg_list$mRNA[41,]
  GeoMean <- GeoMeanFun(x)
  expect_equal(length(GeoMean), 1)
  expect_equal(GeoMean, exp(sum(log(x[x > 0]))/length(x)))
})


test_that("GeoMeans_Lrn_init", {
  YTrg_prep = sapply(views, TargetDataPrefiltering, YTrg_list, Fctrzn, smpls)
  GeoMeans <- GeoMeans_Lrn_init("mRNA", expdat_meta_Lrn, rownames(YTrg_prep$mRNA))
  expect_equal(length(GeoMeans), length(rownames(YTrg_prep$mRNA)))
  expect_in(GeoMeans, expdat_meta_Lrn$GeoMeans_mRNA)
  GeoMeans <- GeoMeans_Lrn_init("miRNA", expdat_meta_Lrn, rownames(YTrg_prep$miRNA))
  expect_equal(length(GeoMeans), length(rownames(YTrg_prep$miRNA)))
  expect_in(GeoMeans, expdat_meta_Lrn$GeoMeans_miRNA)
  GeoMeans <- GeoMeans_Lrn_init("DNAme", expdat_meta_Lrn, rownames(YTrg_prep$DNAme))
  expect_equal(GeoMeans, numeric())
})


test_that("countsTransformation", {
  expdat_trans <- countsTransformation(YTrg_list$mRNA, 50)
  expect_equal(nrow(expdat_trans), 50)
  expect_equal(names(expdat_trans), names(YTrg_list$mRNA))
  expect_in(rownames(expdat_trans), rownames(YTrg_list$mRNA))
})

test_that("countsNormalization_Lrn", {
  expdat <- TargetDataPrefiltering(view = "mRNA", YTrg_list, Fctrzn, smpls)
  expdat <- apply(expdat, MARGIN = 2, round)
  expdat <- SummarizedExperiment(assays = expdat)
  expdat_norm <- countsNormalization(expdat, GeoMeans = "Lrn")
  expect_equal(length(expdat_norm), 2)
  expect_equal(dim(expdat_norm$counts), dim(expdat))
  expect_equal(rownames(expdat_norm$counts), rownames(expdat))
  expect_equal(colnames(expdat_norm$counts), colnames(expdat))
})

test_that("countsNormalization_Trg", {
  expdat = TargetDataPrefiltering(view = "mRNA", YTrg_list, Fctrzn, smpls)
  expdat <- apply(expdat, MARGIN = 2, round)
  expdat <- SummarizedExperiment(assays = expdat)
  expdat_norm <- countsNormalization(expdat, GeoMeans = "Trg")
  expect_equal(length(expdat_norm), 1)
  expect_equal(dim(expdat_norm$counts), dim(expdat))
  expect_equal(rownames(expdat_norm$counts), rownames(expdat))
  expect_equal(colnames(expdat_norm$counts), colnames(expdat))
})

test_that("countsNormalization_num", {
  expdat = TargetDataPrefiltering(view = "mRNA", YTrg_list, Fctrzn, smpls)
  expdat <- apply(expdat, MARGIN = 2, round)
  expdat <- SummarizedExperiment(assays = expdat)
  GeoMeans <- GeoMeans_Lrn_init("mRNA", expdat_meta_Lrn, rownames(expdat))
  expdat_norm <- countsNormalization(expdat, GeoMeans = GeoMeans)
  expect_equal(length(expdat_norm), 1)
  expect_equal(dim(expdat_norm$counts), dim(expdat))
  expect_equal(rownames(expdat_norm$counts), rownames(expdat))
  expect_equal(colnames(expdat_norm$counts), colnames(expdat))
})


test_that("TargetDataPrefiltering", {
  YTrg_prep = sapply(views, TargetDataPrefiltering, YTrg_list, Fctrzn, smpls)
  expect_equal(names(YTrg_prep), names(YTrg_list))
  expect_equal(nrow(YTrg_prep$mRNA), 3043)
  expect_equal(nrow(YTrg_prep$miRNA), 1047)
  expect_equal(nrow(YTrg_prep$DNAme), 2180)
  expect_equal(colnames(YTrg_prep$mRNA), colnames(YTrg_prep$miRNA))
  expect_equal(colnames(YTrg_prep$mRNA), colnames(YTrg_prep$DNAme))
  expect_equal(colnames(YTrg_prep$mRNA), smpls)
  expect_in(rownames(YTrg_prep$mRNA), Fctrzn@features_metadata[Fctrzn@features_metadata$view == "mRNA",1])
  expect_in(rownames(YTrg_prep$miRNA), Fctrzn@features_metadata[Fctrzn@features_metadata$view == "miRNA",1])
  expect_in(rownames(YTrg_prep$DNAme), Fctrzn@features_metadata[Fctrzn@features_metadata$view == "DNAme",1])
})

test_that("preprocessCountsData_SIMPLE", {
  YTrg_pre = sapply(views, TargetDataPrefiltering, YTrg_list, Fctrzn, smpls)
  ## NO TRANSFORMATION NO NORMALIZATION
  YTrg = sapply(viewsTrg, preprocessCountsData, YTrg_pre, normalization = FALSE, expdat_meta_Lrn, transformation = FALSE)
  expect_equal(names(YTrg), names(YTrg_pre))
  expect_equal(YTrg$mRNA, YTrg_pre$mRNA)
  expect_equal(YTrg$miRNA, YTrg_pre$miRNA)
  expect_equal(YTrg$DNAme, YTrg_pre$DNAme)
})

test_that("preprocessCountsData_TRANS", {
  YTrg_pre = sapply(views, TargetDataPrefiltering, YTrg_list, Fctrzn, smpls)
  ## TRANSFORMATION NO NORMALIZATION
  YTrg = sapply(viewsTrg, preprocessCountsData, YTrg_pre, normalization = FALSE, expdat_meta_Lrn, transformation = TRUE)
  expect_equal(names(YTrg), names(YTrg_pre))
  expect_equal(names(YTrg$mRNA), names(YTrg_pre$mRNA))
  expect_equal(rownames(YTrg$mRNA), rownames(YTrg_pre$mRNA))
  expect_equal(names(YTrg$miRNA), names(YTrg_pre$miRNA))
  expect_equal(rownames(YTrg$miRNA), rownames(YTrg_pre$miRNA))
  expect_equal(names(YTrg$DNAme), names(YTrg_pre$DNAme))
  expect_equal(rownames(YTrg$DNAme), rownames(YTrg_pre$DNAme))
})

test_that("preprocessCountsData_NORM_LRN", {
  ## NO TRANSFORMATION NORMALIZATION WITH LRN DATASET
  expdat <- list()
  expdat$mRNA <- TargetDataPrefiltering(view = "mRNA", YTrg_list, Fctrzn, smpls)
  expdat$mRNA <- apply(expdat$mRNA, MARGIN = 2, round)
  expdat_norm <- preprocessCountsData(view = "mRNA", expdat, normalization = "Lrn", expdat_meta_Lrn, transformation = FALSE)
  expect_equal(dim(expdat_norm), dim(expdat$mRNA))
  expect_equal(rownames(expdat_norm), rownames(expdat$mRNA))
  expect_equal(colnames(expdat_norm), colnames(expdat$mRNA))
})


test_that("preprocessCountsData_NORM_TRG", {
  ## NO TRANSFORMATION NORMALIZATION WITH TRG DATASET
  expdat <- list()
  expdat$mRNA <- TargetDataPrefiltering(view = "mRNA", YTrg_list, Fctrzn, smpls)
  expdat$mRNA <- apply(expdat$mRNA, MARGIN = 2, round)
  expdat_norm <- preprocessCountsData(view = "mRNA", expdat, normalization = "Trg", expdat_meta_Lrn, transformation = FALSE)
  expect_equal(dim(expdat_norm), dim(expdat$mRNA))
  expect_equal(rownames(expdat_norm), rownames(expdat$mRNA))
  expect_equal(colnames(expdat_norm), colnames(expdat$mRNA))
})


test_that("TargetDataPreparation", {
  expdat_prep <- TargetDataPreparation(views = c("mRNA", "miRNA"), YTrg_list, Fctrzn, smpls, expdat_meta_Lrn, normalization = FALSE, transformation = FALSE)
  expect_equal(names(expdat_prep), c("mRNA", "miRNA"))
  expect_equal(colnames(expdat_prep$mRNA), colnames(expdat_prep$miRNA))
})

test_that("initTransferLearningParamaters", {
  ## INPUT PREPARATION
  Fctrzn@expectations[["Tau"]] = Tau_init(viewsLrn, Fctrzn, InputModel)
  Fctrzn@expectations[["TauLn"]] = sapply(viewsLrn, TauLn_calculation, likelihoodsLrn, Fctrzn, LrnFctrnDir)
  Fctrzn@expectations[["WSq"]] = sapply(viewsLrn, WSq_calculation, Fctrzn, LrnFctrnDir)
  Fctrzn@expectations[["W0"]] = sapply(viewsLrn, W0_calculation, CenterTrg, Fctrzn, LrnFctrnDir)
  YTrg_prep = sapply(views, TargetDataPrefiltering, YTrg_list, Fctrzn, smpls)
  YTrgFtrs <- lapply(YTrg_prep, rownames)
  ##
  TL_param = initTransferLearningParamaters(YTrg_prep, views, expdat_meta_Lrn, Fctrzn, likelihoods)
  expect_equal(names(TL_param), c("YTrg", "Fctrzn_Lrn_W0", "Fctrzn_Lrn_W", "Fctrzn_Lrn_WSq", "Tau", "TauLn"))
  expect_equal(names(TL_param$Fctrzn_Lrn_W0$mRNA), YTrgFtrs$mRNA)
  expect_equal(names(TL_param$Fctrzn_Lrn_W0$miRNA), YTrgFtrs$miRNA)
  expect_equal(names(TL_param$Fctrzn_Lrn_W0$DNAme), YTrgFtrs$DNAme)
  expect_equal(rownames(TL_param$Fctrzn_Lrn_W$mRNA), YTrgFtrs$mRNA)
  expect_equal(rownames(TL_param$Fctrzn_Lrn_W$miRNA), YTrgFtrs$miRNA)
  expect_equal(rownames(TL_param$Fctrzn_Lrn_W$DNAme), YTrgFtrs$DNAme)
  expect_equal(rownames(TL_param$Fctrzn_Lrn_WSq$mRNA), YTrgFtrs$mRNA)
  expect_equal(rownames(TL_param$Fctrzn_Lrn_WSq$miRNA), YTrgFtrs$miRNA)
  expect_equal(rownames(TL_param$Fctrzn_Lrn_WSq$DNAme), YTrgFtrs$DNAme)
  expect_equal(nrow(TL_param$Tau$mRNA), ncol(YTrg_prep$mRNA))
  expect_equal(ncol(TL_param$Tau$mRNA), nrow(YTrg_prep$mRNA))
  expect_equal(nrow(TL_param$Tau$miRNA), ncol(YTrg_prep$miRNA))
  expect_equal(ncol(TL_param$Tau$miRNA), nrow(YTrg_prep$miRNA))
  expect_equal(nrow(TL_param$Tau$DNAme), ncol(YTrg_prep$DNAme))
  expect_equal(ncol(TL_param$Tau$DNAme), nrow(YTrg_prep$DNAme))
  ## MAKE TEST WITH SNV - ASK_DAVID
  expect_equal(nrow(TL_param$TauLn$mRNA), ncol(YTrg_prep$mRNA))
  expect_equal(ncol(TL_param$TauLn$mRNA), nrow(YTrg_prep$mRNA))
  expect_equal(nrow(TL_param$TauLn$miRNA), ncol(YTrg_prep$miRNA))
  expect_equal(ncol(TL_param$TauLn$miRNA), nrow(YTrg_prep$miRNA))
  expect_equal(nrow(TL_param$TauLn$DNAme), ncol(YTrg_prep$DNAme))
  expect_equal(ncol(TL_param$TauLn$DNAme), nrow(YTrg_prep$DNAme))
  expect_equal(TL_param$YTrg$mRNA, t(YTrg_prep$mRNA))
  expect_equal(TL_param$YTrg$miRNA, t(YTrg_prep$miRNA))
  expect_equal(TL_param$YTrg$DNAme, t(YTrg_prep$DNAme))
})

test_that("Tau_init", {
  Fctrzn@expectations[["Tau"]] <- Tau_init(viewsLrn, Fctrzn, InputModel)
  expect_equal(nrow(Fctrzn@expectations[["Tau"]]$mRNA$group0), Fctrzn@dimensions$D["mRNA"][[1]])
  expect_equal(nrow(Fctrzn@expectations[["Tau"]]$miRNA$group0), Fctrzn@dimensions$D["miRNA"][[1]])
  expect_equal(nrow(Fctrzn@expectations[["Tau"]]$DNAme$group0), Fctrzn@dimensions$D["DNAme"][[1]])
  expect_equal(nrow(Fctrzn@expectations[["Tau"]]$SNV$group0), Fctrzn@dimensions$D["SNV"][[1]])
})

##
test_that("TauLn_calculation_TRUE", {
  Fctrzn@expectations[["Tau"]] <- Tau_init(viewsLrn, Fctrzn, InputModel)
  ## LrnSimple = TRUE / DON'T USE LrnFctrnDir TAU FILES
  Fctrzn@expectations[["TauLn"]] <-  sapply(viewsLrn, TauLn_calculation, likelihoodsLrn, Fctrzn, LrnFctrnTauLn, LrnSimple = TRUE)
  expect_equal(length(Fctrzn@expectations[["TauLn"]]$mRNA), Fctrzn@dimensions$D["mRNA"][[1]])
  expect_equal(Fctrzn@expectations[["TauLn"]]$mRNA, log(Fctrzn@expectations$Tau$mRNA$group0[,1]))
  expect_equal(length(Fctrzn@expectations[["TauLn"]]$miRNA), Fctrzn@dimensions$D["miRNA"][[1]])
  expect_equal(Fctrzn@expectations[["TauLn"]]$miRNA, log(Fctrzn@expectations$Tau$miRNA$group0[,1]))
  expect_equal(length(Fctrzn@expectations[["TauLn"]]$DNAme), Fctrzn@dimensions$D["DNAme"][[1]])
  expect_equal(Fctrzn@expectations[["TauLn"]]$DNAme, log(Fctrzn@expectations$Tau$DNAme$group0[,1]))
  expect_equal(Fctrzn@expectations[["TauLn"]]$SNV, numeric())
})

test_that("TauLn_calculation_FALSE", {
  skip_issue_2_solve("miRNA length is different between model and imported TauLn file")
  ## LrnSimple = FALSE / USE LrnFctrnDir TAU FILES
  Fctrzn@expectations[["TauLn"]]$mRNA <- TauLn_calculation(view = "mRNA", likelihoodsLrn, Fctrzn, LrnFctrnDir, LrnSimple = FALSE)
  expect_equal(length(Fctrzn@expectations[["TauLn"]]$mRNA), Fctrzn@dimensions$D["mRNA"][[1]])
  Fctrzn@expectations[["TauLn"]]$miRNA <- TauLn_calculation(view = "miRNA", likelihoodsLrn, Fctrzn, LrnFctrnDir, LrnSimple = FALSE)
  expect_equal(length(Fctrzn@expectations[["TauLn"]]$miRNA), Fctrzn@dimensions$D["miRNA"][[1]])
  Fctrzn@expectations[["TauLn"]]$DNAme <- TauLn_calculation(view = "DNAme", likelihoodsLrn, Fctrzn, LrnFctrnDir, LrnSimple = FALSE)
  expect_equal(length(Fctrzn@expectations[["TauLn"]]$DNAme), Fctrzn@dimensions$D["DNAme"][[1]])
  Fctrzn@expectations[["TauLn"]]$SNV <- TauLn_calculation(view = "SNV", likelihoodsLrn, Fctrzn, LrnFctrnDir, LrnSimple = FALSE)
  expect_equal(Fctrzn@expectations[["TauLn"]]$SNV, numeric())
})

##
test_that("WSq_calculation_TRUE", {
  ## LrnSimple = TRUE / DON'T USE LrnFctrnDir W FILES
  Fctrzn@expectations[["WSq"]] <- sapply(viewsLrn, WSq_calculation, Fctrzn, LrnFctrnDir, LrnSimple = TRUE)
  expect_equal(nrow(Fctrzn@expectations[["WSq"]]$mRNA), Fctrzn@dimensions$D["mRNA"][[1]])
  expect_equal(Fctrzn@expectations[["WSq"]]$mRNA, Fctrzn@expectations$W[["mRNA"]]^2)
  expect_equal(nrow(Fctrzn@expectations[["WSq"]]$miRNA), Fctrzn@dimensions$D["miRNA"][[1]])
  expect_equal(Fctrzn@expectations[["WSq"]]$miRNA, Fctrzn@expectations$W[["miRNA"]]^2)
  expect_equal(nrow(Fctrzn@expectations[["WSq"]]$DNAme), Fctrzn@dimensions$D["DNAme"][[1]])
  expect_equal(Fctrzn@expectations[["WSq"]]$DNAme, Fctrzn@expectations$W[["DNAme"]]^2)
  expect_equal(nrow(Fctrzn@expectations[["WSq"]]$SNV), Fctrzn@dimensions$D["SNV"][[1]])
  expect_equal(Fctrzn@expectations[["WSq"]]$SNV, Fctrzn@expectations$W[["SNV"]]^2)
})

test_that("WSq_calculation_FALSE", {
  skip_issue_2_solve("W files contain 84 samples and model contains 91 samples ... wrong files")
  ## LrnSimple = FALSE / USE LrnFctrnDir W FILES
  Fctrzn@expectations[["WSq"]] = sapply(viewsLrn, WSq_calculation, Fctrzn, LrnFctrnDir, LrnSimple = FALSE)
  Fctrzn@expectations[["WSq"]]$mRNA <- WSq_calculation(view = "mRNA", Fctrzn, LrnFctrnDir, LrnSimple = FALSE)
  expect_equal(length(Fctrzn@expectations[["WSq"]]$mRNA), Fctrzn@dimensions$D["mRNA"][[1]])
  Fctrzn@expectations[["WSq"]]$miRNA <- WSq_calculation(view = "miRNA", Fctrzn, LrnFctrnDir, LrnSimple = FALSE)
  expect_equal(length(Fctrzn@expectations[["WSq"]]$miRNA), Fctrzn@dimensions$D["miRNA"][[1]])
  Fctrzn@expectations[["WSq"]]$DNAme <- WSq_calculation(view = "DNAme", Fctrzn, LrnFctrnDir, LrnSimple = FALSE)
  expect_equal(length(Fctrzn@expectations[["WSq"]]$DNAme), Fctrzn@dimensions$D["DNAme"][[1]])
  Fctrzn@expectations[["WSq"]]$SNV <- WSq_calculation(view = "SNV", Fctrzn, LrnFctrnDir, LrnSimple = FALSE)
  expect_equal(length(Fctrzn@expectations[["WSq"]]$SNV), Fctrzn@dimensions$D["SNV"][[1]])
})

test_that("W0_calculation_FALSE", {
  ## CenterTrg = FALSE
  Fctrzn@expectations[["W0"]] = sapply(viewsLrn, W0_calculation, CenterTrg = FALSE, Fctrzn, LrnFctrnDir)
  expect_equal(length(Fctrzn@expectations[["W0"]]$mRNA), Fctrzn@dimensions$D["mRNA"][[1]])
  expect_equal(length(Fctrzn@expectations[["W0"]]$miRNA), Fctrzn@dimensions$D["miRNA"][[1]])
  expect_equal(length(Fctrzn@expectations[["W0"]]$DNAme), Fctrzn@dimensions$D["DNAme"][[1]])
  expect_equal(length(Fctrzn@expectations[["W0"]]$SNV), Fctrzn@dimensions$D["SNV"][[1]])
})

test_that("W0_calculation_TRUE", {
  ## CenterTrg = TRUE
  Fctrzn@expectations[["W0"]] = sapply(viewsLrn, W0_calculation, CenterTrg = TRUE, Fctrzn, LrnFctrnDir)
  expect_equal(length(Fctrzn@expectations[["W0"]]$mRNA), Fctrzn@dimensions$D["mRNA"][[1]])
  expect_equal(sum(Fctrzn@expectations[["W0"]]$mRNA), 0)
  expect_equal(length(Fctrzn@expectations[["W0"]]$miRNA), Fctrzn@dimensions$D["miRNA"][[1]])
  expect_equal(sum(Fctrzn@expectations[["W0"]]$miRNA), 0)
  expect_equal(length(Fctrzn@expectations[["W0"]]$DNAme), Fctrzn@dimensions$D["DNAme"][[1]])
  expect_equal(sum(Fctrzn@expectations[["W0"]]$DNAme), 0)
  expect_equal(length(Fctrzn@expectations[["W0"]]$SNV), Fctrzn@dimensions$D["SNV"][[1]])
  expect_equal(sum(Fctrzn@expectations[["W0"]]$SNV), 0)
})

test_that("intercepts_calculation", {
  print("to do")
})

test_that("Zeta_calculation", {
  likelihoods <- list("mRNA" = "gaussian", "miRNA" = "bernoulli")
  E_ZE_W <- list("mRNA" = E_ZE_W_update(view = "mRNA", ZMu_0, ZMu, Fctrzn_Lrn_W0, Fctrzn_Lrn_W))
  E_Z_SqE_W_Sq <- list("mRNA" = E_Z_SqE_W_Sq_update(view = "mRNA", ZMu_0, ZMu, Fctrzn_Lrn_W0, Fctrzn_Lrn_W))
  E_ZSqE_WSq <- list("mRNA" = E_ZSqE_WSq_update(view = "mRNA", ZMu_0, ZMuSq, Fctrzn_Lrn_W0, Fctrzn_Lrn_WSq))
  E_ZWSq <- list("mRNA" = E_ZWSq_update(view = "mRNA", E_ZE_W, ZMuSq, E_Z_SqE_W_Sq, E_ZSqE_WSq))
  E_ZWSq$miRNA <- E_ZWSq$mRNA
  E_ZE_W$miRNA <- E_ZE_W$mRNA
  zeta <- Zeta_calculation(view = "mRNA", likelihoods, E_ZWSq, E_ZE_W)
  expect_equal(zeta, E_ZE_W$mRNA)
  zeta <- Zeta_calculation(view = "miRNA", likelihoods, E_ZWSq, E_ZE_W)
  expect_equal(zeta, sqrt(E_ZWSq$miRNA))
})

test_that("Tau_calculation", {
  likelihoods <- list("mRNA" = "gaussian", "miRNA" = "bernoulli")
  E_ZE_W <- list("mRNA" = E_ZE_W_update(view = "mRNA", ZMu_0, ZMu, Fctrzn_Lrn_W0, Fctrzn_Lrn_W))
  E_Z_SqE_W_Sq <- list("mRNA" = E_Z_SqE_W_Sq_update(view = "mRNA", ZMu_0, ZMu, Fctrzn_Lrn_W0, Fctrzn_Lrn_W))
  E_ZSqE_WSq <- list("mRNA" = E_ZSqE_WSq_update(view = "mRNA", ZMu_0, ZMuSq, Fctrzn_Lrn_W0, Fctrzn_Lrn_WSq))
  E_ZWSq <- list("mRNA" = E_ZWSq_update(view = "mRNA", E_ZE_W, ZMuSq, E_Z_SqE_W_Sq, E_ZSqE_WSq))
  E_ZWSq$miRNA <- E_ZWSq$mRNA
  E_ZE_W$miRNA <- E_ZE_W$mRNA
  Zeta <- list("mRNA" = Zeta_calculation(view = "mRNA", likelihoods, E_ZWSq, E_ZE_W))
  Zeta$miRNA <- Zeta$mRNA
  Tau_m <- Tau_calculation(view = "mRNA", likelihoods, Zeta, Tau)
  expect_equal(Tau_m, Tau$mRNA)
  Tau_m <- Tau_calculation(view = "miRNA", likelihoods, Zeta, Tau)
  expect_equal(Tau_m, (1/2)*(1/Zeta$miRNA)*tanh(Zeta$miRNA/2))
})

test_that("YGauss_calculation", {
  likelihoods <- list("mRNA" = "gaussian", "miRNA" = "bernoulli")
  E_ZE_W <- sapply(c("mRNA", "miRNA"), E_ZE_W_update, ZMu_0, ZMu, Fctrzn_Lrn_W0, Fctrzn_Lrn_W)
  E_Z_SqE_W_Sq <- sapply(c("mRNA", "miRNA"), E_Z_SqE_W_Sq_update, ZMu_0, ZMu, Fctrzn_Lrn_W0, Fctrzn_Lrn_W)
  E_ZSqE_WSq <- sapply(c("mRNA", "miRNA"), E_ZSqE_WSq_update, ZMu_0, ZMuSq, Fctrzn_Lrn_W0, Fctrzn_Lrn_WSq)
  E_ZWSq <- sapply(c("mRNA", "miRNA"), E_ZWSq_update, E_ZE_W, ZMuSq, E_Z_SqE_W_Sq, E_ZSqE_WSq)
  Zeta <- sapply(c("mRNA", "miRNA"), Zeta_calculation, likelihoods, E_ZWSq, E_ZE_W)
  Tau_m <- sapply(c("mRNA", "miRNA"), Tau_calculation, likelihoods, Zeta, Tau)
  YGauss <- YGauss_calculation("mRNA", likelihoods, YTrg, Zeta, Tau_m, CenterTrg, PoisRateCstnt)
  expect_equal(YGauss, YTrg$mRNA)
  YGauss <- YGauss_calculation("miRNA", likelihoods, YTrg, Zeta, Tau_m, CenterTrg, PoisRateCstnt)
  expect_equal(YGauss, (2*YTrg$miRNA - 1) / (2*Tau_m$miRNA))
})

test_that("ZVar_calculation", {
  # ASK_DAVID - matrix looks strange - number inside are same
  ZVar <- ZVar_calculation(view = "mRNA", Tau = TL_param$Tau, Fctrzn_Lrn_WSq = TL_param$Fctrzn_Lrn_WSq)
  expect_equal(nrow(ZVar), nrow(TL_param$Tau$mRNA))
  expect_equal(ncol(ZVar), ncol(TL_param$Fctrzn_Lrn_WSq$mRNA))
  expect_equal(ZVar, TL_param$Tau$mRNA %*% TL_param$Fctrzn_Lrn_WSq$mRNA)
})

test_that("ZMu_calculation", {
  E_ZE_W <- list("mRNA" = E_ZE_W_update(view = c("mRNA"), ZMu_0, ZMu, Fctrzn_Lrn_W0, Fctrzn_Lrn_W))
  E_Z_SqE_W_Sq <- list("mRNA" = E_Z_SqE_W_Sq_update(view = c("mRNA"), ZMu_0, ZMu, Fctrzn_Lrn_W0, Fctrzn_Lrn_W))
  E_ZSqE_WSq <- list("mRNA" = E_ZSqE_WSq_update(view = c("mRNA"), ZMu_0, ZMuSq, Fctrzn_Lrn_W0, Fctrzn_Lrn_WSq))
  E_ZWSq <- list("mRNA" = E_ZWSq_update(view = c("mRNA"), E_ZE_W, ZMuSq, E_Z_SqE_W_Sq, E_ZSqE_WSq))
  Zeta <- Zeta_calculation("mRNA", likelihoods, E_ZWSq, E_ZE_W)
  Tau_m <- Tau_calculation("mRNA", likelihoods, Zeta, Tau)
  YGauss <- YGauss_calculation("mRNA", likelihoods, YTrg, Zeta, Tau_m, CenterTrg, PoisRateCstnt)
  ZMu <- ZMu_calculation(view = "mRNA", k = 1, Fctrzn_Lrn_W, Fctrzn_Lrn_W0, Tau, ZMu_0, ZMu, YGauss)
  # ZVar <- ZVar_calculation(view = "mRNA", Tau, Fctrzn_Lrn_WSq)
  # ZMu_0 <- rep(1,dim(ZVar)[1])
  # ZMu <- ZMu_calculation(view = "mRNA", k = 1, Fctrzn_Lrn_W, Fctrzn_Lrn_W0, Tau, ZMu_0, ZMu, YGauss = )
  print("to do")
})

test_that("ELBO_calculation", {
  print("to do")
})

test_that("E_ZE_W_update", {
  E_ZE_W <- E_ZE_W_update(view = "mRNA", ZMu_0, ZMu, Fctrzn_Lrn_W0, Fctrzn_Lrn_W)
  expect_equal(rownames(E_ZE_W), rownames(ZMu))
  expect_equal(colnames(E_ZE_W), names(Fctrzn_Lrn_W0$mRNA))
})

test_that("E_Z_SqE_W_Sq_update", {
  E_Z_SqE_W_Sq <- E_Z_SqE_W_Sq_update(view = "mRNA", ZMu_0, ZMu, Fctrzn_Lrn_W0, Fctrzn_Lrn_W)
  expect_equal(rownames(E_Z_SqE_W_Sq), rownames(ZMu))
  expect_equal(colnames(E_Z_SqE_W_Sq), names(Fctrzn_Lrn_W0$mRNA))
})

test_that("E_ZSqE_WSq_update", {
  E_ZSqE_WSq <- E_ZSqE_WSq_update(view = "mRNA", ZMu_0, ZMuSq, Fctrzn_Lrn_W0, Fctrzn_Lrn_WSq)
  expect_equal(rownames(E_ZSqE_WSq), rownames(ZMu_0))
  expect_equal(colnames(E_ZSqE_WSq), names(Fctrzn_Lrn_W0$mRNA))
})

test_that("E_ZWSq_update", {
  E_ZE_W <- list("mRNA" = E_ZE_W_update(view = "mRNA", ZMu_0, ZMu, Fctrzn_Lrn_W0, Fctrzn_Lrn_W))
  E_Z_SqE_W_Sq <- list("mRNA" = E_Z_SqE_W_Sq_update(view = "mRNA", ZMu_0, ZMu, Fctrzn_Lrn_W0, Fctrzn_Lrn_W))
  E_ZSqE_WSq <- list("mRNA" = E_ZSqE_WSq_update(view = "mRNA", ZMu_0, ZMuSq, Fctrzn_Lrn_W0, Fctrzn_Lrn_WSq))
  E_ZWSq <- E_ZWSq_update(view = "mRNA", E_ZE_W, ZMuSq, E_Z_SqE_W_Sq, E_ZSqE_WSq)
  expect_equal(dim(E_ZWSq), dim(E_Z_SqE_W_Sq$mRNA))
  expect_equal(dim(E_ZWSq), dim(E_ZSqE_WSq$mRNA))
})

test_that("VarExplFun", {
  print("to do")
})

test_that("transferLearning_function", {
  print("to do")
})
