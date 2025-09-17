skip_issue_2_solve <- function(message) {
    skip(paste0("ISSUE 2 SOLVE - ASK_DAVID - ", message))
}

test_that("TCGATargetDataPrefiltering", {
  YTrg_mRNA <- TCGATargetDataPrefiltering(view = "mRNA", brcds_SS = Trg$brcds_SS, SS = 1, YTrg_list, Lrn_Fctrzn_init)
  expect_equal(class(YTrg_mRNA), class(SummarizedExperiment()))
  expect_equal(dim(YTrg_mRNA), dim(YTrg_list$mRNA))
  expect_equal(colnames(YTrg_mRNA), colnames(YTrg_list$mRNA))
  expect_in(rownames(YTrg_mRNA), rownames(YTrg_list$mRNA))
  YTrg_miRNA <- TCGATargetDataPrefiltering(view = "miRNA", brcds_SS = Trg$brcds_SS, SS = 1, YTrg_list, Lrn_Fctrzn_init)
  expect_equal(class(YTrg_miRNA), class(SummarizedExperiment()))
  expect_equal(dim(YTrg_miRNA)[[1]], sum(rowVars(assay(YTrg_list$miRNA))>0))
  expect_equal(colnames(YTrg_miRNA), colnames(YTrg_list$miRNA))
  expect_in(rownames(YTrg_miRNA), rownames(YTrg_list$miRNA))
  YTrg_DNAme <- TCGATargetDataPrefiltering(view = "DNAme", brcds_SS = Trg$brcds_SS, SS = 1, YTrg_list, Lrn_Fctrzn_init)
  expect_equal(class(YTrg_DNAme), class(SummarizedExperiment()))
  expect_equal(dim(YTrg_DNAme), dim(YTrg_list$DNAme))
  expect_equal(colnames(YTrg_DNAme), colnames(YTrg_list$DNAme))
  expect_in(rownames(YTrg_DNAme), rownames(YTrg_list$DNAme))
  YTrg_SNV <- TCGATargetDataPrefiltering(view = "SNV", brcds_SS = Trg$brcds_SS, SS = 1, YTrg_list, Lrn_Fctrzn_init)
  expect_equal(class(YTrg_SNV), class(matrix()))
  expect_equal(dim(YTrg_SNV)[[1]], sum(rowVars(YTrg_list$SNV)>0))
  expect_equal(colnames(YTrg_SNV), colnames(YTrg_list$SNV))
  expect_in(rownames(YTrg_SNV), rownames(YTrg_list$SNV))
  expect_equal(colnames(YTrg_mRNA), colnames(YTrg_miRNA))
  expect_equal(colnames(YTrg_mRNA), colnames(YTrg_DNAme))
  expect_equal(colnames(YTrg_mRNA), colnames(YTrg_SNV))
})

test_that("TCGATargetDataPreparation", {
  Trg_exp <- TCGATargetDataPreparation(views, YTrg_list, Trg$brcds_SS, SS = 1, Lrn_Fctrzn, smpls, normalization = FALSE, transformation = FALSE, Lrn_meta)
  expect_equal(names(Trg_exp), views)
  expect_equal(Trg_exp, YTrg_prep)
})

test_that("mRNA_addVersion", {
    expdat_mRNA_wov_up = mRNA_addVersion(expdat = df, Lrndat)
    expect_in(rownames(expdat_mRNA_wov_up), rownames(Lrndat))
})

test_that("GeoMeanFun", {
  x <- assay(YTrg_list$mRNA)[41,]
  GeoMean <- GeoMeanFun(x)
  expect_equal(length(GeoMean), 1)
  expect_equal(GeoMean, exp(sum(log(x[x > 0]))/length(x)))
})

test_that("GeoMeans_Lrn_init", {
  # YTrg_list$mRNA <- assay(YTrg_list$mRNA)
  # YTrg_list$miRNA <- assay(YTrg_list$miRNA)
  # YTrg_list$DNAme <- assay(YTrg_list$DNAme)
  # YTrg_prep <- sapply(views, TargetDataPrefiltering, YTrg_list, Lrn_Fctrzn, smpls)
  GeoMeans <- GeoMeans_Lrn_init("mRNA", Lrn_meta, rownames(YTrg_prep$mRNA))
  expect_equal(length(GeoMeans), length(rownames(YTrg_prep$mRNA)))
  expect_in(GeoMeans, Lrn_meta$GeoMeans_mRNA)
  GeoMeans <- GeoMeans_Lrn_init("miRNA", Lrn_meta, rownames(YTrg_prep$miRNA))
  expect_equal(length(GeoMeans), length(rownames(YTrg_prep$miRNA)))
  expect_in(GeoMeans, Lrn_meta$GeoMeans_miRNA)
  GeoMeans <- GeoMeans_Lrn_init("DNAme", Lrn_meta, rownames(YTrg_prep$DNAme))
  expect_equal(GeoMeans, numeric())
  GeoMeans <- GeoMeans_Lrn_init("SNV", Lrn_meta, rownames(YTrg_prep$SNV))
  expect_equal(GeoMeans, numeric())
})


test_that("countsTransformation", {
  mat <- assay(YTrg_list$mRNA)
  expdat_trans <- countsTransformation(mat, 50)
  expect_equal(nrow(expdat_trans), 50)
  expect_equal(names(expdat_trans), names(mat))
  expect_in(rownames(expdat_trans), rownames(mat))
})

test_that("countsNormalization_Lrn", {
  # YTrg_list$mRNA <- assay(YTrg_list$mRNA)
  # expdat <- TargetDataPrefiltering(view = "mRNA", YTrg_list, Lrn_Fctrzn, smpls)
  expdat <- apply(YTrg_prep$mRNA, MARGIN = 2, round)
  expdat <- SummarizedExperiment(assays = expdat)
  expdat_norm <- countsNormalization(expdat, GeoMeans = "Lrn")
  expect_equal(length(expdat_norm), 2)
  expect_equal(dim(expdat_norm$counts), dim(expdat))
  expect_equal(rownames(expdat_norm$counts), rownames(expdat))
  expect_equal(colnames(expdat_norm$counts), colnames(expdat))
})

test_that("countsNormalization_Trg", {
  # YTrg_list$mRNA <- assay(YTrg_list$mRNA)
  # expdat = TargetDataPrefiltering(view = "mRNA", YTrg_list, Lrn_Fctrzn, smpls)
  expdat <- apply(YTrg_prep$mRNA, MARGIN = 2, round)
  expdat <- SummarizedExperiment(assays = expdat)
  expdat_norm <- countsNormalization(expdat, GeoMeans = "Trg")
  expect_equal(length(expdat_norm), 1)
  expect_equal(dim(expdat_norm$counts), dim(expdat))
  expect_equal(rownames(expdat_norm$counts), rownames(expdat))
  expect_equal(colnames(expdat_norm$counts), colnames(expdat))
})

test_that("countsNormalization_num", {
  # YTrg_list$mRNA <- assay(YTrg_list$mRNA)
  # expdat = TargetDataPrefiltering(view = "mRNA", YTrg_list, Lrn_Fctrzn, smpls)
  expdat <- apply(YTrg_prep$mRNA, MARGIN = 2, round)
  expdat <- SummarizedExperiment(assays = expdat)
  GeoMeans <- GeoMeans_Lrn_init("mRNA", Lrn_meta, rownames(expdat))
  expdat_norm <- countsNormalization(expdat, GeoMeans = GeoMeans)
  expect_equal(length(expdat_norm), 1)
  expect_equal(dim(expdat_norm$counts), dim(expdat))
  expect_equal(rownames(expdat_norm$counts), rownames(expdat))
  expect_equal(colnames(expdat_norm$counts), colnames(expdat))
})

test_that("TargetDataPrefiltering", {
  ## SE object into matrix
  YTrg_list$mRNA <- assay(YTrg_list$mRNA)
  YTrg_list$miRNA <- assay(YTrg_list$miRNA)
  YTrg_list$DNAme <- assay(YTrg_list$DNAme)
  ## Change names of some rows - expected removed
  rownames(YTrg_list$mRNA)[c(2,10)] <- "newNames"
  rownames(YTrg_list$DNAme)[c(2,46,50)] <- "newNames"
  ##
  YTrg = sapply(views, TargetDataPrefiltering, YTrg_list, Lrn_Fctrzn, smpls)
  expect_equal(names(YTrg), names(YTrg_list))
  expect_equal(nrow(YTrg$mRNA), nrow(YTrg_list$mRNA) - 2)
  expect_equal(nrow(YTrg$miRNA), sum(rowVars(YTrg_list$miRNA) > 0))
  expect_equal(nrow(YTrg$DNAme), nrow(YTrg_list$DNAme) - 3)
  expect_equal(nrow(YTrg$SNV), sum(rowVars(YTrg_list$SNV) > 0))
  expect_equal(colnames(YTrg$mRNA), smpls)
  expect_equal(colnames(YTrg$miRNA), smpls)
  expect_equal(colnames(YTrg$DNAme), smpls)
  expect_equal(colnames(YTrg$SNV), smpls)
  expect_in(rownames(YTrg$mRNA), Lrn_Fctrzn@features_metadata[Lrn_Fctrzn@features_metadata$view == "mRNA",1])
  expect_in(rownames(YTrg$miRNA), Lrn_Fctrzn@features_metadata[Lrn_Fctrzn@features_metadata$view == "miRNA",1])
  expect_in(rownames(YTrg$DNAme), Lrn_Fctrzn@features_metadata[Lrn_Fctrzn@features_metadata$view == "DNAme",1])
  expect_in(rownames(YTrg$SNV), Lrn_Fctrzn@features_metadata[Lrn_Fctrzn@features_metadata$view == "SNV",1])
})

test_that("preprocessCountsData_SIMPLE", {
  ## NO TRANSFORMATION NO NORMALIZATION
  YTrg = sapply(views, preprocessCountsData, YTrg_prep, normalization = FALSE, Lrn_meta, transformation = FALSE)
  expect_equal(names(YTrg), names(YTrg_prep))
  expect_equal(YTrg$mRNA, YTrg_prep$mRNA)
  expect_equal(YTrg$miRNA, YTrg_prep$miRNA)
  expect_equal(YTrg$DNAme, YTrg_prep$DNAme)
  expect_equal(YTrg$SNV, YTrg_prep$SNV)
})

test_that("preprocessCountsData_TRANS", {
  # YTrg_pre = sapply(views, TargetDataPrefiltering, YTrg_list, Fctrzn, smpls)
  ## TRANSFORMATION NO NORMALIZATION
  YTrg = sapply(views, preprocessCountsData, YTrg_prep, normalization = FALSE, Lrn_meta, transformation = TRUE)
  expect_equal(names(YTrg), names(YTrg_prep))
  expect_equal(names(YTrg$mRNA), names(YTrg_prep$mRNA))
  expect_equal(rownames(YTrg$mRNA), rownames(YTrg_prep$mRNA))
  expect_equal(names(YTrg$miRNA), names(YTrg_prep$miRNA))
  expect_equal(rownames(YTrg$miRNA), rownames(YTrg_prep$miRNA))
  expect_equal(names(YTrg$DNAme), names(YTrg_prep$DNAme))
  expect_equal(rownames(YTrg$DNAme), rownames(YTrg_prep$DNAme))
  expect_equal(YTrg$DNAme, YTrg_prep$DNAme)
  expect_equal(names(YTrg$SNV), names(YTrg_prep$SNV))
  expect_equal(rownames(YTrg$SNV), rownames(YTrg_prep$SNV))
  expect_equal(YTrg$SNV, YTrg_prep$SNV)
})

test_that("preprocessCountsData_NORM_LRN", {
  ## NO TRANSFORMATION NORMALIZATION WITH LRN DATASET
  YTrg <- YTrg_prep
  YTrg$mRNA <- apply(YTrg$mRNA, MARGIN = 2, round)
  YTrg$DNA <- apply(YTrg$DNA, MARGIN = 2, round)
  expdat_norm <- preprocessCountsData(view = "mRNA", YTrg, normalization = "Lrn", Lrn_meta, transformation = FALSE)
  expect_equal(dim(expdat_norm), dim(YTrg$mRNA))
  expect_equal(rownames(expdat_norm), rownames(YTrg$mRNA))
  expect_equal(colnames(expdat_norm), colnames(YTrg$mRNA))
  expdat_norm <- preprocessCountsData(view = "DNAme", YTrg, normalization = "Lrn", Lrn_meta, transformation = FALSE)
  expect_equal(expdat_norm, YTrg$DNAme)
})

test_that("preprocessCountsData_NORM_TRG", {
  ## NO TRANSFORMATION NORMALIZATION WITH TRG DATASET
  YTrg <- YTrg_prep
  YTrg$mRNA <- apply(YTrg$mRNA, MARGIN = 2, round)
  YTrg$DNA <- apply(YTrg$DNA, MARGIN = 2, round)
  expdat_norm <- preprocessCountsData(view = "mRNA", YTrg, normalization = "Trg", Lrn_meta, transformation = FALSE)
  expect_equal(dim(expdat_norm), dim(YTrg$mRNA))
  expect_equal(rownames(expdat_norm), rownames(YTrg$mRNA))
  expect_equal(colnames(expdat_norm), colnames(YTrg$mRNA))
  expdat_norm <- preprocessCountsData(view = "DNAme", YTrg, normalization = "Trg", Lrn_meta, transformation = FALSE)
  expect_equal(expdat_norm, YTrg$DNAme)
})

test_that("TargetDataPreparation", {
  expdat_prep <- TargetDataPreparation(views, YTrg_prep, Lrn_Fctrzn, smpls, Lrn_meta, normalization = FALSE, transformation = FALSE)
  expect_equal(names(expdat_prep), views)
  expect_equal(expdat_prep, YTrg_prep)
})

test_that("initTransferLearningParamaters", {
  YTrgFtrs <- lapply(YTrg_prep, rownames)
  TL_param = initTransferLearningParamaters(YTrg_prep, views, Lrn_meta, Lrn_Fctrzn_init, likelihoods)
  expect_equal(names(TL_param), c("YTrg", "Fctrzn_Lrn_W0", "Fctrzn_Lrn_W", "Fctrzn_Lrn_WSq", "Tau", "TauLn"))
  expect_equal(names(TL_param$Fctrzn_Lrn_W0$mRNA), YTrgFtrs$mRNA)
  expect_equal(names(TL_param$Fctrzn_Lrn_W0$miRNA), YTrgFtrs$miRNA)
  expect_equal(names(TL_param$Fctrzn_Lrn_W0$DNAme), YTrgFtrs$DNAme)
  expect_equal(names(TL_param$Fctrzn_Lrn_W0$SNV), YTrgFtrs$SNV)
  expect_equal(rownames(TL_param$Fctrzn_Lrn_W$mRNA), YTrgFtrs$mRNA)
  expect_equal(rownames(TL_param$Fctrzn_Lrn_W$miRNA), YTrgFtrs$miRNA)
  expect_equal(rownames(TL_param$Fctrzn_Lrn_W$DNAme), YTrgFtrs$DNAme)
  expect_equal(rownames(TL_param$Fctrzn_Lrn_W$SNV), YTrgFtrs$SNV)
  expect_equal(rownames(TL_param$Fctrzn_Lrn_WSq$mRNA), YTrgFtrs$mRNA)
  expect_equal(rownames(TL_param$Fctrzn_Lrn_WSq$miRNA), YTrgFtrs$miRNA)
  expect_equal(rownames(TL_param$Fctrzn_Lrn_WSq$DNAme), YTrgFtrs$DNAme)
  expect_equal(rownames(TL_param$Fctrzn_Lrn_WSq$SNV), YTrgFtrs$SNV)
  expect_equal(nrow(TL_param$Tau$mRNA), ncol(YTrg_prep$mRNA))
  expect_equal(ncol(TL_param$Tau$mRNA), nrow(YTrg_prep$mRNA))
  expect_equal(nrow(TL_param$Tau$miRNA), ncol(YTrg_prep$miRNA))
  expect_equal(ncol(TL_param$Tau$miRNA), nrow(YTrg_prep$miRNA))
  expect_equal(nrow(TL_param$Tau$DNAme), ncol(YTrg_prep$DNAme))
  expect_equal(ncol(TL_param$Tau$DNAme), nrow(YTrg_prep$DNAme))
  expect_equal(nrow(TL_param$Tau$SNV), ncol(YTrg_prep$SNV))
  expect_equal(ncol(TL_param$Tau$SNV), nrow(YTrg_prep$SNV))
  expect_equal(nrow(TL_param$TauLn$mRNA), ncol(YTrg_prep$mRNA))
  expect_equal(ncol(TL_param$TauLn$mRNA), nrow(YTrg_prep$mRNA))
  expect_equal(nrow(TL_param$TauLn$miRNA), ncol(YTrg_prep$miRNA))
  expect_equal(ncol(TL_param$TauLn$miRNA), nrow(YTrg_prep$miRNA))
  expect_equal(nrow(TL_param$TauLn$DNAme), ncol(YTrg_prep$DNAme))
  expect_equal(ncol(TL_param$TauLn$DNAme), nrow(YTrg_prep$DNAme))
  expect_equal(TL_param$TauLn$SNV, t(numeric()))
  expect_equal(TL_param$YTrg$mRNA, t(YTrg_prep$mRNA))
  expect_equal(TL_param$YTrg$miRNA, t(YTrg_prep$miRNA))
  expect_equal(TL_param$YTrg$DNAme, t(YTrg_prep$DNAme))
  expect_equal(TL_param$YTrg$SNV, t(YTrg_prep$SNV))
})

test_that("Tau_init", {
  Tau <- Tau_init(views, Lrn_Fctrzn, Lrn_ModelFile)
  expect_equal(nrow(Tau$mRNA$group0), Lrn_Fctrzn@dimensions$D["mRNA"][[1]])
  expect_equal(nrow(Tau$miRNA$group0), Lrn_Fctrzn@dimensions$D["miRNA"][[1]])
  expect_equal(nrow(Tau$DNAme$group0), Lrn_Fctrzn@dimensions$D["DNAme"][[1]])
  expect_equal(nrow(Tau$SNV$group0), Lrn_Fctrzn@dimensions$D["SNV"][[1]])
})

##
test_that("TauLn_calculation_TRUE", {
  ## LrnSimple = TRUE / DON'T USE LrnFctrnDir TAU FILES
  TauLn <-  sapply(views, TauLn_calculation, likelihoods, Lrn_Fctrzn_init, Lrn_ModelDir, LrnSimple = TRUE)
  expect_equal(length(TauLn$mRNA), Lrn_Fctrzn_init@dimensions$D["mRNA"][[1]])
  expect_equal(TauLn$mRNA, log(Lrn_Fctrzn_init@expectations$Tau$mRNA$group0[,1]))
  expect_equal(length(TauLn$miRNA), Lrn_Fctrzn_init@dimensions$D["miRNA"][[1]])
  expect_equal(TauLn$miRNA, log(Lrn_Fctrzn_init@expectations$Tau$miRNA$group0[,1]))
  expect_equal(length(TauLn$DNAme), Lrn_Fctrzn_init@dimensions$D["DNAme"][[1]])
  expect_equal(TauLn$DNAme, log(Lrn_Fctrzn_init@expectations$Tau$DNAme$group0[,1]))
  expect_equal(TauLn$SNV, numeric())
})

test_that("TauLn_calculation_FALSE", {
  # skip_issue_2_solve("miRNA length is different between model and imported TauLn file")
  ## LrnSimple = FALSE / USE LrnFctrnDir TAU FILES
  mRNA <- TauLn_calculation(view = "mRNA", likelihoods, Lrn_Fctrzn_init, Lrn_ModelDir, LrnSimple = FALSE)
  expect_equal(length(mRNA), Lrn_Fctrzn_init@dimensions$D["mRNA"][[1]])
  miRNA <- TauLn_calculation(view = "miRNA", likelihoods, Lrn_Fctrzn_init, Lrn_ModelDir, LrnSimple = FALSE)
  expect_equal(length(miRNA), Lrn_Fctrzn_init@dimensions$D["miRNA"][[1]])
  DNAme <- TauLn_calculation(view = "DNAme", likelihoods, Lrn_Fctrzn_init, Lrn_ModelDir, LrnSimple = FALSE)
  expect_equal(length(DNAme), Lrn_Fctrzn_init@dimensions$D["DNAme"][[1]])
  SNV <- TauLn_calculation(view = "SNV", likelihoods, Lrn_Fctrzn_init, Lrn_ModelDir, LrnSimple = FALSE)
  expect_equal(SNV, numeric())
})

test_that("WSq_calculation_TRUE", {
  ## LrnSimple = TRUE / DON'T USE LrnFctrnDir W FILES
  WSq <- sapply(views, WSq_calculation, Lrn_Fctrzn_init, Lrn_ModelDir, LrnSimple = TRUE)
  expect_equal(nrow(WSq$mRNA), Lrn_Fctrzn_init@dimensions$D["mRNA"][[1]])
  expect_equal(WSq$mRNA, Lrn_Fctrzn_init@expectations$W[["mRNA"]]^2)
  expect_equal(nrow(WSq$miRNA), Lrn_Fctrzn_init@dimensions$D["miRNA"][[1]])
  expect_equal(WSq$miRNA, Lrn_Fctrzn_init@expectations$W[["miRNA"]]^2)
  expect_equal(nrow(WSq$DNAme), Lrn_Fctrzn_init@dimensions$D["DNAme"][[1]])
  expect_equal(WSq$DNAme, Lrn_Fctrzn_init@expectations$W[["DNAme"]]^2)
  expect_equal(nrow(WSq$SNV), Lrn_Fctrzn_init@dimensions$D["SNV"][[1]])
  expect_equal(WSq$SNV, Lrn_Fctrzn_init@expectations$W[["SNV"]]^2)
})

test_that("WSq_calculation_FALSE", {
  # skip_issue_2_solve("W files contain 84 samples and model contains 91 samples ... wrong files")
  ## LrnSimple = FALSE / USE LrnFctrnDir W FILES
  # WSq = sapply(views, WSq_calculation, Lrn_Fctrzn_init, Lrn_ModelDir, LrnSimple = FALSE)
  mRNA <- WSq_calculation(view = "mRNA", Lrn_Fctrzn_init, Lrn_ModelDir, LrnSimple = FALSE)
  expect_equal(nrow(mRNA), Lrn_Fctrzn_init@dimensions$D["mRNA"][[1]])
  miRNA <- WSq_calculation(view = "miRNA", Lrn_Fctrzn_init, Lrn_ModelDir, LrnSimple = FALSE)
  expect_equal(nrow(miRNA), Lrn_Fctrzn_init@dimensions$D["miRNA"][[1]])
  DNAme <- WSq_calculation(view = "DNAme", Lrn_Fctrzn_init, Lrn_ModelDir, LrnSimple = FALSE)
  expect_equal(nrow(DNAme), Lrn_Fctrzn_init@dimensions$D["DNAme"][[1]])
  SNV <- WSq_calculation(view = "SNV", Lrn_Fctrzn_init, Lrn_ModelDir, LrnSimple = FALSE)
  expect_equal(nrow(SNV), Lrn_Fctrzn_init@dimensions$D["SNV"][[1]])
})

test_that("W0_calculation_FALSE", {
  ## CenterTrg = FALSE
  W0 <- sapply(views, W0_calculation, CenterTrg = FALSE, Lrn_Fctrzn_init, Lrn_ModelDir)
  expect_equal(length(W0$mRNA), Lrn_Fctrzn_init@dimensions$D["mRNA"][[1]])
  expect_equal(length(W0$miRNA), Lrn_Fctrzn_init@dimensions$D["miRNA"][[1]])
  expect_equal(length(W0$DNAme), Lrn_Fctrzn_init@dimensions$D["DNAme"][[1]])
  expect_equal(length(W0$SNV), Lrn_Fctrzn_init@dimensions$D["SNV"][[1]])
})

test_that("W0_calculation_TRUE", {
  ## CenterTrg = TRUE
  W0 = sapply(views, W0_calculation, CenterTrg = TRUE, Lrn_Fctrzn_init, Lrn_ModelDir)
  expect_equal(length(W0$mRNA), Lrn_Fctrzn_init@dimensions$D["mRNA"][[1]])
  expect_equal(sum(W0$mRNA), 0)
  expect_equal(length(W0$miRNA), Lrn_Fctrzn_init@dimensions$D["miRNA"][[1]])
  expect_equal(sum(W0$miRNA), 0)
  expect_equal(length(W0$DNAme), Lrn_Fctrzn_init@dimensions$D["DNAme"][[1]])
  expect_equal(sum(W0$DNAme), 0)
  expect_equal(length(W0$SNV), Lrn_Fctrzn_init@dimensions$D["SNV"][[1]])
  expect_equal(sum(W0$SNV), 0)
})

test_that("intercepts_calculation", {
  print("to do")
})

test_that("Zeta_calculation", {
  likelihoods <- list("mRNA" = "gaussian", "miRNA" = "bernoulli")
  E_ZE_W <- list("mRNA" = E_ZE_W_update(view = "mRNA", TL_param$ZMu_0, TL_param$ZMu, TL_param$Fctrzn_Lrn_W0, TL_param$Fctrzn_Lrn_W))
  E_Z_SqE_W_Sq <- list("mRNA" = E_Z_SqE_W_Sq_update(view = "mRNA", TL_param$ZMu_0, TL_param$ZMu, TL_param$Fctrzn_Lrn_W0, TL_param$Fctrzn_Lrn_W))
  E_ZSqE_WSq <- list("mRNA" = E_ZSqE_WSq_update(view = "mRNA", TL_param$ZMu_0, TL_param$ZMuSq, TL_param$Fctrzn_Lrn_W0, TL_param$Fctrzn_Lrn_WSq))
  E_ZWSq <- list("mRNA" = E_ZWSq_update(view = "mRNA", E_ZE_W, TL_param$ZMuSq, E_Z_SqE_W_Sq, E_ZSqE_WSq))
  E_ZWSq$miRNA <- E_ZWSq$mRNA
  E_ZE_W$miRNA <- E_ZE_W$mRNA
  zeta <- Zeta_calculation(view = "mRNA", likelihoods, E_ZWSq, E_ZE_W)
  expect_equal(zeta, E_ZE_W$mRNA)
  zeta <- Zeta_calculation(view = "miRNA", likelihoods, E_ZWSq, E_ZE_W)
  expect_equal(zeta, sqrt(E_ZWSq$miRNA))
})

test_that("Tau_calculation", {
  likelihoods <- list("mRNA" = "gaussian", "miRNA" = "bernoulli")
  E_ZE_W <- list("mRNA" = E_ZE_W_update(view = "mRNA", TL_param$ZMu_0, TL_param$ZMu, TL_param$Fctrzn_Lrn_W0, TL_param$Fctrzn_Lrn_W))
  E_Z_SqE_W_Sq <- list("mRNA" = E_Z_SqE_W_Sq_update(view = "mRNA", TL_param$ZMu_0, TL_param$ZMu, TL_param$Fctrzn_Lrn_W0, TL_param$Fctrzn_Lrn_W))
  E_ZSqE_WSq <- list("mRNA" = E_ZSqE_WSq_update(view = "mRNA", TL_param$ZMu_0, TL_param$ZMuSq, TL_param$Fctrzn_Lrn_W0, TL_param$Fctrzn_Lrn_WSq))
  E_ZWSq <- list("mRNA" = E_ZWSq_update(view = "mRNA", E_ZE_W, TL_param$ZMuSq, E_Z_SqE_W_Sq, E_ZSqE_WSq))
  E_ZWSq$miRNA <- E_ZWSq$mRNA
  E_ZE_W$miRNA <- E_ZE_W$mRNA
  Zeta <- list("mRNA" = Zeta_calculation(view = "mRNA", likelihoods, E_ZWSq, E_ZE_W))
  Zeta$miRNA <- Zeta$mRNA
  Tau_m <- Tau_calculation(view = "mRNA", likelihoods, Zeta, TL_param$Tau)
  expect_equal(Tau_m, TL_param$Tau$mRNA)
  Tau_m <- Tau_calculation(view = "miRNA", likelihoods, Zeta, TL_param$Tau)
  expect_equal(Tau_m, (1/2)*(1/Zeta$miRNA)*tanh(Zeta$miRNA/2))
})

test_that("YGauss_calculation", {
  likelihoods <- list("mRNA" = "gaussian", "miRNA" = "bernoulli")
  E_ZE_W <- sapply(c("mRNA", "miRNA"), E_ZE_W_update, TL_param$ZMu_0, TL_param$ZMu, TL_param$Fctrzn_Lrn_W0, TL_param$Fctrzn_Lrn_W)
  E_Z_SqE_W_Sq <- sapply(c("mRNA", "miRNA"), E_Z_SqE_W_Sq_update, TL_param$ZMu_0, TL_param$ZMu, TL_param$Fctrzn_Lrn_W0, TL_param$Fctrzn_Lrn_W)
  E_ZSqE_WSq <- sapply(c("mRNA", "miRNA"), E_ZSqE_WSq_update, TL_param$ZMu_0, TL_param$ZMuSq, TL_param$Fctrzn_Lrn_W0, TL_param$Fctrzn_Lrn_WSq)
  E_ZWSq <- sapply(c("mRNA", "miRNA"), E_ZWSq_update, E_ZE_W, TL_param$ZMuSq, E_Z_SqE_W_Sq, E_ZSqE_WSq)
  Zeta <- sapply(c("mRNA", "miRNA"), Zeta_calculation, likelihoods, E_ZWSq, E_ZE_W)
  Tau_m <- sapply(c("mRNA", "miRNA"), Tau_calculation, likelihoods, Zeta, TL_param$Tau)
  YGauss <- YGauss_calculation("mRNA", likelihoods, TL_param$YTrg, Zeta, Tau_m, CenterTrg, PoisRateCstnt)
  expect_equal(YGauss, TL_param$YTrg$mRNA)
  YGauss <- YGauss_calculation("miRNA", likelihoods, TL_param$YTrg, Zeta, Tau_m, CenterTrg, PoisRateCstnt)
  expect_equal(YGauss, (2*TL_param$YTrg$miRNA - 1) / (2*Tau_m$miRNA))
})

test_that("ZVar_calculation", {
  # ASK_DAVID - matrix looks strange - number inside are same
  ZVar <- ZVar_calculation(view = "mRNA", Tau = TL_param$Tau, Fctrzn_Lrn_WSq = TL_param$Fctrzn_Lrn_WSq)
  expect_equal(nrow(ZVar), nrow(TL_param$Tau$mRNA))
  expect_equal(ncol(ZVar), ncol(TL_param$Fctrzn_Lrn_WSq$mRNA))
  expect_equal(ZVar, TL_param$Tau$mRNA %*% TL_param$Fctrzn_Lrn_WSq$mRNA)
})

test_that("ZMu_calculation", {
  E_ZE_W <- list("mRNA" = E_ZE_W_update(view = c("mRNA"), TL_param$ZMu_0, TL_param$ZMu, TL_param$Fctrzn_Lrn_W0, TL_param$Fctrzn_Lrn_W))
  E_Z_SqE_W_Sq <- list("mRNA" = E_Z_SqE_W_Sq_update(view = c("mRNA"), TL_param$ZMu_0, TL_param$ZMu, TL_param$Fctrzn_Lrn_W0, TL_param$Fctrzn_Lrn_W))
  E_ZSqE_WSq <- list("mRNA" = E_ZSqE_WSq_update(view = c("mRNA"), TL_param$ZMu_0, TL_param$ZMuSq, TL_param$Fctrzn_Lrn_W0, TL_param$Fctrzn_Lrn_WSq))
  E_ZWSq <- list("mRNA" = E_ZWSq_update(view = c("mRNA"), E_ZE_W, TL_param$ZMuSq, E_Z_SqE_W_Sq, E_ZSqE_WSq))
  Zeta <- Zeta_calculation("mRNA", likelihoods, E_ZWSq, E_ZE_W)
  Tau_m <- Tau_calculation("mRNA", likelihoods, Zeta, TL_param$Tau)
  YGauss <- list("mRNA" = YGauss_calculation("mRNA", likelihoods, TL_param$YTrg, Zeta, Tau_m, CenterTrg, PoisRateCstnt))
  ZMu_m <- ZMu_calculation(view = "mRNA", k = 1, TL_param$Fctrzn_Lrn_W, TL_param$Fctrzn_Lrn_W0, TL_param$Tau, TL_param$ZMu_0, TL_param$ZMu, YGauss)
  expect_equal(length(ZMu_m), nrow(TL_param$Tau$mRNA))
  expect_equal(length(ZMu_m), length(TL_param$ZMu_0))
})

test_that("ELBO_calculation", {
  likelihoods <- list("mRNA" = "gaussian", "miRNA" = "gaussian", "DNAme" = "poisson", "SNV" = "bernoulli")
  E_ZE_W <- sapply(views, E_ZE_W_update, TL_param$ZMu_0, TL_param$ZMu, TL_param$Fctrzn_Lrn_W0, TL_param$Fctrzn_Lrn_W)
  E_Z_SqE_W_Sq <- sapply(views, E_Z_SqE_W_Sq_update, TL_param$ZMu_0, TL_param$ZMu, TL_param$Fctrzn_Lrn_W0, TL_param$Fctrzn_Lrn_W)
  E_ZSqE_WSq <- sapply(views, E_ZSqE_WSq_update, TL_param$ZMu_0, TL_param$ZMuSq, TL_param$Fctrzn_Lrn_W0, TL_param$Fctrzn_Lrn_WSq)
  E_ZWSq <- sapply(views, E_ZWSq_update, E_ZE_W, TL_param$ZMuSq, E_Z_SqE_W_Sq, E_ZSqE_WSq)
  Zeta <- sapply(views, Zeta_calculation, likelihoods, E_ZWSq, E_ZE_W)
  Tau_m <- sapply(views, Tau_calculation, likelihoods, Zeta, TL_param$Tau)
  YGauss <- sapply(views, YGauss_calculation, likelihoods, TL_param$YTrg, Zeta, Tau_m, CenterTrg, PoisRateCstnt)
  ELBO_L <- ELBO_calculation(view = "mRNA", likelihoods, TL_param$Tau, TL_param$TauLn, E_ZWSq, E_ZE_W, Zeta, TL_param$YTrg, YGauss, PoisRateCstnt)
  expect_contains(ELBO_L, numeric())
  ELBO_L <- ELBO_calculation(view = "miRNA", likelihoods, TL_param$Tau, TL_param$TauLn, E_ZWSq, E_ZE_W, Zeta, TL_param$YTrg, YGauss, PoisRateCstnt)
  expect_contains(ELBO_L, numeric())
  ELBO_L <- ELBO_calculation(view = "DNAme", likelihoods, TL_param$Tau, TL_param$TauLn, E_ZWSq, E_ZE_W, Zeta, TL_param$YTrg, YGauss, PoisRateCstnt)
  expect_contains(ELBO_L, numeric())
})

test_that("E_ZE_W_update", {
  E_ZE_W <- E_ZE_W_update(view = "mRNA", TL_param$ZMu_0, TL_param$ZMu, TL_param$Fctrzn_Lrn_W0, TL_param$Fctrzn_Lrn_W)
  expect_equal(rownames(E_ZE_W), rownames(TL_param$ZMu))
  expect_equal(colnames(E_ZE_W), names(TL_param$Fctrzn_Lrn_W0$mRNA))
})

test_that("E_Z_SqE_W_Sq_update", {
  E_Z_SqE_W_Sq <- E_Z_SqE_W_Sq_update(view = "mRNA", TL_param$ZMu_0, TL_param$ZMu, TL_param$Fctrzn_Lrn_W0, TL_param$Fctrzn_Lrn_W)
  expect_equal(rownames(E_Z_SqE_W_Sq), rownames(TL_param$ZMu))
  expect_equal(colnames(E_Z_SqE_W_Sq), names(TL_param$Fctrzn_Lrn_W0$mRNA))
})

test_that("E_ZSqE_WSq_update", {
  E_ZSqE_WSq <- E_ZSqE_WSq_update(view = "mRNA", TL_param$ZMu_0, TL_param$ZMuSq, TL_param$Fctrzn_Lrn_W0, TL_param$Fctrzn_Lrn_WSq)
  expect_equal(rownames(E_ZSqE_WSq), rownames(TL_param$ZMu_0))
  expect_equal(colnames(E_ZSqE_WSq), names(TL_param$Fctrzn_Lrn_W0$mRNA))
})

test_that("E_ZWSq_update", {
  E_ZE_W <- list("mRNA" = E_ZE_W_update(view = "mRNA", TL_param$ZMu_0, TL_param$ZMu, TL_param$Fctrzn_Lrn_W0, TL_param$Fctrzn_Lrn_W))
  E_Z_SqE_W_Sq <- list("mRNA" = E_Z_SqE_W_Sq_update(view = "mRNA", TL_param$ZMu_0, TL_param$ZMu, TL_param$Fctrzn_Lrn_W0, TL_param$Fctrzn_Lrn_W))
  E_ZSqE_WSq <- list("mRNA" = E_ZSqE_WSq_update(view = "mRNA", TL_param$ZMu_0, TL_param$ZMuSq, TL_param$Fctrzn_Lrn_W0, TL_param$Fctrzn_Lrn_WSq))
  E_ZWSq <- E_ZWSq_update(view = "mRNA", E_ZE_W, TL_param$ZMuSq, E_Z_SqE_W_Sq, E_ZSqE_WSq)
  expect_equal(dim(E_ZWSq), dim(E_Z_SqE_W_Sq$mRNA))
  expect_equal(dim(E_ZWSq), dim(E_ZSqE_WSq$mRNA))
})

test_that("VarExplFun", {
  print("to do")
})

test_that("transferLearning_function", {
  print("to do")
})

