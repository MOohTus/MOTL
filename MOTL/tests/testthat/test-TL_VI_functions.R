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
  print("to do")
})

test_that("GeoMeans_Lrn_init", {
  print("to do")
})

test_that("preprocessCountsData", {
  print("to do")
})

test_that("TargetDataPreparation", {
  print("to do")
})

test_that("initTransferLearningParamaters", {
  print("to do")
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

## ASK_DAVID FAUDRA ESSAYER AVEC LE FICHIER
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

## ASK_DAVID / REVOIR L'APPEL AU FICHIER
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
  print("to do")
})

test_that("Tau_calculation", {
  print("to do")
})

test_that("YGauss_calculation", {
  print("to do")
})

test_that("ZVar_calculation", {
  print("to do")
})

test_that("ZMu_calculation", {
  print("to do")
})

test_that("ELBO_calculation", {
  print("to do")
})

test_that("E_ZE_W_update", {
  print("to do")
})

test_that("E_Z_SqE_W_Sq_update", {
  print("to do")
})

test_that("E_ZWSq_update", {
  print("to do")
})

test_that("VarExplFun", {
  print("to do")
})

test_that("transferLearning_function", {
  print("to do")
})











