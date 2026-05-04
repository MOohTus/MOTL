test_that("Zeta_calculation", {
  likelihoods <- list("mRNA" = "gaussian", "miRNA" = "bernoulli")
  E_ZE_W <- list(
    "mRNA" = E_ZE_W_update(
      view = "mRNA",
      TL_param$ZMu_0,
      TL_param$ZMu,
      TL_param$Fctrzn_Lrn_W0,
      TL_param$Fctrzn_Lrn_W
    )
  )
  E_Z_SqE_W_Sq <- list(
    "mRNA" = E_Z_SqE_W_Sq_update(
      view = "mRNA",
      TL_param$ZMu_0,
      TL_param$ZMu,
      TL_param$Fctrzn_Lrn_W0,
      TL_param$Fctrzn_Lrn_W
    )
  )
  E_ZSqE_WSq <- list(
    "mRNA" = E_ZSqE_WSq_update(
      view = "mRNA",
      TL_param$ZMu_0,
      TL_param$ZMuSq,
      TL_param$Fctrzn_Lrn_W0,
      TL_param$Fctrzn_Lrn_WSq
    )
  )
  E_ZWSq <- list("mRNA" = E_ZWSq_update(view = "mRNA",
                                        E_ZE_W,
                                        TL_param$ZMuSq,
                                        E_Z_SqE_W_Sq,
                                        E_ZSqE_WSq))
  E_ZWSq$miRNA <- E_ZWSq$mRNA
  E_ZE_W$miRNA <- E_ZE_W$mRNA
  zeta <- Zeta_calculation(view = "mRNA", likelihoods, E_ZWSq, E_ZE_W)
  expect_equal(zeta, E_ZE_W$mRNA)
  zeta <- Zeta_calculation(view = "miRNA", likelihoods, E_ZWSq, E_ZE_W)
  expect_equal(zeta, sqrt(E_ZWSq$miRNA))
})

test_that("Tau_calculation", {
  likelihoods <- list("mRNA" = "gaussian", "miRNA" = "bernoulli")
  E_ZE_W <- list(
    "mRNA" = E_ZE_W_update(
      view = "mRNA",
      TL_param$ZMu_0,
      TL_param$ZMu,
      TL_param$Fctrzn_Lrn_W0,
      TL_param$Fctrzn_Lrn_W
    )
  )
  E_Z_SqE_W_Sq <- list(
    "mRNA" = E_Z_SqE_W_Sq_update(
      view = "mRNA",
      TL_param$ZMu_0,
      TL_param$ZMu,
      TL_param$Fctrzn_Lrn_W0,
      TL_param$Fctrzn_Lrn_W
    )
  )
  E_ZSqE_WSq <- list(
    "mRNA" = E_ZSqE_WSq_update(
      view = "mRNA",
      TL_param$ZMu_0,
      TL_param$ZMuSq,
      TL_param$Fctrzn_Lrn_W0,
      TL_param$Fctrzn_Lrn_WSq
    )
  )
  E_ZWSq <- list("mRNA" = E_ZWSq_update(view = "mRNA",
                                        E_ZE_W,
                                        TL_param$ZMuSq,
                                        E_Z_SqE_W_Sq,
                                        E_ZSqE_WSq))
  E_ZWSq$miRNA <- E_ZWSq$mRNA
  E_ZE_W$miRNA <- E_ZE_W$mRNA
  Zeta <- list("mRNA" = Zeta_calculation(view = "mRNA",
                                         likelihoods,
                                         E_ZWSq,
                                         E_ZE_W))
  Zeta$miRNA <- Zeta$mRNA
  Tau_m <- Tau_calculation(view = "mRNA", likelihoods, Zeta, TL_param$Tau)
  expect_equal(Tau_m, TL_param$Tau$mRNA)
  Tau_m <- Tau_calculation(view = "miRNA", likelihoods, Zeta, TL_param$Tau)
  expect_equal(Tau_m, (1 / 2) * (1 / Zeta$miRNA) * tanh(Zeta$miRNA / 2))
})

test_that("YGauss_calculation", {
  views = c("mRNA" = "mRNA", "miRNA" = "miRNA")
  likelihoods <- list("mRNA" = "gaussian", "miRNA" = "bernoulli")
  E_ZE_W <- lapply(
    views,
    E_ZE_W_update,
    TL_param$ZMu_0,
    TL_param$ZMu,
    TL_param$Fctrzn_Lrn_W0,
    TL_param$Fctrzn_Lrn_W
  )
  E_Z_SqE_W_Sq <- lapply(
    views,
    E_Z_SqE_W_Sq_update,
    TL_param$ZMu_0,
    TL_param$ZMu,
    TL_param$Fctrzn_Lrn_W0,
    TL_param$Fctrzn_Lrn_W
  )
  E_ZSqE_WSq <- lapply(
    views,
    E_ZSqE_WSq_update,
    TL_param$ZMu_0,
    TL_param$ZMuSq,
    TL_param$Fctrzn_Lrn_W0,
    TL_param$Fctrzn_Lrn_WSq
  )
  E_ZWSq <- lapply(
    views,
    E_ZWSq_update,
    E_ZE_W,
    TL_param$ZMuSq,
    E_Z_SqE_W_Sq,
    E_ZSqE_WSq
  )
  Zeta <- lapply(
    views,
    Zeta_calculation,
    likelihoods,
    E_ZWSq,
    E_ZE_W
  )
  Tau_m <- lapply(
    views,
    Tau_calculation,
    likelihoods,
    Zeta,
    TL_param$Tau
  )
  YGauss <- YGauss_calculation(
    "mRNA",
    likelihoods,
    TL_param$YTrg,
    Zeta,
    Tau_m,
    CenterTrg,
    PoisRateCstnt
  )
  expect_equal(YGauss, TL_param$YTrg$mRNA)
  YGauss <- YGauss_calculation(
    "miRNA",
    likelihoods,
    TL_param$YTrg,
    Zeta,
    Tau_m,
    CenterTrg,
    PoisRateCstnt
  )
  expect_equal(YGauss, (2 * TL_param$YTrg$miRNA - 1) / (2 * Tau_m$miRNA))
})

test_that("ZVar_calculation", {
  # ASK_DAVID - matrix looks strange - number inside are same
  ZVar <- ZVar_calculation(
    view = "mRNA",
    Tau = TL_param$Tau,
    Fctrzn_Lrn_WSq = TL_param$Fctrzn_Lrn_WSq
  )
  expect_equal(nrow(ZVar), nrow(TL_param$Tau$mRNA))
  expect_equal(ncol(ZVar), ncol(TL_param$Fctrzn_Lrn_WSq$mRNA))
  expect_equal(ZVar, TL_param$Tau$mRNA %*% TL_param$Fctrzn_Lrn_WSq$mRNA)
})

test_that("ZMu_calculation", {
  E_ZE_W <- list(
    "mRNA" = E_ZE_W_update(
      view = c("mRNA"),
      TL_param$ZMu_0,
      TL_param$ZMu,
      TL_param$Fctrzn_Lrn_W0,
      TL_param$Fctrzn_Lrn_W
    )
  )
  E_Z_SqE_W_Sq <- list(
    "mRNA" = E_Z_SqE_W_Sq_update(
      view = c("mRNA"),
      TL_param$ZMu_0,
      TL_param$ZMu,
      TL_param$Fctrzn_Lrn_W0,
      TL_param$Fctrzn_Lrn_W
    )
  )
  E_ZSqE_WSq <- list(
    "mRNA" = E_ZSqE_WSq_update(
      view = c("mRNA"),
      TL_param$ZMu_0,
      TL_param$ZMuSq,
      TL_param$Fctrzn_Lrn_W0,
      TL_param$Fctrzn_Lrn_WSq
    )
  )
  E_ZWSq <- list("mRNA" = E_ZWSq_update(
    view = c("mRNA"),
    E_ZE_W,
    TL_param$ZMuSq,
    E_Z_SqE_W_Sq,
    E_ZSqE_WSq
  ))
  Zeta <- Zeta_calculation("mRNA", likelihoods, E_ZWSq, E_ZE_W)
  Tau_m <- Tau_calculation("mRNA", likelihoods, Zeta, TL_param$Tau)
  YGauss <- list("mRNA" = YGauss_calculation(
    "mRNA",
    likelihoods,
    TL_param$YTrg,
    Zeta,
    Tau_m,
    CenterTrg,
    PoisRateCstnt
  ))
  ZMu_m <- ZMu_calculation(
    view = "mRNA",
    k = 1,
    TL_param$Fctrzn_Lrn_W,
    TL_param$Fctrzn_Lrn_W0,
    TL_param$Tau,
    TL_param$ZMu_0,
    TL_param$ZMu,
    YGauss
  )
  expect_equal(length(ZMu_m), nrow(TL_param$Tau$mRNA))
  expect_equal(length(ZMu_m), length(TL_param$ZMu_0))
})

test_that("ELBO_calculation", {
  likelihoods <- list(
    "mRNA" = "gaussian",
    "miRNA" = "gaussian",
    "DNAme" = "poisson",
    "SNV" = "bernoulli"
  )
  E_ZE_W <- lapply(
    views,
    E_ZE_W_update,
    TL_param$ZMu_0,
    TL_param$ZMu,
    TL_param$Fctrzn_Lrn_W0,
    TL_param$Fctrzn_Lrn_W
  )
  E_Z_SqE_W_Sq <- lapply(
    views,
    E_Z_SqE_W_Sq_update,
    TL_param$ZMu_0,
    TL_param$ZMu,
    TL_param$Fctrzn_Lrn_W0,
    TL_param$Fctrzn_Lrn_W
  )
  E_ZSqE_WSq <- lapply(
    views,
    E_ZSqE_WSq_update,
    TL_param$ZMu_0,
    TL_param$ZMuSq,
    TL_param$Fctrzn_Lrn_W0,
    TL_param$Fctrzn_Lrn_WSq
  )
  E_ZWSq <- lapply(
    views,
    E_ZWSq_update,
    E_ZE_W,
    TL_param$ZMuSq,
    E_Z_SqE_W_Sq,
    E_ZSqE_WSq
  )
  Zeta <- lapply(views, Zeta_calculation, likelihoods, E_ZWSq, E_ZE_W)
  Tau_m <- lapply(views, Tau_calculation, likelihoods, Zeta, TL_param$Tau)
  YGauss <- lapply(
    views,
    YGauss_calculation,
    likelihoods,
    TL_param$YTrg,
    Zeta,
    Tau_m,
    CenterTrg,
    PoisRateCstnt
  )
  ELBO_L <- ELBO_calculation(
    view = "mRNA",
    likelihoods,
    TL_param$Tau,
    TL_param$TauLn,
    E_ZWSq,
    E_ZE_W,
    Zeta,
    TL_param$YTrg,
    YGauss,
    PoisRateCstnt
  )
  expect_contains(ELBO_L, numeric())
  ELBO_L <- ELBO_calculation(
    view = "miRNA",
    likelihoods,
    TL_param$Tau,
    TL_param$TauLn,
    E_ZWSq,
    E_ZE_W,
    Zeta,
    TL_param$YTrg,
    YGauss,
    PoisRateCstnt
  )
  expect_contains(ELBO_L, numeric())
  ELBO_L <- ELBO_calculation(
    view = "DNAme",
    likelihoods,
    TL_param$Tau,
    TL_param$TauLn,
    E_ZWSq,
    E_ZE_W,
    Zeta,
    TL_param$YTrg,
    YGauss,
    PoisRateCstnt
  )
  expect_contains(ELBO_L, numeric())
})

test_that("E_ZE_W_update", {
  E_ZE_W <- E_ZE_W_update(
    view = "mRNA",
    TL_param$ZMu_0,
    TL_param$ZMu,
    TL_param$Fctrzn_Lrn_W0,
    TL_param$Fctrzn_Lrn_W
  )
  expect_equal(rownames(E_ZE_W), rownames(TL_param$ZMu))
  expect_equal(colnames(E_ZE_W), names(TL_param$Fctrzn_Lrn_W0$mRNA))
})

test_that("E_Z_SqE_W_Sq_update", {
  E_Z_SqE_W_Sq <- E_Z_SqE_W_Sq_update(
    view = "mRNA",
    TL_param$ZMu_0,
    TL_param$ZMu,
    TL_param$Fctrzn_Lrn_W0,
    TL_param$Fctrzn_Lrn_W
  )
  expect_equal(rownames(E_Z_SqE_W_Sq), rownames(TL_param$ZMu))
  expect_equal(
    colnames(E_Z_SqE_W_Sq),
    names(TL_param$Fctrzn_Lrn_W0$mRNA)
  )
})

test_that("E_ZSqE_WSq_update", {
  E_ZSqE_WSq <- E_ZSqE_WSq_update(
    view = "mRNA",
    TL_param$ZMu_0,
    TL_param$ZMuSq,
    TL_param$Fctrzn_Lrn_W0,
    TL_param$Fctrzn_Lrn_WSq
  )
  expect_equal(rownames(E_ZSqE_WSq), rownames(TL_param$ZMu_0))
  expect_equal(colnames(E_ZSqE_WSq), names(TL_param$Fctrzn_Lrn_W0$mRNA))
})

test_that("E_ZWSq_update", {
  E_ZE_W <- list(
    "mRNA" = E_ZE_W_update(
      view = "mRNA",
      TL_param$ZMu_0,
      TL_param$ZMu,
      TL_param$Fctrzn_Lrn_W0,
      TL_param$Fctrzn_Lrn_W
    )
  )
  E_Z_SqE_W_Sq <- list(
    "mRNA" = E_Z_SqE_W_Sq_update(
      view = "mRNA",
      TL_param$ZMu_0,
      TL_param$ZMu,
      TL_param$Fctrzn_Lrn_W0,
      TL_param$Fctrzn_Lrn_W
    )
  )
  E_ZSqE_WSq <- list(
    "mRNA" = E_ZSqE_WSq_update(
      view = "mRNA",
      TL_param$ZMu_0,
      TL_param$ZMuSq,
      TL_param$Fctrzn_Lrn_W0,
      TL_param$Fctrzn_Lrn_WSq
    )
  )
  E_ZWSq <- E_ZWSq_update(view = "mRNA",
                          E_ZE_W,
                          TL_param$ZMuSq,
                          E_Z_SqE_W_Sq,
                          E_ZSqE_WSq)
  expect_equal(dim(E_ZWSq), dim(E_Z_SqE_W_Sq$mRNA))
  expect_equal(dim(E_ZWSq), dim(E_ZSqE_WSq$mRNA))
})

test_that("VarExplFun", {
  likelihoods <- list("mRNA" = "gaussian",
                      "miRNA" = "gaussian",
                      "DNAme" = "poisson",
                      "SNV" = "bernoulli")
  YGauss <- YGauss_calculation(view = "mRNA",
                               likelihoods = likelihoods,
                               YTrg = YTrg,
                               Zeta = Zeta,
                               Tau = Tau,
                               CenterTrg = FALSE,
                               PoisRateCstnt = 0.0001)
  expect_equal(YGauss, YTrg$mRNA)
  YGauss <- YGauss_calculation(view = "DNAme",
                               likelihoods = likelihoods,
                               YTrg = YTrg,
                               Zeta = Zeta,
                               Tau = Tau,
                               CenterTrg = FALSE,
                               PoisRateCstnt = 0)
  expect_equal(dim(YGauss), dim(YTrg$DNAme))
  YGauss <- YGauss_calculation(view = "SNV",
                               likelihoods = likelihoods,
                               YTrg = YTrg,
                               Zeta = Zeta,
                               Tau = Tau,
                               CenterTrg = FALSE,
                               PoisRateCstnt = 0.0001)
  expect_equal(dim(YGauss), dim(YTrg$SNV))
})

test_that("transferLearning_function", {
  expected_names <- c("YTrg",
                      "Fctrzn_Lrn_W0",
                      "Fctrzn_Lrn_W",
                      "Fctrzn_Lrn_WSq",
                      "Tau",
                      "TauLn" )
  likelihoods <- c("mRNA" = "gaussian",
                   "miRNA" = "gaussian",
                   "DNAme" = "poisson",
                   "SNV" = "bernoulli")
  TL_param_exp <- initTransferLearningParamaters(
    YTrg = YTrg_prep,
    views = views,
    Fctrzn = Lrn_Fctrzn_init,
    likelihoods = likelihoods)
  expect_equal(names(TL_param_exp), expected_names)
  expect_equal(TL_param_exp$YTrg, YTrg)
  expect_equal(TL_param_exp$Fctrzn_Lrn_W0, Fctrzn_Lrn_W0)
  expect_equal(TL_param_exp$Fctrzn_Lrn_W, Fctrzn_Lrn_W)
  expect_equal(TL_param_exp$Fctrzn_Lrn_WSq, Fctrzn_Lrn_WSq)
  expect_equal(TL_param_exp$Tau, TL_param$Tau)
  expect_equal(TL_param_exp$TauLn, TL_param$TauLn)
})

