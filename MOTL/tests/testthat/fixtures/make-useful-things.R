## Lrn files from
## https://zenodo.org/records/10848217
LrnFctrnDir <- "../00_Ressources/00_DATA/Lrn_5000D_Fctrzn_100K_001TH/"
InputModel <- file.path(LrnFctrnDir, "Model.hdf5")
Fctrzn <- MOFA2::load_model(file = InputModel)
saveRDS(Fctrzn, file = "tests/testthat/fixtures/Lrn_5000D_Fctrzn_model.rds")

## Trg files
TrgDir <- "../00_Ressources/00_DATA/Trg_LAML_PAAD_SKCM_Full_5000D/"
expdat_mRNA <- as.matrix(read.table(file.path(TrgDir, "mRNA.csv"), sep = ","))
expdat_miRNA <- as.matrix(read.table(file.path(TrgDir, "miRNA.csv"), sep = ","))
expdat_DNAme <- as.matrix(read.table(file.path(TrgDir, "DNAme.csv"), sep = ","))
expdat_meta <- readRDS(file.path(TrgDir, "expdat_meta.rds"))
row.names(expdat_mRNA) <- expdat_meta$ftrs_mRNA
row.names(expdat_miRNA) <- expdat_meta$ftrs_miRNA
row.names(expdat_DNAme) <- expdat_meta$ftrs_DNAme
colnames(expdat_mRNA) <- expdat_meta$smpls
colnames(expdat_miRNA) <- expdat_meta$smpls
colnames(expdat_DNAme) <- expdat_meta$smpls

expdat_mRNA[c(1,5,89),] <- 0
expdat_miRNA[c(1,89),] <- 0
expdat_miRNA <- expdat_miRNA[,sample(1:ncol(expdat_miRNA))]
expdat_DNAme <- expdat_DNAme[,sample(1:ncol(expdat_DNAme))]

YTrg_list = list(mRNA = expdat_mRNA, miRNA = expdat_miRNA, DNAme = expdat_DNAme)
saveRDS(YTrg_list, file = "tests/testthat/fixtures/YTrg_list.rds")


expdat_mRNA[c(1:5), c(1:5)]
dim(expdat_mRNA)
expdat_miRNA[c(1:5), c(1:5)]
dim(expdat_miRNA)
expdat_DNAme[c(1:5), c(1:5)]
dim(expdat_DNAme)


## 
Fctrzn@expectations[["Tau"]] = Tau_init(viewsLrn, Fctrzn, InputModel)
Fctrzn@expectations[["TauLn"]] = sapply(viewsLrn, TauLn_calculation, likelihoodsLrn, Fctrzn, LrnFctrnDir)
Fctrzn@expectations[["WSq"]] = sapply(viewsLrn, WSq_calculation, Fctrzn, LrnFctrnDir)
Fctrzn@expectations[["W0"]] = sapply(viewsLrn, W0_calculation, CenterTrg, Fctrzn, LrnFctrnDir)
YTrg_prep = sapply(views, TargetDataPrefiltering, YTrg_list, Fctrzn, smpls)
YTrgFtrs <- lapply(YTrg_prep, rownames)
##
TL_param = initTransferLearningParamaters(YTrg_prep, views, expdat_meta_Lrn, Fctrzn, likelihoods)
TL_param$ZMu <- matrix(data = as.vector(colMeans(Fctrzn@expectations$Z$group0)),
                            nrow = dim(ZVar)[1], ncol = dim(ZVar)[2], byrow = TRUE)
rownames(TL_param$ZMu) <- smpls
colnames(TL_param$ZMu) <- colnames(TL_param$Fctrzn_Lrn_W[[1]])
TL_param$ZVar <- ZVar_calculation(view = "mRNA", Tau, Fctrzn_Lrn_WSq)
TL_param$ZMu_0 <- rep(1,dim(ZVar)[1])
TL_param$ZMuSq <- ZVar + ZMu^2


saveRDS(TL_param, file = "tests/testthat/fixtures/TL_param.rds")

