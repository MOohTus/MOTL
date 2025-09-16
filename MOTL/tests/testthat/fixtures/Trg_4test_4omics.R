## MT - 16/09/2025
##

## CREATE TARGET DATASET TOY DATA FOR TESTS

## LIBRARIES

## ENVIRONMENT
wd <- "/home/morgane/Documents/01_Projects/03_OtherProjects/05_David_Thesis/test_package_aout25/MOTL/tests/testthat/fixtures/"
Trg_dir <- file.path(wd, "Trg_4test_4omics")
Lrn_dir <- file.path(wd, "Lrn_4test_5000D")
Lrn_FctrnDir <- file.path(Lrn_dir, "Fctrzn_50K_01TH")

## LOAD FACTORIZED LEARNING DATA SET
Lrn_meta <- readRDS(file.path(Lrn_dir, "Lrn_meta.rds"))
InputModel = file.path(Lrn_FctrnDir, "Model.hdf5")
Fctrzn = load_model(file = InputModel)

## INIT VALUES FROM FACTORIZATION
Lrn_views <- Fctrzn@data_options$views
Lrn_likelihoods <- Fctrzn@model_options$likelihoods
Lrn_M <- Fctrzn@dimensions$M

Fctrzn@expectations[["Tau"]] = Tau_init(Lrn_views, Fctrzn, InputModel)
Fctrzn@expectations[["TauLn"]] = sapply(Lrn_views, TauLn_calculation, Lrn_likelihoods, Fctrzn, Lrn_FctrnDir)
Fctrzn@expectations[["WSq"]] = sapply(Lrn_views, WSq_calculation, Fctrzn, Lrn_FctrnDir)
Fctrzn@expectations[["W0"]] = sapply(Lrn_views, W0_calculation, CenterTrg, Fctrzn, Lrn_FctrnDir)

## INIT PARAMETERS
views = c("mRNA", "miRNA", "DNAme", "SNV")
Lrn_ftrs_mRNA <- Fctrzn@features_metadata[Fctrzn@features_metadata$view == "mRNA", "feature"]
Lrn_ftrs_miRNA <- Fctrzn@features_metadata[Fctrzn@features_metadata$view == "miRNA", "feature"]
Lrn_ftrs_DNAme <- Fctrzn@features_metadata[Fctrzn@features_metadata$view == "DNAme", "feature"]
Lrn_ftrs_SNV <- Fctrzn@features_metadata[Fctrzn@features_metadata$view == "SNV", "feature"]

## METADATA OF TARGET DATA
expdat_meta <- readRDS(file.path(Trg_dir, "expdat_meta.rds"))
smpls <- c("TCGA-D9-A4Z2-01A", "TCGA-EE-A2MN-06A", "TCGA-IB-A5ST-01A",
           "TCGA-XD-AAUG-01A", "TCGA-GN-A26D-06A", "TCGA-IB-7887-01A",
           "TCGA-FB-AAPP-01A", "TCGA-D9-A1X3-06A", "TCGA-2J-AABV-01A",
           "TCGA-FS-A1Z7-06A")
# smpls <- sample(expdat_meta$smpls, 10)

## mRNA
expdat_mRNA <- as.matrix(read.table(file = file.path(Trg_dir, "mRNA.csv"), sep = ","))
rownames(expdat_mRNA) <- expdat_meta$ftrs_mRNA
colnames(expdat_mRNA) <- expdat_meta$smpls
expdat_mRNA_ftrs <- sample(rownames(expdat_mRNA), 500)
expdat_mRNA_ftrs <- expdat_mRNA_ftrs[expdat_mRNA_ftrs %in% Lrn_ftrs_mRNA]
expdat_mRNA <- expdat_mRNA[expdat_mRNA_ftrs, smpls]
## miRNA
expdat_miRNA <- as.matrix(read.table(file = file.path(Trg_dir, "miRNA.csv"), sep = ","))
rownames(expdat_miRNA) <- expdat_meta$ftrs_miRNA
colnames(expdat_miRNA) <- expdat_meta$smpls
expdat_miRNA_ftrs <- sample(rownames(expdat_miRNA), 500)
expdat_miRNA_ftrs <- expdat_miRNA_ftrs[expdat_miRNA_ftrs %in% Lrn_ftrs_miRNA]
expdat_miRNA <- expdat_miRNA[expdat_miRNA_ftrs, smpls]
## DNAme
expdat_DNAme <- as.matrix(read.table(file = file.path(Trg_dir, "DNAme.csv"), sep = ","))
rownames(expdat_DNAme) <- expdat_meta$ftrs_DNAme
colnames(expdat_DNAme) <- expdat_meta$smpls
expdat_DNAme_ftrs <- sample(rownames(expdat_DNAme), 500)
expdat_DNAme_ftrs <- expdat_DNAme_ftrs[expdat_DNAme_ftrs %in% Lrn_ftrs_DNAme]
expdat_DNAme <- expdat_DNAme[expdat_DNAme_ftrs, smpls]
## SNV
expdat_SNV <- as.matrix(read.table(file = file.path(Trg_dir, "SNV.csv"), sep = ","))
rownames(expdat_SNV) <- expdat_meta$ftrs_SNV
colnames(expdat_SNV) <- expdat_meta$smpls
expdat_SNV_ftrs <- sample(rownames(expdat_SNV), 500)
expdat_SNV_ftrs <- expdat_SNV_ftrs[expdat_SNV_ftrs %in% Lrn_ftrs_SNV]
expdat_SNV <- expdat_SNV[expdat_SNV_ftrs, smpls]

## CREATE LIST OF TARGET DATASETS
YTrg_list <- list("mRNA" = SummarizedExperiment(list(expdat_mRNA)),
                  "miRNA" = SummarizedExperiment(list(expdat_miRNA)),
                  "DNAme" = SummarizedExperiment(list(expdat_DNAme)),
                  "SNV" = expdat_SNV)

## UPDATE METADATA OF TARGET DATA
expdat_meta$smpls <- smpls
expdat_meta$ftrs_mRNA <- rownames(expdat_mRNA)
expdat_meta$ftrs_miRNA <- rownames(expdat_miRNA)
expdat_meta$ftrs_DNAme <- rownames(expdat_DNAme)
expdat_meta$ftrs_SNV <- rownames(expdat_SNV)

## CREATE LIST OF TABLES WITH SAMPLE NAMES
brcds_SS <- list("brcds_mRNA_SS" = list(data.frame("brcds" = smpls)),
                 "brcds_miRNA_SS" = list(data.frame("brcds" = smpls)),
                 "brcds_DNAme_SS" = list(data.frame("brcds" = smpls)),
                 "brcds_SNV_SS" = list(data.frame("brcds" = smpls)))

## SAVE INTO RDS OBJECTF
Trg <- list("YTrg_list" = YTrg_list, "Trg_meta" = expdat_meta, "brcds_SS" = brcds_SS)
saveRDS(Trg, file = file.path(wd, "Trg_4test_4omics.rds"))

## PREPARE TARGET DATASET LIST
YTrg_prep <- TCGATargetDataPreparation(views = views,
                                       YTrgFull = YTrg_list,
                                       brcds_SS = brcds_SS,
                                       SS = 1,
                                       Fctrzn = Fctrzn,
                                       smpls = smpls,
                                       normalization = FALSE,
                                       expdat_meta_Lrn = Lrn_meta,
                                       transformation = FALSE)
saveRDS(YTrg_prep, file = file.path(wd, "Trg_prep_4test_4omics.rds"))

