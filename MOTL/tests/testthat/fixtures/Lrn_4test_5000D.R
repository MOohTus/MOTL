## MT - 12/09/2025
##

## CREATE LEARNING DATASET TOY DATA FOR TESTS

## expdat_meta.json is coming from https://zenodo.org/records/10848217
## mRNA_4test.csv / miRNA_4test.csv / DNAme_4test.csv / SNV_4test.csv were renamed without _4test

## LIBRARIES
library("MOFA2")
library("rjson")

## ENVIRONEMENT
Seed <- 1234567
mode(Seed) <- "integer"
set.seed(Seed)
wd <- "/home/morgane/Documents/01_Projects/03_OtherProjects/05_David_Thesis/test_package_aout25/00_Ressources/01_dataOfPackage/"
data_rep <- "00_originalData/"

## LEARNING DATASET
Prjct_dir <- "Lrn_4test_5000D/"
Lrn_meta <- readRDS(file = file.path(wd, data_rep, "expdat_meta_Lrn.rds"))
Lrn_mRNA_ftrs <- sample(Lrn_meta$ftrs_mRNA, 1000)
Lrn_miRNA_ftrs <- sample(Lrn_meta$ftrs_miRNA, 250)
Lrn_DNAme_ftrs <- sample(Lrn_meta$ftrs_DNAme, 1000)
Lrn_SNV_ftrs <- sample(Lrn_meta$ftrs_SNV, 1000)
Lrn_smpls <- sample(Lrn_meta$smpls, 250)

## mRNA
Lrn_mRNA <- read.table(file = file.path(wd, data_rep, "mRNA.csv"), sep = ",")
rownames(Lrn_mRNA) <- Lrn_meta$ftrs_mRNA
colnames(Lrn_mRNA) <- Lrn_meta$smpls
Lrn_mRNA <- Lrn_mRNA[Lrn_mRNA_ftrs, Lrn_smpls]
write.table(x = Lrn_mRNA, file = file.path(wd, Prjct_dir, "mRNA_4test.csv"), sep = ",", row.names = FALSE, col.name = FALSE)

## miRNA
Lrn_miRNA <- read.table(file = file.path(wd, data_rep, "miRNA.csv"), sep = ",")
rownames(Lrn_miRNA) <- Lrn_meta$ftrs_miRNA
colnames(Lrn_miRNA) <- Lrn_meta$smpls
Lrn_miRNA <- Lrn_miRNA[Lrn_miRNA_ftrs, Lrn_smpls]
write.table(x = Lrn_miRNA, file = file.path(wd, Prjct_dir, "miRNA_4test.csv"), sep = ",", row.names = FALSE, col.name = FALSE)

## DNAme
Lrn_DNAme <- read.table(file = file.path(wd, data_rep, "DNAme.csv"), sep = ",")
rownames(Lrn_DNAme) <- Lrn_meta$ftrs_DNAme
colnames(Lrn_DNAme) <- Lrn_meta$smpls
Lrn_DNAme <- Lrn_DNAme[Lrn_DNAme_ftrs, Lrn_smpls]
write.table(x = Lrn_DNAme, file = file.path(wd, Prjct_dir, "DNAme_4test.csv"), sep = ",", row.names = FALSE, col.name = FALSE)

## SNV
Lrn_SNV <- read.table(file = file.path(wd, data_rep, "SNV.csv"), sep = ",")
rownames(Lrn_SNV) <- Lrn_meta$ftrs_SNV
colnames(Lrn_SNV) <- Lrn_meta$smpls
Lrn_SNV <- Lrn_SNV[Lrn_SNV_ftrs, Lrn_smpls]
write.table(x = Lrn_SNV, file = file.path(wd, Prjct_dir, "SNV_4test.csv"), sep = ",", row.names = FALSE, col.name = FALSE)

## METADATA
Lrn_meta$ftrs_mRNA <- rownames(Lrn_mRNA)
Lrn_meta$ftrs_miRNA <- rownames(Lrn_miRNA)
Lrn_meta$ftrs_DNAme <- rownames(Lrn_DNAme)
Lrn_meta$ftrs_SNV <- rownames(Lrn_SNV)
Lrn_meta$smpls <- Lrn_smpls
Lrn_meta.json <- toJSON(Lrn_meta)
write(Lrn_meta.json, file.path(wd, Prjct_dir, "Lrn_meta.json"))
saveRDS(Lrn_meta, file.path(wd, Prjct_dir, "Lrn_meta.rds"))

## CREATE MODEL
FctrznDir <- file.path(wd, Prjct_dir, "Fctrzn_50K_01TH")
Lrn_meta <- readRDS(file.path(wd, Prjct_dir, "Lrn_meta.rds"))
InputModel <- file.path(FctrznDir, "Model.hdf5")
Fctrzn <- load_model(file = InputModel)
Lrn <- list("Fctrzn" = Fctrzn)
saveRDS(Lrn, file.path(FctrznDir, "Lrn_Fctrzn_4omics.rds"))

intercepts_calculation(
    expdat_meta = Lrn_meta,
    Fctrzn = Fctrzn,
    FctrznDir = FctrznDir,
    ExpDataDir = file.path(wd, Prjct_dir),
    Seed = Seed
)
