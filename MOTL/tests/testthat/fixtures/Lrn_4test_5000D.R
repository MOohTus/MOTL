## MT - 12/09/2025
##

## CREATE LEARNING DATASET TOY DATA FOR TESTS

## LIBRARIES
library("rjson")
library("MOFA2")

## ENVIRONEMENT
Seed = 1234567
mode(Seed) = 'integer'
set.seed(Seed)
wd <- "/home/morgane/Documents/01_Projects/03_OtherProjects/05_David_Thesis/test_package_aout25/MOTL/tests/testthat/fixtures/"

## LEARNING DATASET
Prjct_dir <- "Lrn_4test_5000D/"
Lrn_meta <- fromJSON(file = file.path(wd, Prjct_dir, "expdat_meta.json"))
Lrn_mRNA_ftrs <- sample(Lrn_meta$ftrs_mRNA, 2000)
Lrn_miRNA_ftrs <- sample(Lrn_meta$ftrs_miRNA, 500)
Lrn_DNAme_ftrs <- sample(Lrn_meta$ftrs_DNAme, 2000)
Lrn_SNV_ftrs <- sample(Lrn_meta$ftrs_SNV, 2000)
Lrn_smpls <- sample(Lrn_meta$smpls, 500)

## mRNA
Lrn_mRNA <- read.table(file = file.path(wd, Prjct_dir, "mRNA.csv"), sep = ",")
rownames(Lrn_mRNA) <- Lrn_meta$ftrs_mRNA
colnames(Lrn_mRNA) <- Lrn_meta$smpls
Lrn_mRNA <- Lrn_mRNA[Lrn_mRNA_ftrs, Lrn_smpls]
write.table(x = Lrn_mRNA, file = file.path(wd, Prjct_dir, "mRNA_4test.csv"), sep = ",", row.names = FALSE, col.name = FALSE)

## miRNA
Lrn_miRNA <- read.table(file = file.path(wd, Prjct_dir, "miRNA.csv"), sep = ",")
rownames(Lrn_miRNA) <- Lrn_meta$ftrs_miRNA
colnames(Lrn_miRNA) <- Lrn_meta$smpls
Lrn_miRNA <- Lrn_miRNA[Lrn_miRNA_ftrs, Lrn_smpls]
write.table(x = Lrn_miRNA, file = file.path(wd, Prjct_dir, "miRNA_4test.csv"), sep = ",", row.names = FALSE, col.name = FALSE)

## DNAme
Lrn_DNAme <- read.table(file = file.path(wd, Prjct_dir, "DNAme.csv"), sep = ",")
rownames(Lrn_DNAme) <- Lrn_meta$ftrs_DNAme
colnames(Lrn_DNAme) <- Lrn_meta$smpls
Lrn_DNAme <- Lrn_DNAme[Lrn_DNAme_ftrs, Lrn_smpls]
write.table(x = Lrn_DNAme, file = file.path(wd, Prjct_dir, "DNAme_4test.csv"), sep = ",", row.names = FALSE, col.name = FALSE)

## SNV
Lrn_SNV <- read.table(file = file.path(wd, Prjct_dir, "SNV.csv"), sep = ",")
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
FctrznDir = file.path(wd, Prjct_dir, 'Fctrzn_50K_01TH')
Lrn_meta = readRDS(file.path(wd, Prjct_dir, 'Lrn_meta.rds'))
InputModel = file.path(FctrznDir, "Model.hdf5")
Fctrzn = load_model(file = InputModel)

intercepts_calculation(expdat_meta = Lrn_meta,
                       Fctrzn = Fctrzn,
                       FctrznDir = FctrznDir,
                       ExpDataDir = file.path(wd, Prjct_dir),
                       Seed = Seed)
