## MT - 12/09/2025
##

## LIBRARIES
library("rjson")

## ENVIRONEMENT
Seed = 1234567
mode(Seed) = 'integer'
set.seed(Seed)
wd <- "/home/morgane/Documents/01_Projects/03_OtherProjects/05_David_Thesis/test_package_aout25/00_Ressources/00_DATA/4Package/"

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

## CREATE MODEL WITH makeFiles4Package.py
FctrznDir = file.path(wd, Prjct_dir,'Fctrzn_50K_01TH')
Lrn_meta = readRDS(file.path(wd, Prjct_dir,'Lrn_meta.rds'))
InputModel = file.path(FctrznDir,"Model.hdf5")
Fctrzn = load_model(file = InputModel)

intercepts_calculation(expdat_meta = Lrn_meta, 
                       Fctrzn = Fctrzn,
                       FctrznDir = FctrznDir, 
                       ExpDataDir = file.path(wd, Prjct_dir), 
                       Seed = Seed)


## MAKE INPUT FILES FOR MOTL PACKAGE

## ENVIRONMENT 
wd <- "/home/morgane/Documents/01_Projects/03_OtherProjects/05_David_Thesis/test_package_aout25/00_Ressources/00_DATA/Trg_PAAD_SKCM_Full_5000D"

## INIT PARAMETERS
views = c("mRNA", "miRNA", "DNAme", "SNV")

## METADATA OF TARGET DATA
expdat_meta <- readRDS(file.path(wd, "expdat_meta.rds"))
smpls <- sample(expdat_meta$smpls, 10)

## mRNA
expdat_mRNA <- read.table(file = file.path(wd, "mRNA.csv"), sep = ",")
rownames(expdat_mRNA) <- expdat_meta$ftrs_mRNA
colnames(expdat_mRNA) <- expdat_meta$smpls
expdat_mRNA <- expdat_mRNA[c(1:(nrow(expdat_mRNA)/10)),smpls]
## miRNA
expdat_miRNA <- read.table(file = file.path(wd, "miRNA.csv"), sep = ",")
rownames(expdat_miRNA) <- expdat_meta$ftrs_miRNA
colnames(expdat_miRNA) <- expdat_meta$smpls
expdat_miRNA <- expdat_miRNA[c(1:round(nrow(expdat_miRNA)/10)),smpls]
## DNAme
expdat_DNAme <- read.table(file = file.path(wd, "DNAme.csv"), sep = ",")
rownames(expdat_DNAme) <- expdat_meta$ftrs_DNAme
colnames(expdat_DNAme) <- expdat_meta$smpls
expdat_DNAme <- expdat_DNAme[c(1:(nrow(expdat_DNAme)/10)),smpls]
## SNV
expdat_SNV <- read.table(file = file.path(wd, "SNV.csv"), sep = ",")
rownames(expdat_SNV) <- expdat_meta$ftrs_SNV
colnames(expdat_SNV) <- expdat_meta$smpls
expdat_SNV <- expdat_SNV[c(1:(nrow(expdat_SNV)/10)),smpls]

## CREATE LIST OF TARGET DATASETS
YTrg_list <- list("mRNA" = expdat_mRNA,
                  "miRNA" = expdat_miRNA,
                  "DNAme" = expdat_DNAme,
                  "SNV" = expdat_SNV)

## PREPARE TARGET DATASET LIST
YTrg_prep <- TCGATargetDataPreparation(views = views, 
                                       YTrgFull = YTrg_list, 
                                       brcds_SS = smpls, 
                                       SS = 1, Fctrzn = , 
                                       smpls = smpls, 
                                       normalization = FALSE, 
                                       expdat_meta_Lrn = , 
                                       transformation = FALSE)




TCGATargetDataPreparation(views = views, YTrgFull = YTrg_list, brcds_SS = )

