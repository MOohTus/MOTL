## Lrn files from
## https://zenodo.org/records/10848217
LrnFctrnDir = "../00_Ressources/00_DATA/Lrn_5000D_Fctrzn_100K_001TH/"
InputModel = file.path(LrnFctrnDir, "Model.hdf5")
Fctrzn = MOFA2::load_model(file = InputModel)
saveRDS(Fctrzn, file = "tests/testthat/fixtures/Lrn_5000D_Fctrzn_model.rds")
