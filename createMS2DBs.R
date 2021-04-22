library(wormLipidBlastR)
library(MsBackendMassbank)
library(Spectra)

lipid_file <- "F:/lipids.txt"

lipid_spectra <- predict_spectra(lipid_file)

export(lipid_spectra[[1]],
       MsBackendMassbank(),
       file = "F:/CelegansLipids_neg.mb")

export(lipid_spectra[[2]],
       MsBackendMassbank(),
       file = "F:/CelegansLipids_pos.mb")
