# ==============================================================================
# parameter
# 
# Here all parameter for the execution of the annotation workflow are set.
# ==============================================================================
# set project folder -----------------------------------------------------------
project_root_pos <- NA
project_root_neg <- NA

# set tolerances ---------------------------------------------------------------
tolerance = 0.005
ppm = 0
rtTolerance = 0.1

# set allowed adducts ----------------------------------------------------------
adducts_pos <- c("[M+H]+", "[M+Na]+", "[M+NH4]+")
adducts_neg <- c("[M-H]-", "[M+FA-H]-")

# ==============================================================================
# Setup 
# 
# Different libraries are required for the execution of the annotation workflow.
# ==============================================================================
# load required libraries ------------------------------------------------------
library(Spectra)
library(MsBackendMgf)
library(MsBackendMassbank)
library(MsBackendMsp)
library(MsCoreUtils)
library(tidyverse)
library(genedataRutils)

#register(bpstart(SnowParam(6, progressbar = TRUE)))
script_root <- dirname(rstudioapi::getActiveDocumentContext()$path)

# ==============================================================================
# Perform annotation for positive mode data
# 
# This block performs all annotation steps on the positive mode data.
# ==============================================================================
if(!is.na(project_root_pos)) {
  
  # set working directory to positive project folder ---------------------------
  setwd(project_root_pos)
  
  # ============================================================================
  # project sanity checks
  # 
  # Here different sanity checks are performed to ensure that the project folder
  # contains the correct information and data required for the annotation 
  # process. Minimally required are following folders and files
  # - file: xyz_Cluster.gda (Genedata .gda file containing all detected cluster)
  # - file: xyz_Peak.gda (Genedata .gda file containing all detected peaks)
  # - folder: MS1_libraries (this folder contains potential MS1 libraries)
  # - folder: MS2_libraries (this folder contains potential MS2 libraries in MSP
  # or MassBank record format)
  #
  # Missing result folders are generated if required.
  # ============================================================================
  # check if data is in place --------------------------------------------------
  # folder
  if(!dir.exists("MS2_consolidated")) {stop("Missing MS2 spectra")}
  if(!dir.exists("MS1_libraries")) {stop("Missing MS1 libraries folder")}
  if(!dir.exists("MS2_libraries")) {stop("Missing MS2 libraries folder")}
  
  # files
  if(length(list.files(project_root_pos,
                       pattern = "_Cluster.gda$",
                       full.names = TRUE)) == 0) {stop("Missing _Cluster.gda file")}
  
  if(length(list.files(project_root_pos,
                       pattern = "_Peak.gda$",
                       full.names = TRUE)) == 0) {stop("Missing _Peak.gda file")}
  
  if(length(list.files(paste0(project_root_pos, "/MS1_libraries"),
                       full.names = TRUE)) == 0) {stop("No MS1 libraries defined")}
  
  if(length(list.files(paste0(project_root_pos, "/MS2_libraries"),
                       full.names = TRUE)) == 0) {stop("No MS2 libraries defined")}
  
  # create folders for results -------------------------------------------------
  if(!dir.exists("MS1_consolidated_R")) {dir.create("MS1_consolidated_R")}
  if(!dir.exists("MS1_results")) {dir.create("MS1_results")}
  if(!dir.exists("MS2_consolidated_R")) {dir.create("MS2_consolidated_R")}
  if(!dir.exists("MS2_results")) {dir.create("MS2_results")}
  if(!dir.exists("Sirius_ms")) {dir.create("Sirius_ms")}
  if(!dir.exists("Sirius_ms_consolidated")) {dir.create("Sirius_ms_consolidated")}
  
  # ============================================================================
  # MS1 processing
  #
  # In this block first the isotopic pattern are reconstructed from the cluster 
  # and peak information provided in the two .gda files. All reconstructed
  # isotope pattern are afterwards combined into a consensus isotope pattern.
  # ============================================================================
  # load MS1 data --------------------------------------------------------------
  ms1_cluster <- readGda(list.files(project_root,
                                    pattern = "_Cluster.gda$",
                                    full.names = TRUE))
  
  ms1_peak <- readGda(list.files(project_root,
                                 pattern = "_Peak.gda$",
                                 full.names = TRUE))
  
  # reconstruct MS1 spectra ----------------------------------------------------
  ms1_spectra <- reconstructIsoPattern(ms1_peak, ms1_cluster)
  
  # consolidation of MS1 spectra -----------------------------------------------
  # TODO add final consolidation strategy
  ms1_spectra_comb <- ms1_spectra %>% 
    combineSpectra(f = ms1_spectra$CLUSTER_ID,
                   #p = ms1_spectra$CLUSTER_ID,
                   intensityFun = base::mean,
                   mzFun = base::mean,
                   tolerance = tolerance,
                   minProp = 0.9,
                   peaks = "intersect")
  
  # ============================================================================
  # MS2 processing
  #
  # MS2 spectra are read using the Spectra package and used CLUSTER_IDs from the
  # MS1 data is added in order to allow mapping between the different MS level.
  # After this spectra are consolidated two combine all spectra from the same
  # MS1 feature into a single consensus spectrum.
  # ============================================================================
  # load MS2 data --------------------------------------------------------------
  mgf_files <- list.files(paste0(project_root, "/MS2_consolidated"),
                          pattern = ".mgf$",
                          full.names = TRUE)
  
  ms2_spectra <- Spectra(mgf_files,
                         source = MsBackendMgf(),
                         backend = MsBackendDataFrame())
  
  # add MS1 id -----------------------------------------------------------------
  ms2_spectra <- ms2AddId(ms1_cluster, ms2_spectra)
  ms2_spectra <- setBackend(ms2_spectra, MsBackendDataFrame())
  
  # consolidation of MS2 spectra -----------------------------------------------
  # TODO add final consolidation strategy
  ms2_spectra_comb <- ms2_spectra %>% 
    combineSpectra(f = ms2_spectra$CLUSTER_ID,
                   #p = ms2_spectra$CLUSTER_ID,
                   intensityFun = base::mean,
                   mzFun = base::mean,
                   tolerance = tolerance,
                   minProp = 0.9,
                   peaks = "intersect")
  
  # ============================================================================
  # export of processed spectra
  # 
  # All processed spectra are stored in MGF files in distinct folders. The
  # combined spectra are as well exported as Sirius .ms files to enable 
  # annotation using the external Sirius tool.
  # ============================================================================
  # export spectra to .mgf file ------------------------------------------------
  export(ms1_spectra,
         MsBackendMgf(),
         file = "MS1_consolidated_R/ms1_spectra.mgf")
  export(ms1_spectra_comb,
         MsBackendMgf(),
         file = "MS1_consolidated_R/ms1_spectra_combined.mgf")
  export(ms2_spectra,
         MsBackendMgf(),
         file = "MS2_consolidated_R/ms2_spectra.mgf")
  export(ms2_spectra_comb,
         MsBackendMgf(),
         file = "MS2_consolidated_R/ms2_spectra_combined.mgf")
  
  # write Sirius files ---------------------------------------------------------
  # write to contain all MS2 spectra
  writeSiriusFile(ms1_spectra_comb,
                  ms2_spectra,
                  folder = "Sirius_ms")
  
  # write to contain only single consolidated MS2 spectra
  writeSiriusFile(ms1_spectra_comb,
                  ms2_spectra_comb,
                  folder = "Sirius_ms_consolidated")
  
  # ============================================================================
  # MS1 annotation
  #
  # Annotation is performed by comparing a list of precalculated m/z values of 
  # different metabolites or lipids. If the name of the library contains in the 
  # name "inhouse", then also the retention times are compared.  Retention times
  # in the dataset as well in the library have to be on the same scale, e.g.
  # both in minutes or seconds. Libraries are filtered to only contain the
  # adducts defined in parameters above.
  # ============================================================================
  # read MS1 libraries ---------------------------------------------------------
  ms1_libraries <- list.files(paste0(project_root, 
                                     "/MS1_libraries"),
                              full.names = TRUE)
  
  for(ms1_library in ms1_libraries) {
    
    if(grepl("inhouse", ms1_library)) {
      
      ms1_library <- ms1_library %>% filter(adduct %in% adducts_pos)
      
      if(nrow(ms1_library) > 0) {
        
        result <- annotateMz()
        write_tsv(result, paste0("MS1_results/result_",
                                 basename(ms1_library),
                                 ".tsv"))
      }
      
    } else {
      
    }
  }
  
  # ============================================================================
  # MS2 annotation
  #
  # MS2 annotation is performed by spectra matching against reference libraries.
  # If the name of the library contains in the name "inhouse", then also the 
  # retention times are compared. Retention times in the dataset as well in the 
  # library have to be on the same scale, e.g. both in minutes or seconds. 
  # Libraries are filtered to only contain the adducts defined in parameters
  # above.
  # ============================================================================
  # read MS2 libraries ---------------------------------------------------------
  ms2_libraries <- list.files(paste0(project_root, 
                                     "/MS2_libraries"),
                              full.names = TRUE)
  
  for(ms2_library in ms2_libraries) {
    
    # read library
    if(grepl(".msp$", ms2_library)) {
      
      # read msp
      ms2_library_spectra <- Spectra(ms2_library,
                                     source = MsBackendMsp(),
                                     backend = MsBackendDataFrame())
      
    } else if(grepl(".mb$", ms2_library)) {
      
      # read massbank
      ms2_library_spectra <- Spectra(ms2_library,
                                     source = MsBackendMassBank(),
                                     backend = MsBackendDataFrame())
      
    }
    
    # filter based on adducts
    ms2_library_spectra <- ms2_library_spectra[ms2_library_spectra$adduct %in% adducts_pos]
    
    if(length(ms2_library_spectra) > 0) {
      
      result <-compareSpectraLibrary(ms2_spectra_comb,
                                     ms2_library_spectra,
                                     tolerance = 0.005)
      
      write_tsv(result, paste0("MS2_results/result_",
                               basename(ms2_library),
                               "_combinedSpectra.tsv"))
      
      result <-compareSpectraLibrary(ms2_spectra,
                                     ms2_library_spectra,
                                     tolerance = 0.005)
      
      write_tsv(result, paste0("MS2_results/result_",
                               basename(ms2_library),
                               "_singleSpectra.tsv"))
      
    }
  }
}
