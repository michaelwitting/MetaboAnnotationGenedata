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

# ==============================================================================
# parameter
# 
# Here all parameter for the execution of the annotation workflow are set.
# ==============================================================================
# set project folder -----------------------------------------------------------
project_root_pos <- "K:/Data_processed/20210421_PathogenMetabolomicsNew/20210219_CelegansPathogenMicrobiome_RP_pos"
project_root_neg <- "K:/Data_processed/20210421_PathogenMetabolomicsNew/20210222_CelegansPathogenMicrobiome_RP_neg"

# set tolerances ---------------------------------------------------------------
tolerance = 0.005
ppm = 0
rtTolerance = 0.1

# set allowed adducts ----------------------------------------------------------
# The adduct information has to conform the notation used in MetaboCoreUtils
# check MetaboCoreUtils::adductNames("positive") or 
# MetaboCoreUtils::adductNames("negative")
adducts_pos <- c("[M+H]+", "[M+Na]+", "[M+NH4]+")
adducts_neg <- c("[M-H]-", "[M+CHO2]-")

# filtering function -----------------------------------------------------------
filter <- TRUE
my_filter <- function(x) {
  x > max(x, na.rm = TRUE) / 10
}

# some housekeeping
script_root <- dirname(rstudioapi::getActiveDocumentContext()$path)
parameter <- ls()

# parallelization of processes, if wanted --------------------------------------
register(bpstart(SnowParam(12, progressbar = TRUE)))

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
  ms1_cluster <- readGda(list.files(project_root_pos,
                                    pattern = "_Cluster.gda$",
                                    full.names = TRUE))
  
  ms1_peak <- readGda(list.files(project_root_pos,
                                 pattern = "_Peak.gda$",
                                 full.names = TRUE))
  
  # reconstruct MS1 spectra ----------------------------------------------------
  ms1_spectra <- reconstructIsoPattern(ms1_peak, ms1_cluster)
  
  # consolidation of MS1 spectra -----------------------------------------------
  # TODO add final consolidation strategy
  ms1_spectra_comb <- ms1_spectra %>% 
    combineSpectra(f = ms1_spectra$CLUSTER_ID,
                   p = ms1_spectra$CLUSTER_ID,
                   intensityFun = base::sum,
                   mzFun = base::mean,
                   tolerance = tolerance,
                   ppm = ppm,
                   minProp = 0.9,
                   peaks = "intersect",
                   weighted = TRUE)
  
  # ============================================================================
  # MS2 processing
  #
  # MS2 spectra are read using the Spectra package and used CLUSTER_IDs from the
  # MS1 data is added in order to allow mapping between the different MS level.
  # After this spectra are consolidated two combine all spectra from the same
  # MS1 feature into a single consensus spectrum.
  # ============================================================================
  # load MS2 data --------------------------------------------------------------
  mgf_files <- list.files(paste0(project_root_pos, "/MS2_consolidated"),
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
                   p = ms2_spectra$CLUSTER_ID,
                   intensityFun = base::sum,
                   mzFun = base::mean,
                   tolerance = tolerance,
                   ppm = ppm,
                   minProp = 0.9,
                   peaks = "intersect",
                   weighted = TRUE)
  
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
  
  # write overview on MS2 data -------------------------------------------------
  table(ms2_spectra$CLUSTER_ID) %>%
    as.data.frame() %>% as_tibble() %>% 
    write_tsv("MS2_consolidated_R_/ms2_spectra_overview.tsv")
  
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
  ms1_libraries <- list.files(paste0(project_root_pos, 
                                     "/MS1_libraries"),
                              full.names = TRUE)
  
  for(ms1_library in ms1_libraries) {
    
    ms1_library_data <- read_tsv(ms1_library)
    
    # correction of adducts
    # no correction required at the moment
    
    ms1_library_data <- ms1_library_data %>% filter(adduct %in% adducts_pos)
    
    
    if(nrow(ms1_library_data) > 0 && grepl("inhouse", ms1_library)) {

        
        result <- annotateMz(ms1_cluster,
                             ms1_library_data,
                             tolerance = tolerance,
                             ppm = ppm,
                             adducts = adducts_pos)
        
        write_tsv(result, paste0("MS1_results/result_",
                                 basename(ms1_library),
                                 ".tsv"))

      
    } else if(nrow(ms1_library_data) > 0) {

        
        result <- annotateMz(ms1_cluster,
                             ms1_library_data,
                             tolerance = tolerance,
                             ppm = ppm,
                             adducts = adducts_pos)
        
        write_tsv(result, paste0("MS1_results/result_",
                                 basename(ms1_library),
                                 ".tsv"))
      
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
  ms2_libraries <- list.files(paste0(project_root_pos, 
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
    
    # correction of adducts
    # no correction required in pos mode so far
    
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
      
      if(filter) {
        
        result <-compareSpectraLibrary(filterIntensity(ms2_spectra, my_filter),
                                       ms2_library_spectra,
                                       tolerance = 0.005)
        
        write_tsv(result, paste0("MS2_results/result_",
                                 basename(ms2_library),
                                 "_singleSpectra_filtered.tsv"))
        
      }
    }
  }
  
  # ============================================================================
  # Aggregation of results
  #
  # Several parts of the workflow create results that have one ore more result
  # for one single peak/cluster etc. This block is aggregating the results for
  # a better overview. In case of MS1 data, results are concatenated by a "_".
  # In regards to MS2 data multiple MS2 spectra with results may exist for one
  # MS1 feature. Aggregation is performed by combining results for the same MS1
  # and compound annotation. The aggregated data contains then min, max, mean,
  # median for the forward and reverse score.
  # ============================================================================
  # aggregate MS1 results ------------------------------------------------------
  # MS1 results files
  ms1_result_files <- list.files(paste0(project_root_pos,
                                        "/MS1_results"),
                                 full.names = TRUE)
  
  # aggregate all result files
  for(ms1_result_file in ms1_result_files) {
    
    ms1_data <- read_tsv(ms1_result_file)
    
    ms1_data_cond <- aggregate(ms1_data,
                               list(ms1_data$id),
                               function(x) {
                                 paste0(unique(x), collapse ="_")
                                 })
    
    write_tsv(ms1_data_cond, paste0("MS1_results/",
                                    basename(ms1_result_file),
                                    "_aggregated.tsv"))
  }
  
  
  # aggregate MS2 results ------------------------------------------------------
  # MS2 results files
  ms2_result_files <- list.files(paste0(project_root_pos,
                                        "/MS2_results"),
                                 full.names = TRUE,
                                 pattern = ".*singleSpectra.tsv|.*singleSpectra_filtered.tsv")
  
  # aggrate all result files
  for(ms2_result_file in ms2_result_files) {
    
    ms2_data <- read_tsv(ms2_result_file)
    
    ms2_data_cond <- ms2_data %>%
      group_by(CLUSTER_ID, lib_name) %>% 
      summarise(no_of_Spectra = n(),
                forwardMin = min(forward),
                forwardMax = max(forward),
                forwardMean = mean(forward),
                forwardMedian = median(forward),
                backwardMin = min(backward),
                backwardMax = max(backward),
                backwardMean = mean(backward),
                backwardMedian = median(backward),
                countMin = min(count),
                countMax = max(count),
                countMean = mean(count),
                countMedian = median(count))
    
    write_tsv(ms2_data_cond, paste0("MS2_results/",
                                    basename(ms2_result_file),
                                    "_aggregated.tsv"))
    
  }
}

  
# clear environment to avoid collisions
rm(list=ls()[!ls() %in% parameter])
parameter <- ls()

# ==============================================================================
# Perform annotation for negative mode data
# 
# This block performs all annotation steps on the negative mode data.
# ==============================================================================
if(!is.na(project_root_neg)) {
  
  # set working directory to negative project folder ---------------------------
  setwd(project_root_neg)
  
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
  if(length(list.files(project_root_neg,
                       pattern = "_Cluster.gda$",
                       full.names = TRUE)) == 0) {stop("Missing _Cluster.gda file")}
  
  if(length(list.files(project_root_neg,
                       pattern = "_Peak.gda$",
                       full.names = TRUE)) == 0) {stop("Missing _Peak.gda file")}
  
  if(length(list.files(paste0(project_root_neg, "/MS1_libraries"),
                       full.names = TRUE)) == 0) {stop("No MS1 libraries defined")}
  
  if(length(list.files(paste0(project_root_neg, "/MS2_libraries"),
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
  ms1_cluster <- readGda(list.files(project_root_neg,
                                    pattern = "_Cluster.gda$",
                                    full.names = TRUE))
  
  ms1_peak <- readGda(list.files(project_root_neg,
                                 pattern = "_Peak.gda$",
                                 full.names = TRUE))
  
  # reconstruct MS1 spectra ----------------------------------------------------
  ms1_spectra <- reconstructIsoPattern(ms1_peak, ms1_cluster)
  
  # consolidation of MS1 spectra -----------------------------------------------
  # TODO add final consolidation strategy
  ms1_spectra_comb <- ms1_spectra %>% 
    combineSpectra(f = ms1_spectra$CLUSTER_ID,
                   p = ms1_spectra$CLUSTER_ID,
                   intensityFun = base::sum,
                   mzFun = base::mean,
                   tolerance = tolerance,
                   ppm = ppm,
                   minProp = 0.9,
                   peaks = "intersect",
                   weighted = TRUE)
  
  # ============================================================================
  # MS2 processing
  #
  # MS2 spectra are read using the Spectra package and used CLUSTER_IDs from the
  # MS1 data is added in order to allow mapping between the different MS level.
  # After this spectra are consolidated two combine all spectra from the same
  # MS1 feature into a single consensus spectrum.
  # ============================================================================
  # load MS2 data --------------------------------------------------------------
  mgf_files <- list.files(paste0(project_root_neg, "/MS2_consolidated"),
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
                   p = ms2_spectra$CLUSTER_ID,
                   intensityFun = base::sum,
                   mzFun = base::mean,
                   tolerance = tolerance,
                   ppm = ppm,
                   minProp = 0.9,
                   peaks = "intersect",
                   weighted = TRUE)
  

  
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
  
  # write overview on MS2 data -------------------------------------------------
  table(ms2_spectra$CLUSTER_ID) %>%
    as.data.frame() %>% as_tibble() %>% 
    write_tsv("MS2_consolidated_R_/ms2_spectra_overview.tsv")
  
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
  ms1_libraries <- list.files(paste0(project_root_neg, 
                                     "/MS1_libraries"),
                              full.names = TRUE)
  
  for(ms1_library in ms1_libraries) {
    
    ms1_library_data <- read_tsv(ms1_library)
    
    # correction of adducts
    ms1_library_data$adduct <- str_replace(ms1_library_data$adduct, "\\[M\\+FA\\-H\\]\\-", "[M+CHO2]-")
    ms1_library_data$adduct <- str_replace(ms1_library_data$adduct, "\\[M\\+Hac\\-H\\]\\-", "[M+C2H3O2]-")
    ms1_library_data$adduct <- str_replace(ms1_library_data$adduct, "\\[M\\+HAc\\-H\\]\\-", "[M+C2H3O2]-")
    
    ms1_library_data <- ms1_library_data %>% filter(adduct %in% adducts_neg)
    
    if(nrow(ms1_library_data) > 0 && grepl("inhouse", ms1_library)) {
      
      
      result <- annotateMz(ms1_cluster,
                           ms1_library_data,
                           tolerance = tolerance,
                           ppm = ppm,
                           adducts = adducts_neg)
      
      write_tsv(result, paste0("MS1_results/result_",
                               basename(ms1_library),
                               ".tsv"))
      
      
    } else if(nrow(ms1_library_data) > 0) {
      
      
      result <- annotateMz(ms1_cluster,
                           ms1_library_data,
                           tolerance = tolerance,
                           ppm = ppm,
                           adducts = adducts_neg)
      
      write_tsv(result, paste0("MS1_results/result_",
                               basename(ms1_library),
                               ".tsv"))
      
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
  ms2_libraries <- list.files(paste0(project_root_neg, 
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
    
    # correction of adducts
    ms2_library_spectra$adduct <- str_replace(ms2_library_spectra$adduct, "\\[M\\+FA\\-H\\]\\-", "[M+CHO2]-")
    ms2_library_spectra$adduct <- str_replace(ms2_library_spectra$adduct, "\\[M\\+Hac\\-H\\]\\-", "[M+C2H3O2]-")
    ms2_library_spectra$adduct <- str_replace(ms2_library_spectra$adduct, "\\[M\\+HAc\\-H\\]\\-", "[M+C2H3O2]-")
    
    # filter based on adducts
    ms2_library_spectra <- ms2_library_spectra[ms2_library_spectra$adduct %in% adducts_neg]
    
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
      if(filter) {
        
        result <-compareSpectraLibrary(filterIntensity(ms2_spectra, my_filter),
                                       ms2_library_spectra,
                                       tolerance = 0.005)
        
        write_tsv(result, paste0("MS2_results/result_",
                                 basename(ms2_library),
                                 "_singleSpectra_filtered.tsv"))
        
      }
    }
  }
  
  # ============================================================================
  # Aggregation of results
  #
  # Several parts of the workflow create results that have one ore more result
  # for one single peak/cluster etc. This block is aggregating the results for
  # a better overview. In case of MS1 data, results are concatenated by a "_".
  # In regards to MS2 data multiple MS2 spectra with results may exist for one
  # MS1 feature. Aggregation is performed by combining results for the same MS1
  # and compound annotation. The aggregated data contains then min, max, mean,
  # median for the forward and reverse score.
  # ============================================================================
  # aggregate MS1 results ------------------------------------------------------
  # MS1 results files
  ms1_result_files <- list.files(paste0(project_root_neg,
                                        "/MS1_results"),
                                 full.names = TRUE)
  
  # aggregate all result files
  for(ms1_result_file in ms1_result_files) {
    
    ms1_data <- read_tsv(ms1_result_file)
    
    ms1_data_cond <- aggregate(ms1_data,
                               list(ms1_data$id),
                               function(x) {
                                 paste0(unique(x), collapse ="_")
                               })
    
    write_tsv(ms1_data_cond, paste0("MS1_results/",
                                    basename(ms1_result_file),
                                    "_aggregated.tsv"))
  }
  
  
  # aggregate MS2 results ------------------------------------------------------
  # MS2 results files
  ms2_result_files <- list.files(paste0(project_root_neg,
                                        "/MS2_results"),
                                 full.names = TRUE,
                                 pattern = ".*singleSpectra.tsv|.*singleSpectra_filtered.tsv")
  
  # aggrate all result files
  for(ms2_result_file in ms2_result_files) {
    
    ms2_data <- read_tsv(ms2_result_file)
    
    ms2_data_cond <- ms2_data %>%
      group_by(CLUSTER_ID, lib_name) %>% 
      summarise(no_of_Spectra = n(),
                forwardMin = min(forward),
                forwardMax = max(forward),
                forwardMean = mean(forward),
                forwardMedian = median(forward),
                backwardMin = min(backward),
                backwardMax = max(backward),
                backwardMean = mean(backward),
                backwardMedian = median(backward),
                countMin = min(count),
                countMax = max(count),
                countMean = mean(count),
                countMedian = median(count))
    
    write_tsv(ms2_data_cond, paste0("MS2_results/",
                                    basename(ms2_result_file),
                                    "_aggregated.tsv"))
    
  }
}

# clear environment to avoid collisions
rm(list=ls()[!ls() %in% parameter])
parameter <- ls()
  
# ==============================================================================
# Perform matching of positive and negative mode data
# 
# This block performs matching of positive and negative mode data based on RT.
# ==============================================================================
if(!is.na(project_root_pos) & !is.na(project_root_neg)) {
  
  # load MS1 data (positive) ---------------------------------------------------
  ms1_cluster_pos <- readGda(list.files(project_root_pos,
                                        pattern = "_Cluster.gda$",
                                        full.names = TRUE))
  
  # load MS1 data (negative) ---------------------------------------------------
  ms1_cluster_neg <- readGda(list.files(project_root_neg,
                                        pattern = "_Cluster.gda$",
                                        full.names = TRUE))
  
  match_df <- matchIonMode(ms1_cluster_pos,
                           ms1_cluster_neg,
                           pos_adducts = adducts_pos,
                           neg_adducts = adducts_neg,
                           tolerance = tolerance,
                           ppm = ppm,
                           rtOffset = 0,
                           rtimeTolerance = rtTolerance)
  
  
}

# stop parallel backend --------------------------------------------------------
bpstop()