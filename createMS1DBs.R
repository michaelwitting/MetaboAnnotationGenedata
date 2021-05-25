# load required libraries ------------------------------------------------------
library(tidyverse)
library(dplyr)
library(readxl)
library(reshape2)
library(MetaboCoreUtils)

# read sheets from Excel file --------------------------------------------------
db_xlsx <- "DBs/CelegansMS1DBs.xlsx"

dbs <- excel_sheets(db_xlsx)
dbs <- dbs[!dbs %in% c("template")]

metabo_adducts_pos <- c("[M+H]+", "[M+Na]+")
metabo_adducts_neg <- c("[M-H]-")
lipid_adducts_pos <- c("[M+H]+", "[M+Na]+", "[M+NH4]+")
lipid_adducts_neg <- c("[M-H]-", "[M+CHO2]-")

# iterate over adducts ---------------------------------------------------------
for(db in dbs) {
  
  # read data ------------------------------------------------------------------
  data <- read_xlsx(db_xlsx,
                    db,
                    skip = 2,
                    col_names = TRUE) %>% 
    filter_all(any_vars(!is.na(.)))

  date <- read_xlsx(db_xlsx,
                    db,
                    range = "B1",
                    col_names = FALSE) %>% 
    as.character()
  
  type <- read_xlsx(db_xlsx,
                    db,
                    range = "B2",
                    col_names = FALSE) %>% 
    as.character()
  
  # create positive data -------------------------------------------------------
  if(type == "metabo") {
    
    adductMzs <- mass2mz(data$exactmass, metabo_adducts_pos)
    
  } else if(type == "lipid") {
    
    adductMzs <- mass2mz(data$exactmass, lipid_adducts_pos)
    
  } else {
    
    message("No type defined, defaulting to metabolomics adducts")
    adductMzs <- mass2mz(data$exactmass, metabo_adducts_pos)
    
  }
  
  pos_mz <- melt(cbind.data.frame(data,
                                 adductMzs,
                                 stringsAsFactors = FALSE),
                id.vars = colnames(data),
                value.name = "mz",
                variable.name = "adduct")
  
  # create positive data -------------------------------------------------------
  if(type == "metabo") {
    
    adductMzs <- mass2mz(data$exactmass, metabo_adducts_neg)
    
  } else if(type == "lipid") {
    
    adductMzs <- mass2mz(data$exactmass, lipid_adducts_neg)
    
  } else {
    
    message("No type defined, defaulting to metabolomics adducts")
    adductMzs <- mass2mz(data$exactmass, metabo_adducts_neg)
    
  }
  
  neg_mz <- melt(cbind.data.frame(data,
                                  adductMzs,
                                  stringsAsFactors = FALSE),
                 id.vars = colnames(data),
                 value.name = "mz",
                 variable.name = "adduct")
  
  # write data to file
  write_tsv(pos_mz, file = paste0("DBs/MS1/", type, "_", db, "_", date, "_pos.tsv"))
  write_tsv(neg_mz, file = paste0("DBs/MS1/", type, "_", db, "_", date, "_neg.tsv"))
  
}
