#############################
# Project: HDPS
# Programmer name: Marleen Bokern
# Date Started:  06/24
#############################
# Description:
# This script processes Hospital Episode Statistics (HES) data for HDPS feature generation.
# Steps:
# 1. Imports the HES data from a parquet file.
# 2. Trims the ICD-10 code to the first 3 characters.
# 3. Saves the cleaned HES data to a new parquet file for use in downstream analysis.
#############################

packages <- c("tidyverse", "arrow")
installed_packages <- packages %in% rownames(installed.packages())
if (any(installed_packages == FALSE)) {
  install.packages(packages[!installed_packages])
}

invisible(lapply(packages, function(pkg) {
  suppressMessages(library(pkg, character.only = TRUE, verbose = FALSE))
}))


#unzip files
setwd(Copd_aurum_extract)

#import mapping file from stata as string variables
file_path <- file.path(Linkage, "HES_for_hdps.parquet")
hes <- read_parquet(file_path)
hes <- setDT(hes)

#remove spno, epikey and epiend from hes
hes <- hes %>% 
  dplyr::select(-spno, -epikey, -epiend)

#create new variable denoting icd, trimmed to 3 characters
hes <- hes %>% 
  mutate(icd_trim = substr(icd, 1, 3)) %>%
  dplyr::select(-icd, -icdx, -d_order)

#rename icd_trim to code
setnames(hes, "icd_trim", "code")

file_path <- file.path(Datadir_copd, "hes_for_hdps_mapped.parquet")
write_parquet(hes, file_path)