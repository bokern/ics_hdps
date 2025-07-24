#######################################
# Project: HDPS
# Programmer name: Marleen Bokern
# Date Started:  06/24
#######################################
# Map SNOMED TO ICD-10 for observation files
# Description:
# This script prepares primary care observation data for HDPS analysis by mapping SNOMED CT codes to ICD-10 codes. 
# The steps include:
#
# 1. imports patient ID lists.
# 2. Reads SNOMED-to-ICD10 mapping files and CPRD Aurum code browser.
# 3. Merging the CPRD browser with the SNOMED-to-ICD10 map.
# 4. Importing and processing CPRD observational data stored as Parquet files:
#    - Filtering to events before 01-Mar-2020 and after 01-Mar-2019 (lookback window)
#    - Keeping only observations with valid `medcodeid`s and matching patient IDs
#    - Mapping `medcodeid` to ICD-10 via SNOMED using the merged lookup
#    - Removing COVID-19-specific codes (U07) as this denotes the outcome
#    - Separately saving matched and unmatched observations 
#
# 5. Creating an "ever-mapped" dataset to capture each patient's unique set of ICD-10 codes.
# 6. Tabulating and exporting unmatched SNOMED CT codes to identify mapping gaps.
#
# Input files:
# - Patient inclusion list: 'copd_Patid_list_included_all.dta'
# - SNOMED-to-ICD10 map: 'snomed_to_icd10.dta'
# - Code browser: 'CPRDAurumMedical.txt'
# - Observations Parquet files: 'Observations_for_hdps_1.parquet', 'Observations_for_hdps_2.parquet'
#
# Output files:
# - Mapped observations: 'observations_for_hdps_mapped.parquet'
# - Unmatched observations: 'observations_for_hdps_unmatched.parquet'
# - Ever-mapped observations: 'observations_for_hdps_ever_mapped.parquet'
# - Unmatched SNOMED CT codes summary: 'snomed_counts_unmatched.csv'

#######################################

packages <- c("tidyverse",
              "arrow",
              "readstata13",
              "data.table")

installed_packages <- packages %in% rownames(installed.packages())
if (any(installed_packages == FALSE)) {
  install.packages(packages[!installed_packages])
}

invisible(lapply(packages, function(pkg) {
  suppressMessages(library(pkg, character.only = TRUE, verbose = FALSE))
}))

setwd(Datadir_copd)
patid_list <- read.dta13("copd_Patid_list_included_all.dta")

#unzip files
setwd(Copd_aurum_extract)

#import mapping file from stata as string variables
file_path <- file.path(Projectdir_stata, "snomed_to_icd10.dta")
code_map <- read.dta13(file_path, convert.factors = FALSE)

# import CPRD code_browser as string variables
file_path <- file.path(Denominator,
                       "CPRD Aurum",
                       "Code browsers",
                       "2022_05" ,
                       "CPRDAurumMedical.txt")
code_browser <- read.delim(file_path,
                           sep = "\t",
                           header = TRUE,
                           colClasses = "character") %>%
  dplyr::select(MedCodeId, SnomedCTConceptId, Term)

#join code browser and code map on snomed
codes <- code_browser %>%
  left_join(code_map, by = c("SnomedCTConceptId" = "snomed"))

#count NA values
codes %>%
  summarise_all( ~ sum(is.na(.)))

#remove codes with NA values
codes_no_icd <- codes %>% filter(!is.na(icd10))

# Import data as a loop
obs_files <- c("Observations_for_hdps_1.parquet",
               "Observations_for_hdps_2.parquet")

obs_mapped <- NULL
obs_unmatched <- NULL

obs_ever_mapped <- NULL
obs_ever_unmatched <- NULL

for (file in obs_files) {

  print(file)
  obs <- arrow::read_parquet(file)
  obs <- setDT(obs)

  #remove codes after index date or that don't belong to patients in the patid list
  obs <- obs %>% filter(eventdate <= as.Date("2020-03-01"))
  
  gc()
  #remove any empty medcodeids
  obs <- obs[!is.na(medcodeid) & medcodeid != "", ]
  obs <- inner_join(obs, patid_list, by = "patid")
  obs <- obs[, c("patid", "eventdate", "medcodeid")]
  
  #merge with code lookup on medcodeid
  obs <- merge(obs,
               codes,
               by.x = "medcodeid",
               by.y = "MedCodeId",
               all.x = TRUE)
  
  #remove U07 codes from obs
  obs <- obs[!obs$icd10 %in% c("U07"), ]
  
  #to create obs_ever, remove any duplicates in terms of patid and medcodeid to check for unmatched medcodes
  obs_ever <- obs[, c("patid", "medcodeid", "icd10")]
  obs_ever <- obs_ever[!duplicated(obs_ever[, c("patid", "medcodeid")]), ]
  
  #keep only deduplicated patid and icd10
  obs_ever <- obs[!duplicated(obs[, c("patid", "icd10")]), ]
  
  obs_ever <- obs_ever[!is.na(icd10) & icd10 != ""]
  obs_ever_mapped <- rbind(obs_ever_mapped, obs_ever[!is.na(icd10) & icd10 != ""])
  
  #for obs, remove codes more than 1 year before the eventdate
  obs <- obs %>% filter(eventdate >= as.Date("2019-03-01"))
  
  # Remove unmatched codes from obs
  obs_mapped <- rbindlist(list(obs_mapped, obs[!is.na(icd10) & icd10 != ""]), use.names = TRUE, fill = TRUE)
  
  # Take all unmatched codes from obs
  obs_unmatched <- rbindlist(list(obs_unmatched, obs[is.na(icd10) | icd10 == ""]), use.names = TRUE, fill = TRUE)
  
  rm(obs)
  rm(obs_ever)
  gc()
  
}

#for all datasets, keep patid and icd10, renamed as code
obs_mapped <- obs_mapped[, c("patid", "icd10")]
obs_mapped <- obs_mapped %>% rename(code = icd10)
obs_ever_mapped <- obs_ever_mapped[, c("patid", "icd10")]
obs_ever_mapped <- obs_ever_mapped %>% rename(code = icd10)

setwd(Datadir_copd)
#save obs_matched and obs_unmatched
write_parquet(obs_mapped, "observations_for_hdps_mapped.parquet")
write_parquet(obs_unmatched, "observations_for_hdps_unmatched.parquet")
write_parquet(obs_ever_mapped,
              "observations_for_hdps_ever_mapped.parquet")

snomed_counts_unmatched <- obs_unmatched[, .N, by = SnomedCTConceptId]
snomed_counts_unmatched <- snomed_counts_unmatched %>%
  left_join(codes, by = c("SnomedCTConceptId" = "SnomedCTConceptId")) %>%
  arrange(desc(N))

#save top 100 unmatched codes
write.csv(snomed_counts_unmatched[1:100, ], 
          "snomed_counts_unmatched.csv", 
          row.names = FALSE, 
          quote = TRUE)
