#############################
# Project: HDPS
# Programmer name: Marleen Bokern
# Date Started:  06/24
#############################  
# Map SNOMED TO ICD-10 for observation files
#############################
packages <- c("tidyverse", "arrow", "readstata13", "data.table", "autoCovariateSelection")
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
file_path <- file.path(Denominator, "CPRD Aurum", "Code browsers", "2022_05" , "CPRDAurumMedical.txt")
code_browser <- read.delim(file_path, sep = "\t", header = TRUE, colClasses = "character") %>% 
  select(MedCodeId, SnomedCTConceptId, Term)

#join code browser and code map on snomed
codes <- code_browser %>% 
  left_join(code_map, by = c("SnomedCTConceptId" = "snomed"))

#count NA values
codes %>% 
  summarise_all(~sum(is.na(.)))

#remove codes with NA values
codes_no_icd <- codes %>% filter(!is.na(icd10))

# Import data as a loop
obs_files <- c("Observations_for_hdps_1.parquet", "Observations_for_hdps_2.parquet")

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
  obs <- inner_join(obs, patid_list, by = "patid")
  obs <- obs[, c("patid", "eventdate", "medcodeid")]
  
  #merge with code lookup on medcodeid
  obs <- merge(obs, codes, by.x = "medcodeid", by.y = "MedCodeId", all.x = TRUE)

  #to create obs_ever, remove any duplicates in terms of patid and medcodeid to check for unmatched medcodes
  obs_ever <- obs[, c("patid", "medcodeid", "icd10")]
  obs_ever <- obs_ever[!duplicated(obs_ever[, c("patid", "medcodeid")]), ]
  
  if (is.null(obs_ever_unmatched)) {
    obs_ever_unmatched <- obs_ever[is.na(icd10)]
  } else {
    obs_ever_unmatched <- rbind(obs_ever_unmatched, obs_ever[is.na(icd10)])
  }
  
  #keep only deduplicated patid and icd10
  obs_ever <- obs[!duplicated(obs[, c("patid", "icd10")]), ]

  if (is.null(obs_ever_mapped)) {
    obs_ever_mapped <- obs_ever[!is.na(icd10)]
  } else {
    obs_ever_mapped <- rbind(obs_ever_mapped, obs_ever[!is.na(icd10)])
  }
  gc()
  
  #for obs, remove codes more than 1 year before the eventdate
  obs <- obs %>% filter(eventdate >= as.Date("2019-03-01"))
  if (is.null(obs_mapped)) {
    obs_mapped <- obs[!is.na(icd10)]
    obs_unmatched <- obs[is.na(icd10)]
  } else {
    obs_mapped <- rbind(obs_mapped, obs[!is.na(icd10)])
    obs_unmatched <- rbind(obs_unmatched, obs[is.na(icd10)])
  }
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
write_parquet(obs_ever_mapped, "observations_for_hdps_ever_mapped.parquet")
write_parquet(obs_ever_unmatched, "observations_for_hdps_ever_unmatched.parquet")

snomed_counts_unmatched <- obs_unmatched[, .N, by = SnomedCTConceptId]
snomed_counts_unmatched <- snomed_counts_unmatched %>% 
  left_join(codes, by = c("SnomedCTConceptId" = "SnomedCTConceptId")) %>% 
  arrange(desc(N))

#save top 100 unmatched codes.
#most unmatched codes are related to 
write.csv(snomed_counts_unmatched[1:100,], "snomed_counts_unmatched.csv")

