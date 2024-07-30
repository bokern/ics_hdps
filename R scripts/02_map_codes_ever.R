#############################
# Project: HDPS
# Programmer name: Marleen Bokern
# Date Started:  05/24
#############################  
# 
#############################

packages <- c("tidyverse", "MetBrewer", "arrow", "readstata13", "bench")
installed_packages <- packages %in% rownames(installed.packages())
if (any(installed_packages == FALSE)) {
  install.packages(packages[!installed_packages])
}

invisible(lapply(packages, function(pkg) {
  suppressMessages(library(pkg, character.only = TRUE, verbose = FALSE))
}))

palette <- met.brewer("Cassatt2")

setwd(Datadir_copd)

n_obs_files <- 48

#read patient file copd_Patid_list_included_all.dta
patid_list <- read.dta13("copd_Patid_list_included_all.dta")

setwd(Copd_aurum_extract)

# define list of observation files as all files in the directory that have the string "copd_Extract_Observation_" and ".dta"in their name
obs_files <- list.files(pattern = "copd_Extract_Observation_.*\\.dta")

read_observation_files <- function(file, all_obs) {
  print(file)
  gc()
  obs <- read.dta13(file)
  #drop obs after 01 march 2020
  obs <- obs %>% filter(eventdate <= as.Date("2020-03-01"))
  #drop numrangelow and numrangehigh variables
  obs <- obs %>% select(patid, medcodeid)
  
  #drop observations that are duplicated in terms of patid and icd10
  obs <- obs[!duplicated(obs[, c("patid", "medcodeid")]), ]
  
  obs <- inner_join(obs, patid_list, by = "patid")
  
  all_obs <- rbind(all_obs, obs)
  
  return(all_obs)
}

# assuming all_obs is defined elsewhere in your script
all_obs <- data.frame()
start_time1 <- Sys.time()

# Split the files into two batches
batch_size <- 24
file_batches <- split(obs_files, ceiling(seq_along(obs_files) / batch_size))

for (i in seq_along(file_batches)) {
  batch_obs <- data.frame() # Initialize an empty data frame for the current batch
  
  # Loop through the files within the current batch
  for (file in file_batches[[i]]) {
    # Read the observation files for the current file and append to batch_obs
    batch_obs <- read_observation_files(file, batch_obs)
    
  }
  gc()
  # Save the batch as a Parquet file
  output_file <- paste0("Observations_for_hdps_ever", i, ".parquet")
  write_parquet(batch_obs, output_file)
}

end_time1 <- Sys.time()
run_time1 <- end_time1 - start_time1
print(run_time1)

rm(all_obs)
gc()

setwd(Copd_aurum_extract)

#read in first observations_for_hdps_ever file
obs1 <- read_parquet("Observations_for_hdps_ever1.parquet")

#append the other files
obs2 <- read_parquet(paste0("Observations_for_hdps_ever2.parquet"))
obs_all <- rbind(obs1, obs2)
#obs_all <- obs_all[1:100000,]

#map this to icd10
#import mapping file from stata as string variables
file_path <- file.path(Projectdir_stata, "snomed_to_icd10.dta")
code_map <- read.dta13(file_path, convert.factors = FALSE)

#code_browser as string variables
file_path <- file.path(Denominator, "CPRD Aurum", "Code browsers", "2022_05" , "CPRDAurumMedical.txt")
code_browser <- read.delim(file_path, sep = "\t", header = TRUE, colClasses = "character") %>% 
  select(MedCodeId, SnomedCTConceptId)
setwd(Copd_aurum_extract)

#join code browser and code map
codes <- code_browser %>% 
  left_join(code_map, by = c("SnomedCTConceptId" = "snomed"))
#remove snomed
codes <- codes %>% select(-SnomedCTConceptId)

#left join on medcodeid = MedCodeId
obs_all <- left_join(obs_all, codes, by = c("medcodeid" = "MedCodeId"))
#remove duplicates in terms of patid and icd10
obs_all <- obs_all[!duplicated(obs_all[, c("patid", "icd10")]), ]
#remove icd NAs
obs_all <- obs_all %>% filter(!is.na(icd10))

obs_all <- obs_all %>% rename(code = icd10)
obs_all <- obs_all %>% select(-medcodeid)
#write to parquet
write_parquet(obs_all, "observations_for_hdps_ever_mapped.parquet")







#create a dataset with patid as column one, and then each icd-10 code as a column, with 0s and 1s indicating whether the patient had the code
obs_all_wide <- obs_all %>% 
  select(patid, icd10) %>% 
  mutate(value = 1) %>% 
  pivot_wider(names_from = icd10, values_from = value, values_fill = 0)
 



#load patient baseline characteristics for main analysis
setwd(Datadir_copd)
pat_baseline <- read_parquet("copd_wave1_60d.parquet")
#keep patid, treatgroup, treat, covid_hes_present, covid_death_present, any_death_present
pat_baseline <- pat_baseline %>% select(patid, treatgroup, covid_hes_present, covid_death_present, any_death_present)

#merge in those baseline characteristics to the obs_all_wide dataset
obs_all_wide <- inner_join(obs_all_wide, pat_baseline, by = "patid")
#rearrage column so that patid, treatgroup, covid_hes_present, covid_death_present, any_death_present, treat are first
obs_all_wide <- obs_all_wide %>% select(patid, covid_hes_present, covid_death_present, any_death_present, treatgroup, everything())
#remove rows with treatgroup = NA
obs_all_wide <- obs_all_wide[!is.na(obs_all_wide$treatgroup), ]

# Function to calculate prevalence, risk ratio, and bias
calculate_measures <- function(data, outcome, treatment, covariate) {
  # Extract necessary data
  #
  e1c1 <- sum(data[[covariate]][data[[treatment]] == 1 & data[[outcome]] == 1], na.rm = TRUE)
  e1 <- sum(data[[treatment]] == 1, na.rm = TRUE)
  e0c1 <- sum(data[[covariate]][data[[treatment]] == 0 & data[[outcome]] == 1], na.rm = TRUE)
  e0 <- sum(data[[treatment]] == 0, na.rm = TRUE)
  d1c1 <- sum(data[[outcome]][data[[covariate]] == 1 & data[[treatment]] == 1], na.rm = TRUE)
  c1 <- sum(data[[covariate]] == 1, na.rm = TRUE)
  d1c0 <- sum(data[[outcome]][data[[covariate]] == 0 & data[[treatment]] == 1], na.rm = TRUE)
  c0 <- sum(data[[covariate]] == 0, na.rm = TRUE)
  
  # Calculate measures
  pc1 <- e1c1 / e1
  pc0 <- e0c1 / e0
  rr_ce <- pc1 / pc0
  rr_ce[rr_ce == 0] <- NA # Replace 0 with NA
  rr_cd <- (d1c1 / c1) / (d1c0 / c0)
  bias <- (pc1 * (rr_cd - 1) + 1) / (pc0 * (rr_cd - 1) + 1)
  abs_log_bias <- abs(log(bias))
  
  return(data.frame(prevalence_treatment = pc1, prevalence_control = pc0, rr_ce = rr_ce, rr_cd = rr_cd, bias = bias, abs_log_bias = abs_log_bias))
}
#check if there are any character columns
obs_all_wide %>% select_if(is.character) %>% colnames()
# Get the column names for ICD-10 codes
icd10_cols <- names(obs_all_wide)[-(1:5)]  # Exclude the first 5 columns
#ensure that 
# Apply the function to each ICD-10 code column
measures <- lapply(icd10_cols, function(covariate) {
  calculate_measures(obs_all_wide, "covid_death_present", "treatgroup", covariate)
})

# Combine the results into a single data frame
measures_df <- do.call(rbind, measures)
rownames(measures_df) <- icd10_cols

#convert to numeric, not scientific notation
measures_df <- measures_df %>% mutate_all(as.numeric) %>% format(digits = 3)

