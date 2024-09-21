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
n_drug_files <- 49
  
#read patient file copd_Patid_list_included_all.dta
patid_list <- read.dta13("copd_Patid_list_included_all.dta")

setwd(Copd_aurum_extract)

# define list of observation files as all files in the directory that have the string "copd_Extract_Observation_" and ".dta"in their name
obs_files <- list.files(pattern = "copd_Extract_Observation_.*\\.dta")

read_observation_files <- function(file, all_obs) {
  print(file)
  obs <- read.dta13(file)
  #drop numrangelow and numrangehigh variables
  obs <- obs %>% select(-numrangelow, -numrangehigh, -obstypeid, -numunitid, -value)
  
  obs <- inner_join(obs, patid_list, by = "patid")
  #drop observations before 01 march 2020
  obs <- obs[obs$eventdate <= as.Date("2020-03-01"), ]
  
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
    gc()
  }
  
  # Save the batch as a Parquet file
  output_file <- paste0("Observations_for_hdps_", i, ".parquet")
  write_parquet(batch_obs, output_file)
}

end_time1 <- Sys.time()
run_time1 <- end_time1 - start_time1
print(run_time1)

rm(all_obs)
gc()
# Prescription files ------------------------------------------------------
start_time2 <- Sys.time()

# define list of drug files as all files in the directory that have the string "copd_Extract_DrugIssue" and ".dta"in their name
drug_files <- list.files(pattern = "copd_Extract_DrugIssue.*\\.dta")
 file <- "copd_Extract_DrugIssue_1.dta"
read_drug_files <- function(file, all_drugs) {
  print(file)
  drugs <- read.dta13(file)
  
  drugs <- inner_join(drugs, patid_list, by = "patid")
  drugs <- drugs[drugs$issuedate <= as.Date("2020-03-01"), ]
  drugs <- drugs[drugs$issuedate >= as.Date("2019-03-01"),]
  drugs <- drugs %>% select(-estnhscost, -quantity, -quantunitid, -duration, -dosageid)
  
  all_drugs <- rbind(all_drugs, drugs)
  
  return(all_drugs)
}

# assuming all_obs is defined elsewhere in your script
all_drugs <- data.frame()
start_time <- Sys.time()

# Split the files into two batches
batch_size <- 20
file_batches <- split(drug_files, ceiling(seq_along(drug_files) / batch_size))

for (i in seq_along(file_batches)) {
  batch_drugs <- data.frame() # Initialize an empty data frame for the current batch
  
  # Loop through the files within the current batch
  for (file in file_batches[[i]]) {
    # Read the observation files for the current file and append to batch_obs
    batch_drugs <- read_drug_files(file, batch_drugs)
    gc()
  }
  
  # Save the batch as a Parquet file
  output_file <- paste0("Drugs_for_hdps_", i, ".parquet")
  write_parquet(batch_drugs, output_file)
}

end_time2 <- Sys.time()
run_time2 <- end_time2 - start_time2

rm(all_drugs)
gc()
# HES ---------------------------------------------------------------------

setwd(Linkage)

# define list of drug files as all files in the directory that have the string "copd_Extract_DrugIssue" and ".dta"in their name
hes_files <- list.files(pattern = "hes_diagnosis_epi.*\\.dta")

read_hes_files <- function(file, all_hes) {
  print(file)
  
  # Read in dta file
  hes <- read.dta13(file)
  
  hes <- inner_join(hes, patid_list, by = "patid")
  
  hes$epistart <- as.Date(hes$epistart, format = "%d/%m/%Y")
  
  hes <- hes[hes$epistart <= as.Date("2020-03-01") & hes$epistart >= as.Date("2019-03-01"), ]
  
  all_hes <- rbind(all_hes, hes)
  
  return(all_hes)
}

start_time3 <- Sys.time()

#loop through files and append them

# Initialize all_hes as an empty data frame
all_hes <- data.frame()

# Loop through files and append them
for (file in hes_files) {
  all_hes <- read_hes_files(file, all_hes)
  
  # Write the combined data frame to a Parquet file
  output_file <- paste0("HES_for_hdps.parquet")
  write_parquet(all_hes, output_file)
  gc()
}

end_time3 <- Sys.time()
run_time3 <- end_time3 - start_time3
