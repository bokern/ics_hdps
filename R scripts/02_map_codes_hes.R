#############################
# Project: HDPS
# Programmer name: Marleen Bokern
# Date Started:  06/24
#############################  
# Look at HES codes from year before index
# Process: trim icd-10 code to 3 characters. Count number of codes per patient, so that the result is a df with variables patid, code, and counts of that individual code per patient. From that, calculate median, p75 and max number of codes per patient for each icd10. If median is equal to 1, set to NA, if 75 is equal to median, set to p75 to NA. From that, create a df with 1 row per patient and up to 3 columns per icd10 code. The columns are numbered A00_1, A00_2 and A00_3. They are binary, indicating whether the patient has ever had the code, whether they have had the code more than the median number of times, and whether they have had the code more than the 75th percentile number of times. If the median or p75 are set to NA, the corresponding columns should be skipped. 
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