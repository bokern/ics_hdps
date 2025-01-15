# Author: Marleen Bokern
# Date: 10/2023
# Purpose: Master script for analysis of hospitalisations and deaths (wave 1) in COPD population
# suppress warnings

#set filepaths using filepaths.R
setwd(Github_folder)

#import and format analytic file
source("01_extract_codes.R")
setwd(HDPS_github)
#make baseline tables
Sys.time()
setwd(HDPS_github)
source("02_map_codes.R")
Sys.time()
setwd(HDPS_github)
source("02_map_codes_drugs.R")
Sys.time()
setwd(HDPS_github)
#source("02_map_codes_ever.R")
Sys.time()
setwd(HDPS_github)
source("02_map_codes_hes.R")
Sys.time()
setwd(HDPS_github)
source("03_assess_recurrence_multioutcome.R")
Sys.time()
setwd(HDPS_github)
source("03_assess_recurrence_multioutcome_no_triple.R")
Sys.time()

setwd(HDPS_github)
source("04_cov_weighting_parallel.R")
Sys.time()
setwd(HDPS_github)
source("05_effect_estimation.R")

Sys.time()
setwd(HDPS_github)
source("combined_forest.R")
Sys.time()
setwd(HDPS_github)
source("combined_forest_by_outcome.R")


setwd(paste0(HDPS_github, "/hdps-diagnostics"))
source("chapters_prevalence.R")
setwd(paste0(HDPS_github, "/hdps-diagnostics"))
source("chapters_prevalence_extra_codes.R")
setwd(paste0(HDPS_github, "/hdps-diagnostics"))
source("stdDiffsPredefinedPlusHDPS.R")

setwd(HDPS_folder)
source("03_assess_recurrence_multioutcome_no_triple_check.R")
Sys.time()
setwd(HDPS_folder)
source("04_cov_weighting_parallel_check.R")
Sys.time()
setwd(HDPS_folder)
source("05_effect_estimation_check.R")

