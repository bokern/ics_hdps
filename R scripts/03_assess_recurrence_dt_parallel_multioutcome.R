# -----------------------------------------------------------------------------
# PROGRAM NAME:  501_HDPS_ACNU_Proc
# PROJECT:      xzx26005
# AUTHOR:        John Tazare, adapted by Marleen Bokern
# DATE CREATED:   03 Sept 2020
# NOTES:       Performs hd-PS procedure for the ACNU
#
# REQUIRES:
# -----------------------------------------------------------------------------

packages <- c(
  "tidyverse",
  "arrow",
  "dbplyr",
  "lubridate",
  "scales",
  "DBI",
  "sparklyr",
  "data.table",
  "furrr",
  "future.apply",
  "progress")

installed_packages <- packages %in% rownames(installed.packages())
if (any(installed_packages == FALSE)) {
  install.packages(packages[!installed_packages])
}

invisible(lapply(packages, function(pkg) {
  suppressMessages(library(pkg, character.only = TRUE, verbose = FALSE))
}))

set.seed(123)
# Flag to enable debugging; if TRUE, script will only be run on 1,000 patients
# selected at random
debugMode = F

sample_size <- 10000
exclude_triple <- F
gc()

if (exclude_triple == TRUE) {
  output_ext <- "_no_triple"
} else {
  output_ext <- ""
}

# Project Setup----------------------------------------------------------------
setwd(Datadir_copd)

# Assign priority of select to dplyr
select <- dplyr::select

# Load in patient summary
pat_summary_data <- read_parquet("copd_wave1_60d.parquet") %>%
  mutate_at(c("patid"), as.character) %>%
  # keep only people whose treatgroup is not na
  filter(!is.na(treatgroup))
if (debugMode == TRUE) {
  # Load in patient summary
  pat_summary_data <- pat_summary_data %>%
    filter(patid %in% sample(unique(patid), sample_size))
}

if (exclude_triple == TRUE) {
  pat_summary_data <- pat_summary_data %>%
    filter(baseline_triple == 0)
}
    
# Clinical
dimension <- "observations"
clinicalDim <- read_parquet(paste0(dimension, "_for_HDPS_mapped.parquet")) %>%
  mutate_at(c("patid"), as.character) %>%
  mutate(n = 1) %>%
  filter(patid %in% pat_summary_data$patid)

if (debugMode == TRUE) {
  unique_patids <- unique(clinicalDim$patid)
  sample_size <- min(5000, length(unique_patids))
  clinicalDim <- clinicalDim %>%
    filter(patid %in% sample(unique_patids, sample_size))
}

clinicalEvs <- read_parquet(paste0(dimension, "_for_HDPS_ever_mapped.parquet")) %>%
  mutate_at(c("patid"), as.character) %>%
  mutate(n = 1) %>%
  filter(patid %in% pat_summary_data$patid)

if (debugMode == TRUE) {
  unique_patids <- unique(clinicalEvs$patid)
  sample_size <- min(5000, length(unique_patids))
  clinicalEvs <- clinicalEvs %>%
    filter(patid %in% sample(unique_patids, sample_size))
}

# Therapy
dimension <- "drugs"
therapyDim <- read_parquet(paste0(dimension, "_for_HDPS_mapped.parquet")) %>%
  mutate_at(c("patid"), as.character) %>%
  mutate(n = 1) %>%
  filter(patid %in% pat_summary_data$patid)

if (debugMode == TRUE) {
  unique_patids <- unique(therapyDim$patid)
  sample_size <- min(5000, length(unique_patids))
  therapyDim <- therapyDim %>%
    filter(patid %in% sample(unique_patids, sample_size))
}

# Hospital
dimension <- "HES"
hospDim <- read_parquet(paste0(dimension, "_for_HDPS_mapped.parquet")) %>%
  mutate_at(c("patid"), as.character) %>%
  mutate(n = 1) %>%
  filter(patid %in% pat_summary_data$patid)

if (debugMode == TRUE) {
  unique_patids <- unique(hospDim$patid)
  sample_size <- min(10000, length(unique_patids))
  hospDim <- hospDim %>%
    filter(patid %in% sample(unique_patids, sample_size))
}

#specify variables that are not covariates
exclude_vars <- c(
  "pracid",
  "linkdate",
  "ltra",
  "asthma",
  "other_resp_disease",
  "emis_ddate",
  "cprd_ddate",
  "regstart",
  "regend",
  "death_date_ons",
  "death_date_cprd",
  "death_date",
  "covid_death_date",
  "death",
  "enddate",
  "copd_date",
  "ics_ever",
  "control_ever",
  "triple_ever",
  "eth5",
  "bmi",
  "timeout1",
  "timeout2",
  "timeout3",
  "timeout_death_any",
  "timeout_any_covid_hes",
  "pos_covid_test_present",
  "time_origin",
  "timeinstudy1",
  "timeinstudy2",
  "timeinstudy3",
  "timeinstudy_death_any",
  "timeinstudy_covid_hes_any",
  "hos_pre_death",
  "pos_covid_test_date",
  "covid_hes_date",
  "any_covid_hes_date",
  "prim_cause_covid",
  "missing_ons",
  "baseline_ics",
  "baseline_control",
  "treatgroup",
  "treat",
  "baseline_triple",
  "covid_hes_present",
  "covid_death_present",
  "any_death_present",
  "any_covid_hes_present",
  "covid_present",
  "exacerbations", 
  "obese4cat", 
  "bmi",
  "smokstatus")

# Remove the specified variables from the dataframe
hdpsCovariates <- pat_summary_data %>% select(-all_of(exclude_vars))

# Save hdpsCovariates to a Parquet file in the HDPS_folder/outputs directory
write_parquet(hdpsCovariates, 
              file.path(HDPS_folder, "outputs", paste0("covariates_no_hdps", output_ext, ".parquet")))

hdpsCohort <- pat_summary_data %>% select(
  "patid",
  "treatgroup",
  "covid_hes_present",
  "covid_death_present",
  "baseline_triple"
)


# #### Step 2: Sort codes by prevalence within each dimension -------------

# Calculate denominator for prevalence calculation
denom <- nrow(hdpsCohort)

# Define dimensions
dimNum <- c(1, 2, 3)

for (x in dimNum) {
  #  browser()
  if (x == 1) {
    dimension <- clinicalDim
    nam <- paste("prevClinical")
    nam1 <- paste("codeDistClinical")
    nam2 <- paste("codeTotsClinical")
  }
  if (x == 2) {
    dimension <- therapyDim
    nam <- paste("prevTherapy")
    nam1 <- paste("codeDistTherapy")
    nam2 <- paste("codeTotsTherapy")
  }
  if (x == 3) {
    dimension <- hospDim
    nam <- paste("prevHosp")
    nam1 <- paste("codeDistHosp")
    nam2 <- paste("codeTotsHosp")
  }
  
  # Calculate prevalence for each code
  prev <- dimension %>%
    group_by(patid, code) %>%
    slice(1) %>%
    ungroup() %>%
    group_by(code) %>%
    count() %>%
    ungroup() %>%
    mutate(prev = n / denom) %>%
    mutate(prev = if_else(prev > 0.5, 1 - prev, prev)) %>%
    mutate(rank = dense_rank(-prev)) %>%
    arrange(rank) %>%
    mutate(dim = x) %>%
    # restrict to top 500 in each dimension
    filter(rank <= 500) #!!!!!!!!!!!!!!!!!!!! Do we need?
  
  assign(nam, prev)
  
  #### Step 3: Assess recurrence of codes
  # Calculate code totals for each code
  codeTots <- dimension %>%
    filter(code %in% !!prev$code) %>% # restrict to those select
    group_by(patid, code) %>%
    count() %>%
    ungroup() %>%
    mutate(dim = x)
  
  #calculate the quantiles for each code
  codeDists <- codeTots %>%
    group_by(code) %>%
    summarise(q2 = quantile(n, c(0.5)), q3 = quantile(n, c(0.75))) %>%
    mutate(dropQ2 = ifelse(q2 == 1.00 , TRUE, FALSE)) %>%
    ungroup() %>%
    group_by(code) %>%
    mutate(dropQ3 = ifelse(q3 == q2, TRUE, FALSE)) %>%
    ungroup() %>%
    mutate(dim = x)
  
  assign(nam1, codeDists)
  assign(nam2, codeTots)
  
}

# Filter the ever dims
prev <- rbind(prevClinical, prevTherapy, prevHosp)

clinicalEvs <- clinicalEvs %>%
  filter(code %in% !!prev$code)

#generate dim = 4 for ever codes
clinicalEvs <- clinicalEvs %>%
  mutate(dim = 4)

# Incorporate ever information
codeTots <- rbind(codeTotsClinical, codeTotsTherapy, codeTotsHosp, clinicalEvs) %>%
  group_by(patid, code, dim) %>%
  arrange(patid, code, -n) %>%
  #pick the dimension with the most occurrences
  slice(1) %>%
  ungroup()

# save codeTots to parquet
output_path <- file.path(HDPS_folder, "outputs", paste0("codeTots_debug", output_ext, ".parquet"))

if (debugMode == TRUE) {
  write_parquet(codeTots, output_path)
} else {
  write_parquet(codeTots,
                file.path(HDPS_folder, "outputs", paste0("codeTots", output_ext, ".parquet")))
}

# Combine code distributions from all dimensions
codeDists <- rbind(codeDistClinical, codeDistTherapy, codeDistHosp)

# Clean up
rm(list = ls()[!ls() %in% c(
  "Datadir_copd",
  "HDPS_folder",
  "debugMode",
  "hdpsCohort",
  "output_ext",
  "exclude_triple",
  "hdpsCovariates",
  "codeTots",
  "codeDists",
  "prev",
  "period"
)])

# Initialize codes_dim data table
codes_dim <- data.table()

# Convert to data.tables for faster processing
setDT(codeDists)
setDT(codeTots)
setDT(hdpsCohort)

# Function to process each code
process_code <- function(code_row) {
  codeInfo <- code_row
  dim_suffix <- paste0("d", codeInfo$dim, "_")
  
  if (codeInfo$dim != "4") {
    codeOnce <- paste0(dim_suffix, codeInfo$code, "_once")
  } else {
    codeOnce <- paste0(dim_suffix, codeInfo$code, "_ever")
  }
  
  codeSpor <- paste0(dim_suffix, codeInfo$code, "_spor")
  codeFreq <- paste0(dim_suffix, codeInfo$code, "_freq")
  
  # Filter codeTots for the specific code
  codeTemp <- codeTots[code == codeInfo$code]
  
  # Calculate variables
  result <- codeTemp[, .(
    temp_once = as.integer(.N >= 1),
    temp_spor = if (!codeInfo$dropQ2)
      as.integer(.N >= codeInfo$q2)
    else
      NA,
    temp_freq = if (!codeInfo$dropQ3)
      as.integer(.N >= codeInfo$q3)
    else
      NA
  ), by = patid]
  
  # Rename columns
  setnames(result,
           c("temp_once", "temp_spor", "temp_freq"),
           c(codeOnce, codeSpor, codeFreq))
  
  # Remove NA columns
  result <- result[, which(unlist(lapply(result, function(x)
    ! all(is.na(
      x
    ))))), with = FALSE]
  
  # Create new rows for codes_dim
  new_rows <- data.table(code = setdiff(names(result), "patid"), dim = codeInfo$dim)
  
  list(result = result, new_rows = new_rows)
}

# Set up parallel processing
plan(multisession)
options(future.globals.maxSize = 1000 * 1024 ^ 2) # 1GB

# Process codes in parallel
results <- future_lapply(seq_len(nrow(codeDists)), function(i) {
  process_code(codeDists[i, ])
}, future.seed = TRUE)

# Remove any NULL results
results <- results[!sapply(results, is.null)]

starttime1 <- Sys.time()
# Combine results
processed_data <- Reduce(function(x, y)
  merge(x, y$result, by = "patid", all = TRUE),
  c(list(data.table(
    patid = unique(codeTots$patid)
  )), results))

# Update hdpsCohort
hdpsCohort <- merge(hdpsCohort, processed_data, by = "patid", all.x = TRUE)
# Fill NA values with 0
hdpsCohort[is.na(hdpsCohort)] <- 0

# Create codes_dim
codes_dim <- rbindlist(lapply(results, function(x)
  x$new_rows))

endtime1 <- Sys.time()
runtime1 <- endtime1 - starttime1
# Clean up
rm(results)

rm(list = setdiff(
  ls(),
  c(
    "Datadir_copd",
    "HDPS_folder",
    "debugMode",
    "output_ext",
    "exclude_triple",
    "hdpsCohort",
    "hdpsCovariates",
    "codeTots",
    "codeDists",
    "prev",
    "period",
    "codes_dim")
))

gc()

# Ranking confounders for potential for causing bias

columns_to_include <- setdiff(names(hdpsCohort), "treatgroup")

exposed <- "treatgroup"
outcome <- c("covid_hes_present", "covid_death_present")
vars <- colnames(hdpsCohort[, !names(hdpsCohort) %in% c("patid",
                                                        "indexdate",
                                                        "baseline_triple",
                                                        "dim",
                                                        exposed,
                                                        outcome)])

calculate_bias <- function(hdpsCohort, exposed, outcome, v) {
  n <- nrow(hdpsCohort)
  e1 <- sum(hdpsCohort[[exposed]] == 1)
  e0 <- sum(hdpsCohort[[exposed]] == 0)
  d1 <- sum(hdpsCohort[[outcome]] == 1)
  d0 <- sum(hdpsCohort[[outcome]] == 0)
  
  # Group by 'exposed' and summarize
  tempFrameEx <- hdpsCohort %>%
    group_by(!!sym(exposed)) %>%
    summarize(sum = sum(as.numeric(as.character(!!sym(
      v)))))
  
  c1 <- sum(tempFrameEx$sum[tempFrameEx[[exposed]] == 1])
  c0 <- n - c1
  e1c1 <- tempFrameEx$sum[tempFrameEx[[exposed]] == 1]
  e0c1 <- tempFrameEx$sum[tempFrameEx[[exposed]] == 0]
  e1c0 <- e1 - e1c1
  e0c0 <- e0 - e0c1
  
  # Group by 'outcome' and summarize
  tempFrameOut <- hdpsCohort %>%
    group_by(!!sym(outcome)) %>%
    summarize(sum = sum(as.numeric(as.character(!!sym(
      v)))))
  
  d1c1 <- tempFrameOut$sum[tempFrameOut[[outcome]] == 1]
  d0c1 <- tempFrameOut$sum[tempFrameOut[[outcome]] == 0]
  d1c0 <- d1 - d1c1
  d0c0 <- d0 - d0c1
  
  pc1 <- e1c1 / e1
  pc0 <- e0c1 / e0
  
  rrCE <- pc1 / pc0
  
  rrCE <- if (is.na(rrCE) | rrCE == 0 | is.nan(rrCE))
    NA
  else
    rrCE
  
  rrCD <- (d1c1 / c1) / (d1c0 / c0)
  rrCD <- if (is.na(rrCD) | rrCD == 0 | is.nan(rrCD))
    NA
  else
    rrCD
  
  bias <- (pc1 * (rrCD - 1) + 1) / (pc0 * (rrCD - 1) + 1)
  absLogBias <- abs(log(bias))
  ce_strength <- abs(rrCE - 1)
  cd_strength <- abs(rrCD - 1)
  
  data.frame(
    code = v,
    e1 = e1,
    e0 = e0,
    d1 = d1,
    d0 = d0,
    c1 = c1,
    c0 = c0,
    e1c1 = e1c1,
    e0c1 = e0c1,
    e1c0 = e1c0,
    e0c0 = e0c0,
    d1c1 = d1c1,
    d0c1 = d0c1,
    d1c0 = d1c0,
    d0c0 = d0c0,
    pc1 = pc1,
    pc0 = pc0,
    rrCE = rrCE,
    rrCD = rrCD,
    bias = bias,
    absLogBias = absLogBias,
    ceStrength = ce_strength,
    cdStrength = cd_strength
  )
}

run_analysis_for_outcome <- function(outcome) {
  total_codes <- nrow(codes_dim)
  
  # Create a progress bar
  pb <- progress_bar$new(
    format = "[:bar] :percent Elapsed: :elapsed ETA: :eta",
    total = total_codes,
    clear = FALSE,
    width = 100
  )
  
  # biasInfo <- lapply(seq_len(total_codes), function(i) {
  #   v <- codes_dim$code[i]
  #   result <- calculate_bias(hdpsCohort, exposed, outcome, v)
  #   pb$tick()  # Update the progress bar
  #   return(result)
  # })
  
  # biasInfo <- do.call(rbind, biasInfo)
  # biasInfo <- biasInfo %>% 
  #   arrange(-absLogBias) %>% 
  #   mutate(rank = row_number())
  # # Return the biasInfo dataset
  # return(biasInfo)
  
  # Return the biasInfo dataset
  purrr::map_dfr(codes_dim$code, function(v) {
    result <- calculate_bias(hdpsCohort, exposed, outcome, v)
    pb$tick()  # Update the progress bar
    return(result)
  }) %>% 
    arrange(-absLogBias) %>% 
    mutate(rank = row_number())
}

# Run the analysis for both outcomes
results_outcome1 <- run_analysis_for_outcome("covid_hes_present")

#save the results
saveRDS(results_outcome1, "results_outcome1.rds")
#read in results
results_outcome1 <- readRDS("results_outcome1.rds")

results_outcome2 <- run_analysis_for_outcome("covid_death_present")
saveRDS(results_outcome2, "results_outcome2.rds")
results_outcome2 <- readRDS("results_outcome2.rds")

create_and_save_hdps_results <- function(results_outcome, HDPS_folder, cohort, output_ext, debugMode = FALSE, outcome) {

  # Function to create top N lists
  create_top_k <- function(k) {
    top_k <- results_outcome %>% filter(rank <= k)
    selectedVars <- c("patid", names(hdpsCovariates)[-1], top_k$code)
    cohort_k <- cohort %>% 
      dplyr::select(patid, all_of(top_k$code)) %>%
      left_join(hdpsCovariates, by = "patid") %>%
      dplyr::select(all_of(selectedVars))
    list(top_k = top_k, cohort = cohort_k)
  }
  
  top_codes <- c(10, 100, 250, 500, 750, 1000)
  top_codes <- c(10)
  
  # Create top lists for different N values
  top_lists <- lapply(top_codes, create_top_k)
  names(top_lists) <- paste0("top", top_codes)
  
  # Save cohorts and top lists for different N values
  for (k in top_codes) {
    
    file_name_cohort <- paste0(ifelse(debugMode, "debug_", ""),
                               "covariates_HDPS_",
                               k,
                               outcome,
                               output_ext,
                               ".parquet")
    write_parquet(top_lists[[paste0("top", k)]]$cohort,
              file.path(HDPS_folder, "outputs", file_name_cohort))
    
    file_name_top <- paste0(ifelse(debugMode, "debug_", ""),
                            "top_",
                            k,
                            "_HDPS_",
                            outcome,
                            output_ext,
                            ".csv")
    write_csv(top_lists[[paste0("top", k)]]$top_k,
              file.path(HDPS_folder, "outputs", file_name_top))
  }
  
  # Save full biasInfo
  file_name_biasInfo <- paste0(ifelse(debugMode, "debug_", ""),
                               "HDPS_biasInfo_",
                               outcome,
                               output_ext,
                               ".csv")
  write_csv(results_outcome,
            file.path(HDPS_folder, "outputs", file_name_biasInfo))
  
  # Return results directly
  results_outcome
  
  hdpsCohort
}


# For the first outcome
df_outcome1 <- create_and_save_hdps_results(
  results_outcome = results_outcome1, 
  HDPS_folder = HDPS_folder, 
  cohort = hdpsCohort,
  debugMode = debugMode,
  output_ext = output_ext,
  outcome = "covid_hes_present"
)

# For the second outcome
df_outcome2 <- create_and_save_hdps_results(
  results_outcome = results_outcome2, 
  HDPS_folder = HDPS_folder, 
  cohort = hdpsCohort,
  debugMode = debugMode,
  output_ext = output_ext,
  outcome = "covid_death_present"
)

