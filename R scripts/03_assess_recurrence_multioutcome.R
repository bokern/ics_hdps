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
# Flag to enable debugging; if TRUE, script will only be run on 1,000 patients selected at random
debugMode = F

sample_size <- 1000
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
  sample_size <- min(1000, length(unique_patids))
  clinicalDim <- clinicalDim %>%
    filter(patid %in% sample(unique_patids, sample_size))
}

clinicalEvs <- read_parquet(paste0(dimension, "_for_HDPS_ever_mapped.parquet")) %>%
  mutate_at(c("patid"), as.character) %>%
  mutate(n = 1) %>%
  filter(patid %in% pat_summary_data$patid)

if (debugMode == TRUE) {
  unique_patids <- unique(clinicalEvs$patid)
  sample_size <- min(1000, length(unique_patids))
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
  sample_size <- min(1000, length(unique_patids))
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
  sample_size <- min(1000, length(unique_patids))
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
hdpsCovariates <- pat_summary_data %>% dplyr::select(-all_of(exclude_vars))

# Save hdpsCovariates to a Parquet file in the HDPS_folder/outputs directory
write_parquet(hdpsCovariates, 
              file.path(HDPS_folder, "outputs", paste0("covariates_no_hdps", output_ext, ".parquet")))

hdpsCohort <- pat_summary_data %>% dplyr::select(
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
    filter(rank <= 500)
  
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

clinicalEvs <- clinicalEvs %>%
  filter(code %in% !!prevClinical$code) %>%
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
  "period",
  "HDPS_github"
)])

# Convert to data.tables for faster processing
setDT(codeDists)
setDT(codeTots)
setDT(hdpsCohort)

process_code <- function(code_row) {
  
  codeInfo <- code_row
  dim_suffix <- paste0("d", codeInfo$dim, "_")
  
  codeOnce <- paste0(dim_suffix, codeInfo$code, "_once")
  codeSpor <- paste0(dim_suffix, codeInfo$code, "_spor")
  codeFreq <- paste0(dim_suffix, codeInfo$code, "_freq")
  
  # Filter codeTots for the specific code
  codeTemp <- codeTots[code == codeInfo$code & dim == codeInfo$dim,]
  
  result <- codeTemp[, .(
    temp_once = as.integer(n >= 1),
    temp_spor = if (!codeInfo$dropQ2) as.integer(n >= codeInfo$q2) else NA_integer_,
    temp_freq = if (!codeInfo$dropQ3) as.integer(n >= codeInfo$q3) else NA_integer_
  ), by = patid]
  
  # Remove temp columns that contain only NAs
  result <- result %>%
    dplyr::select(where(~!all(is.na(.))))
  
  # Create new_names vector, excluding empty entries
  new_names <- c(codeOnce, 
                 if (!codeInfo$dropQ2 && !all(is.na(result$temp_spor))) codeSpor else character(0), 
                 if (!codeInfo$dropQ3 && !all(is.na(result$temp_freq))) codeFreq else character(0))
  new_names <- new_names[new_names != ""]
  
  # Rename columns, matching the number of columns in result (excluding patid)
  setnames(result, c("patid", new_names))
  
  # If this is a dimension 4 code, update or add d1_xxx_once
  if (codeInfo$dim == "4") {
    d1_code_once <- paste0("d1_", codeInfo$code, "_once")
    if (d1_code_once %in% names(result)) {
      # If d1_xxx_once already exists, update it
      result[, (d1_code_once) := pmax(get(d1_code_once), get(codeOnce))]
    } else {
      # If d1_xxx_once doesn't exist, add it
      result[, (d1_code_once) := get(codeOnce)]
      new_names <- c(d1_code_once, new_names)
    }
  }
  
  return(result)
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

processed_data <- rbindlist(results, use.names = TRUE, fill = TRUE)
# Ensure one row per patient by taking the maximum value for each column
processed_data <- processed_data[, lapply(.SD, max, na.rm = TRUE), by = patid]

#export the results to parquet
write_parquet(processed_data, file.path(HDPS_folder, "outputs", paste0("cohort_hdps_covariates", output_ext, ".parquet")))
processed_data <- read_parquet(file.path(HDPS_folder, "outputs", paste0("cohort_hdps_covariates", output_ext, ".parquet")))

#check if there are any columns that are all NA in processed_data
if (any(colSums(is.na(processed_data)) == nrow(processed_data))) {
  stop("There are columns that are all NA in processed_data")
}

# Update hdpsCohort
hdpsCohort <- merge(hdpsCohort, processed_data, by = "patid", all.x = TRUE)

# Fill NA values with 0
hdpsCohort[is.na(hdpsCohort)] <- 0

# Clean up
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
    "HDPS_github")
))

gc()

# Ranking confounders for potential for causing bias
columns_to_include <- setdiff(names(hdpsCohort), "treatgroup")

calculate_bias <- function(hdpsCohort, v) {
  n <- nrow(hdpsCohort)
  e1 <- sum(hdpsCohort[[exposed]] == 1, na.rm = TRUE)
  e0 <- sum(hdpsCohort[[exposed]] == 0, na.rm = TRUE)
  d1 <- sum(hdpsCohort[[outcome]] == 1, na.rm = TRUE)
  d0 <- sum(hdpsCohort[[outcome]] == 0, na.rm = TRUE)

  c1 <- sum(hdpsCohort[[v]] == 1, na.rm = TRUE)
  c0 <- n - c1 # number of people without the covariate
  
  # Use the full variable name (including dimension and frequency) instead of just the code
  tempFrameEx <- hdpsCohort %>%
    group_by(!!sym(exposed)) %>%
    summarize(sum = sum(as.numeric(as.character(!!sym(v))), na.rm = TRUE), .groups = 'drop')
  
  e1c1 <- ifelse(1 %in% tempFrameEx[[exposed]], tempFrameEx$sum[tempFrameEx[[exposed]] == 1], NA) # number of people with the covariate and exposed
  e0c1 <- ifelse(0 %in% tempFrameEx[[exposed]], tempFrameEx$sum[tempFrameEx[[exposed]] == 0], NA) # number of people with the covariate and not exposed
  e1c0 <- e1 - e1c1
  e0c0 <- e0 - e0c1
  
  tempFrameOut <- hdpsCohort %>%
    group_by(!!sym(outcome)) %>%
    summarize(sum = sum(as.numeric(as.character(!!sym(v))), na.rm = TRUE), .groups = 'drop')
  
  d1c1 <- ifelse(1 %in% tempFrameOut[[outcome]], tempFrameOut$sum[tempFrameOut[[outcome]] == 1], NA)
  d0c1 <- ifelse(0 %in% tempFrameOut[[outcome]], tempFrameOut$sum[tempFrameOut[[outcome]] == 0], NA)
  d1c0 <- d1 - d1c1
  d0c0 <- d0 - d0c1
  
  #if any of the dc or ec cells are 0, add 0.1 to all cells
  if (d1c1 == 0 | d1c0 == 0 | d0c1 == 0 | d0c0 == 0| e1c1 == 0 | e1c0 == 0 | e0c1 == 0 | e0c0 == 0) {
    d1c1 <- d1c1 + 0.1
    d1c0 <- d1c0 + 0.1
    d0c1 <- d0c1 + 0.1
    d0c0 <- d0c0 + 0.1
    e1c1 <- e1c1 + 0.1
    e1c0 <- e1c0 + 0.1
    e0c1 <- e0c1 + 0.1
    e0c0 <- e0c0 + 0.1
  }
  
  pc1 <- e1c1 / e1
  pc0 <- e0c1 / e0
  
  rrCE <- pc1 / pc0
  rrCD <- (d1c1 / c1) / (d1c0 / c0)

  bias <- (pc1 * (rrCD - 1) + 1) / (pc0 * (rrCD - 1) + 1)
  absLogBias <- abs(log(bias))
  ce_strength <- abs(rrCE - 1)
  cd_strength <- abs(rrCD - 1)
  
  # Calculate log risk ratios for IV checks
  logRRCE <- log(rrCE)
  logRRCD <- log(rrCD)
  
  # Add IV flags based on the specified criteria
  iv_flag_1 <- abs(logRRCE) > 1.5 & abs(logRRCD) < 0.5
  iv_flag_2 <- abs(logRRCE) > 1.1 & abs(logRRCD) < 0.5
  
  return(list(
    variable = v,  # Use the full variable name instead of just the code
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
    cdStrength = cd_strength,
    logRRCE = logRRCE,
    logRRCD = logRRCD,
    iv_flag_1 = iv_flag_1,
    iv_flag_2 = iv_flag_2
  ))
}

# Main workflow
exposed <- "treatgroup"
outcomes <- c("covid_hes_present", "covid_death_present")
covariates <- setdiff(names(hdpsCohort), c("patid", "baseline_triple", exposed, outcomes))

# Initialize a list to store results for each outcome
results_list <- list()

# Analyze for both outcomes
for (outcome in outcomes) {
  cat("Processing outcome:", outcome, "\n")
  
  # Create progress bar for covariates
  pb <- progress_bar$new(
    format = "[:bar] :percent Covariate: :current/:total Elapsed: :elapsed ETA: :eta",
    total = length(covariates),
    clear = FALSE,
    width = 100
  )
  
  results <- map_dfr(covariates, function(v) {
    bias_results <- calculate_bias(hdpsCohort, v)
    pb$tick()  # Update progress bar
    # Convert the list to a data frame
    as.data.frame(bias_results, stringsAsFactors = FALSE) %>%
      tibble::as_tibble()
  })
  
  results <- results %>%
    arrange(desc(absLogBias)) %>%
    mutate(rank = row_number())
  
  # Save results to HDPS folder/output
  file_path <- file.path(HDPS_folder, "/outputs", paste0("HDPS_biasInfo_", outcome, output_ext, ".csv"))
  write_csv(results, file_path)

  # Create and save top K lists, not including codes where the iv_flag_2 is TRUE
  top_codes <- c(100, 250, 500, 750, 1000)
  top_lists <- map(top_codes, function(k) {
    top_k <- results %>%
      arrange(rank) %>%
      filter(iv_flag_2 == FALSE) %>%
      slice_head(n = k) %>%
      dplyr::select(variable, rank)
    
    cohort <- hdpsCohort %>%
      dplyr::select(c("patid", "baseline_triple", exposed, outcome, top_k$variable))
    
    list(cohort = cohort, top_k = top_k)
  })
  
  names(top_lists) <- paste0("top", top_codes)
  
  walk2(top_lists, top_codes, function(x, y) {
    write_parquet(x$cohort, file.path(HDPS_folder, "outputs", paste0("covariates_HDPS_", y, "_", outcome, output_ext, ".parquet")))
    write_csv(x$top_k, file.path(HDPS_folder, "outputs", paste0("top_", y, "_HDPS_", outcome, output_ext, ".csv")))
  })
  
  # Store results for this outcome
  results_list[[outcome]] <- list(results = results, top_lists = top_lists)
}

# Access and save results for all outcomes
for (outcome in outcomes) {
  results_outcome <- results_list[[outcome]]$results
  write_csv(results_outcome, paste0("results_", outcome, output_ext, ".csv"))
}