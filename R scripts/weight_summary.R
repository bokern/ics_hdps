# -----------------------------------------------------------------------------
# PROGRAM NAME:  summarize_iptw_weights.R
# PROJECT:       COPD Study Weight Analysis
# AUTHOR:        
# DATE CREATED:  April 23, 2025
# NOTES:         Summarizes IPTW weights (predefined and HDPS) before and after trimming
# -----------------------------------------------------------------------------

# Load libraries -------------------------------------------------------------
packages <- c(
  "tidyverse",
  "arrow",
  "readr"
)

# Install missing packages
installed_packages <- packages %in% rownames(installed.packages())
if (any(installed_packages == FALSE)) {
  install.packages(packages[!installed_packages])
}

invisible(lapply(packages, function(pkg) {
  suppressMessages(library(pkg, character.only = TRUE, verbose = FALSE))
}))

# Set your HDPS folder path here
setwd(Datadir_copd)

# Parameters for analyses
outcomes <- c("covid_hes_present", "covid_death_present")
cohort_exts <- c("", "_no_triple")
topVars <- c("100", "250", "500", "750", "1000")

# Unified function to get weight summaries (both pre and post-trimming)
get_weight_summary <- function(cohort_ext, outcome, k, is_trimmed = FALSE) {
  # Determine file path and weight column based on parameters
  if (k == "predefined") {
    # For predefined weights
    if (!is_trimmed) {
      file_path <- file.path(
        HDPS_folder,
        "outputs",
        paste0("HDPS_predefined_", outcome, cohort_ext, ".parquet")
      )
    } else {
      file_path <- file.path(
        HDPS_folder,
        "outputs",
        paste0("HDPS_predefined_", outcome, cohort_ext, "_trimmed.parquet")
      )
    }
    weight_col <- "iptw_weight_predefined"
    weight_type <- "predefined"
    k_value <- "0" # Using 0 to represent predefined in results
  } else {
    # For HDPS weights
    if (!is_trimmed) {
      file_path <- file.path(
        HDPS_folder,
        "outputs",
        paste0("covariates_HDPS_", k, "_", outcome, cohort_ext, ".parquet")
      )
    } else {
      file_path <- file.path(
        HDPS_folder,
        "outputs",
        paste0("HDPS_", k, "_", outcome, cohort_ext, "_trimmed.parquet")
      )
    }
    weight_col <- "iptw_weight"
    weight_type <- "hdps"
    k_value <- k
  }
  
  stage_suffix <- ifelse(is_trimmed, "after", "before")
  
  if (file.exists(file_path)) {
    # Load data
    data <- read_parquet(file_path)
    
    # Get summary of weights by treatment group
    weight_summary <- data %>%
      group_by(treatgroup) %>%
      summarise(
        mean_weight = mean(!!sym(weight_col), na.rm = TRUE),
        median_weight = median(!!sym(weight_col), na.rm = TRUE),
        n = n()
      ) %>%
      pivot_wider(
        names_from = treatgroup,
        values_from = c(mean_weight, median_weight, n),
        names_prefix = "group_"
      ) %>%
      mutate(
        weight_type = weight_type,
        cohort = ifelse(cohort_ext == "", "normal", "no_triple"),
        outcome = outcome,
        k = k_value,
        stage = stage_suffix
      )
    
    return(weight_summary)
  } else {
    # Warning for missing files
    warning(paste("File not found:", file_path))
    return(NULL)
  }
}

# Function to calculate excluded counts from trimming results files
get_exclusion_counts <- function(cohort_ext, outcome, k) {
  if (k == "predefined") {
    file_path <- file.path(
      HDPS_folder,
      "outputs",
      paste0("HDPS_", outcome, "_predefined", cohort_ext, "_trimming_results.csv")
    )
    k_value <- "0"
  } else {
    file_path <- file.path(
      HDPS_folder,
      "outputs",
      paste0("HDPS_", outcome, "_", k, cohort_ext, "_trimming_results.csv")
    )
    k_value <- k
  }
  
  if (file.exists(file_path)) {
    trimming_data <- read_csv(file_path, show_col_types = FALSE)
    
    # Extract exclusion counts
    exclusion_summary <- trimming_data %>%
      mutate(
        weight_type = ifelse(grepl("predefined", weights), "predefined", "hdps"),
        cohort = ifelse(cohort == "", "normal", "no_triple"),
        n_excluded_ics = n_ics,
        n_excluded_labalama = n_labalama,
        k = k_value
      ) %>%
      select(weight_type, cohort, outcome, k, n_excluded_ics, n_excluded_labalama)
    
    return(exclusion_summary)
  } else {
    warning(paste("Trimming results file not found:", file_path))
    return(NULL)
  }
}

# Create parameter combinations including predefined
param_combinations <- bind_rows(
  # For HDPS analyses
  expand_grid(
    cohort_ext = cohort_exts,
    outcome = outcomes,
    k = topVars
  ),
  # For predefined analyses
  expand_grid(
    cohort_ext = cohort_exts,
    outcome = outcomes,
    k = "predefined"
  )
)

# Collect weight summaries for before trimming
weight_summaries_pre <- map_df(1:nrow(param_combinations), function(i) {
  get_weight_summary(
    param_combinations$cohort_ext[i],
    param_combinations$outcome[i],
    param_combinations$k[i],
    is_trimmed = FALSE
  )
})

# Collect weight summaries for after trimming
weight_summaries_post <- map_df(1:nrow(param_combinations), function(i) {
  get_weight_summary(
    param_combinations$cohort_ext[i],
    param_combinations$outcome[i],
    param_combinations$k[i],
    is_trimmed = TRUE
  )
})

# Rename columns for consistent joining
weight_summaries_pre <- weight_summaries_pre %>%
  rename_with(~paste0(., "_before"), matches("^group_"))

weight_summaries_post <- weight_summaries_post %>%
  rename_with(~paste0(., "_after"), matches("^group_"))

# Get exclusion counts
exclusion_counts <- map_df(1:nrow(param_combinations), function(i) {
  get_exclusion_counts(
    param_combinations$cohort_ext[i],
    param_combinations$outcome[i],
    param_combinations$k[i]
  )
})

# Combine all summaries
weight_summary_combined <- weight_summaries_pre %>%
  full_join(weight_summaries_post, 
            by = c("weight_type", "cohort", "outcome", "k")) %>%
  full_join(exclusion_counts, 
            by = c("weight_type", "cohort", "outcome", "k")) %>%
  mutate(
    # Calculate excluded counts if not available from trimming results
    n_excluded_ics = coalesce(n_excluded_ics, as.integer(group_1_n_before - group_1_n_after)),
    n_excluded_labalama = coalesce(n_excluded_labalama, as.integer(group_0_n_before - group_0_n_after))
  ) %>%
  # Organize columns as requested
  select(
    cohort,
    outcome,
    k,
    weight_type,
    mean_weight_ics_before = group_1_mean_weight_before,
    median_weight_ics_before = group_1_median_weight_before,
    mean_weight_ics_after = group_1_mean_weight_after,
    median_weight_ics_after = group_1_median_weight_after,
    mean_weight_labalama_before = group_0_mean_weight_before,
    median_weight_labalama_before = group_0_median_weight_before,
    mean_weight_labalama_after = group_0_mean_weight_after,
    median_weight_labalama_after = group_0_median_weight_after,
    n_excluded_ics,
    n_excluded_labalama
  )

# Calculate additional statistics
weight_summary_final <- weight_summary_combined %>%
  mutate(
    total_excluded = n_excluded_ics + n_excluded_labalama,
    
    # Calculate percentage excluded, handling NA and zero denominators
    pct_excluded_ics = case_when(
      is.na(n_excluded_ics) | is.na(group_1_n_after) ~ NA_real_,
      (n_excluded_ics + group_1_n_after) == 0 ~ 0,
      TRUE ~ (n_excluded_ics / (n_excluded_ics + group_1_n_after) * 100)
    ),
    
    pct_excluded_labalama = case_when(
      is.na(n_excluded_labalama) | is.na(group_0_n_after) ~ NA_real_,
      (n_excluded_labalama + group_0_n_after) == 0 ~ 0,
      TRUE ~ (n_excluded_labalama / (n_excluded_labalama + group_0_n_after) * 100)
    ),
    
    pct_excluded_total = case_when(
      is.na(total_excluded) | is.na(group_0_n_after) | is.na(group_1_n_after) ~ NA_real_,
      (total_excluded + group_0_n_after + group_1_n_after) == 0 ~ 0,
      TRUE ~ (total_excluded / (total_excluded + group_0_n_after + group_1_n_after) * 100)
    )
  ) %>%
  select(
    cohort,
    outcome,
    k,
    weight_type,
    mean_weight_ics_before,
    median_weight_ics_before,
    mean_weight_ics_after,
    median_weight_ics_after,
    mean_weight_labalama_before,
    median_weight_labalama_before,
    mean_weight_labalama_after,
    median_weight_labalama_after,
    n_excluded_ics,
    n_excluded_labalama,
    total_excluded,
    pct_excluded_ics,
    pct_excluded_labalama,
    pct_excluded_total
  )

# Write output
write_csv(
  weight_summary_final,
  file.path(HDPS_folder, "outputs", "iptw_weight_summary.csv")
)

# Create a more readable summary table for presentation
readable_summary <- weight_summary_final %>%
  mutate(
    analysis_type = case_when(
      k == "0" ~ "Predefined covariates",
      TRUE ~ paste0("HDPS: Top ", k, " covariates")
    ),
    cohort_type = case_when(
      cohort == "normal" ~ "Full cohort",
      cohort == "no_triple" ~ "No triple therapy"
    ),
    outcome_type = case_when(
      outcome == "covid_hes_present" ~ "COVID hospitalization",
      outcome == "covid_death_present" ~ "COVID mortality"
    )
  ) %>%
  select(
    analysis_type,
    cohort_type, 
    outcome_type,
    mean_weight_ics_before,
    mean_weight_ics_after,
    mean_weight_labalama_before,
    mean_weight_labalama_after,
    n_excluded_ics,
    n_excluded_labalama,
    pct_excluded_total
  ) %>%
  arrange(outcome_type, cohort_type, analysis_type)

write_csv(
  readable_summary,
  file.path(HDPS_folder, "outputs", "iptw_weight_summary_readable.csv")
)

cat("Weight summary analysis complete. Results saved to:", 
    file.path(HDPS_folder, "outputs", "iptw_weight_summary.csv"), "\n")