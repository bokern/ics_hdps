# ==============================================================================
# PROPENSITY SCORE ANALYSIS WITH MULTIPLE IMPUTATION FOR MISSING ETHNICITY
# ==============================================================================
# Author: Marleen Bokern
# Date: 03/2023
# Purpose: Generate propensity scores and perform IPTW analysis on COPD data
#          with multiple imputation for missing ethnicity data

# Load required packages ------------------------------------------------------
packages <- c("tidyverse", "arrow", "mice", "MatchThem", "parallelly", 
              "furrr", "survey", "cobalt", "ggplot2", "survival")
installed_packages <- packages %in% rownames(installed.packages())
if (any(installed_packages == FALSE)) {
  install.packages(packages[!installed_packages])
}
lapply(packages, library, character.only = TRUE)

# Set working directory and load data ----------------------------------------
setwd(Datadir_copd)
df <- read_parquet("copd_wave1_60d.parquet")

start_time <- Sys.time()

# DATA PREPARATION ============================================================

# Filter to treatment groups and remove baseline triple therapy patients

cohort_exts <- c("", "_no_triple")

for (cohort_ext in cohort_exts) {
  
  if (cohort_ext == "_no_triple") {
    subset_df <- df[!is.na(df$treatgroup) & df$baseline_triple != 1, ]
  }
  else {
    subset_df <- df[!is.na(df$treatgroup), ]
  }
  
  subset_df$treat <- factor(subset_df$treat)
  subset_df$treat <- relevel(subset_df$treat, ref = "LABA/LAMA")
  
  # Select relevant variables for analysis
  vars_to_keep <- c("patid", "age_index", "gender", "imd", "bmicat", "eth", 
                    "diabetes_present", "hypertension_present", "cvd_present", 
                    "past_asthma_present", "allcancers_present", "kidney_present", 
                    "immunosuppression_present", "flu_vacc_present", "pneumo_vacc_present", 
                    "exacerb_present", "smok", "timeinstudy1", "timeinstudy2", 
                    "timeinstudy3", "pos_covid_test_present", "covid_hes_present", 
                    "covid_death_present", "treat", "treatgroup")
  
  subset_df <- subset_df[, vars_to_keep]
  
  # Prepare ethnicity variable for imputation ----------------------------------
  # Set ethnicity "unknown" (5) as missing for imputation
  subset_df$eth <- ifelse(subset_df$eth == 5, NA, subset_df$eth)
  # Convert to factor with meaningful labels
  subset_df$eth <- factor(subset_df$eth, 
                          levels = c(1, 2, 3, 4), 
                          labels = c("White", "Asian", "Black", "Mixed"))
  
  cat("Data preparation complete. Sample size:", nrow(subset_df), "\n")
  
  # MISSING DATA ASSESSMENT ====================================================
  
  # Examine missing data patterns
  cat("\nMissing data patterns:\n")
  md.pattern(subset_df, rotate.names = TRUE)
  
  # Count missing values for model variables
  model_vars <- c("age_index", "gender", "imd", "bmicat", "eth", "diabetes_present", 
                  "hypertension_present", "cvd_present", "allcancers_present", 
                  "past_asthma_present", "kidney_present", "immunosuppression_present", 
                  "flu_vacc_present", "pneumo_vacc_present", "exacerb_present", "smok")
  
  missing_counts <- colSums(is.na(subset_df[, model_vars]))
  cat("\nMissing counts by variable:\n")
  print(missing_counts)
  
  # MULTIPLE IMPUTATION ========================================================
  
  cat("\nStarting multiple imputation...\n")
  imp_start_time <- Sys.time()
  
  # Set up imputation method (only ethnicity needs imputation via polytomous regression)
  method_vector <- make.method(subset_df)
  method_vector["eth"] <- "polyreg"  # Polytomous regression for categorical ethnicity
  
  # Perform multiple imputation using parallel processing
  set.seed(123)  # For reproducibility
  imp_mids <- futuremice(data = subset_df, 
                         method = method_vector, 
                         m = 10,           # 10 imputed datasets
                         maxit = 10,       # 10 iterations
                         parallelseed = 123, 
                         print = FALSE)
  
  #imp_mids is the mids object containing the imputed datasets
  
  imp_end_time <- Sys.time()
  cat("Multiple imputation completed in:", format(imp_end_time - imp_start_time), "\n")
  
  # Check imputation results
  print(imp_mids)
  cat("\nEthnicity imputation check:\n")
  print(imp_mids$imp$eth)
  
  # Convert to long format for analysis
  data <- complete(imp_mids, action = "long") #data is then the full dataset with all imputations (m=10) appended
  cat("Imputed datasets created. Total observations:", nrow(data), "\n")
  
  # PROPENSITY SCORE ESTIMATION ================================================
  
  cat("\nEstimating propensity scores and calculating weights...\n")
  
  # Define the PS model formula
  ps_formula <- treatgroup ~ age_index + gender + eth + bmicat + diabetes_present +
    hypertension_present + cvd_present + allcancers_present +
    past_asthma_present + kidney_present + immunosuppression_present +
    flu_vacc_present + pneumo_vacc_present + exacerb_present + smok
  
  # Use weightthem() to perform weighting across all imputations
  set.seed(123)
  weighted_mids <- weightthem(ps_formula,
                              datasets = imp_mids,
                              approach = 'within', # Estimate PS model within each imputation
                              method = 'ps',      # Use propensity score weighting
                              stabilize = TRUE)   # Use stabilized weights (your original approach)
  
  cat("Weighting complete.\n")
  
  # You can now check balance across all imputed datasets
  bal.tab(weighted_mids)
  
  # OUTCOME ANALYSIS ============================================================
  
  cat("\nStarting outcome analysis with pooling...\n")
  
  # Define outcomes
  outcomes_config <- list(
    list(outcome = "covid_hes_present", time_var = "timeinstudy2", name = "COVID Hospitalisation"),
    list(outcome = "covid_death_present", time_var = "timeinstudy3", name = "COVID Death")
  )
  
  results <- list()
  for (config in outcomes_config) {
    cat("Analyzing", config$name, "...\n")
    
    # Build the formula for the outcome model
    outcome_formula <- as.formula(paste0("Surv(", config$time_var, ", ", config$outcome, ") ~ treat"))
    
    # Run the weighted Cox model on the weighted imputed datasets
    model_fits <- with(weighted_mids,
                       svycoxph(outcome_formula, design = svydesign(ids = ~1, weights = ~weights, data = .data)))
    
    # Pool the results
    pooled_summary <- summary(pool(model_fits), exponentiate = TRUE)
    
    results[[config$name]] <- list(iptw = pooled_summary)
  }
  
  print(results)
  
  
  # EXTRACT RESULTS=================================
  
  # CREATE POOLED RESULTS DATAFRAME
  if (!exists("pooled_df")) {
  pooled_df <- data.frame()
  }
  
  for (outcome_name in names(results)) {
    iptw_results <- results[[outcome_name]]$iptw
    
    if (!is.null(iptw_results) && nrow(iptw_results) > 0) {
      # Since you used exponentiate = TRUE, the estimate is already the HR
      hr <- iptw_results$estimate[1]
      std_error <- iptw_results$std.error[1]
      p_value <- iptw_results$p.value[1]
      
      # Calculate 95% CI manually
      # Since estimate is already exponentiated, we need to work on log scale for CI
      log_hr <- log(hr)
      log_lower <- log_hr - 1.96 * std_error
      log_upper <- log_hr + 1.96 * std_error
      
      # Convert back to HR scale
      lower_ci <- exp(log_lower)
      upper_ci <- exp(log_upper)
      
      pooled_df <- rbind(pooled_df, data.frame(
        outcome = outcome_name,
        cohort = cohort_ext,
        hr = hr,
        lower_ci = lower_ci,
        upper_ci = upper_ci,
        p_value = p_value
      ))
    }
  }
}

# Save results to CSV
write.csv(pooled_df, file = file.path(HDPS_folder, "outputs", "Results", "copd_mi_iptw_results.csv"), row.names = FALSE)
