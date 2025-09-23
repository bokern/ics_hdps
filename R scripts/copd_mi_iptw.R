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
subset_df <- df[!is.na(df$treatgroup) & df$baseline_triple != 1, ]
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

imp_end_time <- Sys.time()
cat("Multiple imputation completed in:", format(imp_end_time - imp_start_time), "\n")

# Check imputation results
print(imp_mids)
cat("\nEthnicity imputation check:\n")
print(imp_mids$imp$eth)

# Convert to long format for analysis
data <- complete(imp_mids, action = "long")
cat("Imputed datasets created. Total observations:", nrow(data), "\n")

# PROPENSITY SCORE ESTIMATION ================================================

cat("\nEstimating propensity scores...\n")

# Initialize columns for propensity scores and weights
data$ps <- NA_real_
data$iptw <- NA_real_

# Define propensity score model formula
ps_formula <- treatgroup ~ age_index + gender + eth + bmicat + diabetes_present + 
  hypertension_present + cvd_present + allcancers_present + 
  past_asthma_present + kidney_present + immunosuppression_present + 
  flu_vacc_present + pneumo_vacc_present + exacerb_present + smok

# Estimate propensity scores for each imputed dataset
for (imp_num in 1:10) {
  cat("Processing imputation", imp_num, "of 10\n")
  
  # Filter data for current imputation
  imp_data <- data[data$.imp == imp_num, ]
  
  # Fit propensity score model
  ps_model <- glm(ps_formula, data = imp_data, family = "binomial")
  
  # Generate propensity scores
  imp_data$new_ps <- predict(ps_model, type = "response")
  
  # Create propensity score overlap plot
  overlap_plot <- imp_data %>% 
    ggplot(aes(x = new_ps, linetype = factor(treat))) +
    geom_density(alpha = 0.7) +
    labs(x = 'Probability of receiving treatment',
         y = 'Density',
         title = paste('Propensity Score Overlap - Imputation', imp_num)) +
    scale_x_continuous(breaks = seq(0, 1, 0.1)) +
    theme_minimal() +
    theme(legend.title = element_blank(),
          legend.position = "bottom")
  
  # Save overlap plot
  ggsave(filename = paste0(Graphdir, "ps_overlap_imp", imp_num, ".png"),
         plot = overlap_plot, width = 8, height = 5)
  
  # Update main dataset with propensity scores
  data <- data %>%
    left_join(imp_data %>% select(patid, .imp, new_ps), by = c("patid", ".imp")) %>%
    mutate(ps = ifelse(is.na(ps) & !is.na(new_ps), new_ps, ps)) %>%
    select(-new_ps)
  
  # Calculate propensity score trimming bounds (common support)
  ps_bounds <- imp_data %>% 
    group_by(treat) %>% 
    summarise(min_ps = min(new_ps), max_ps = max(new_ps), .groups = 'drop') %>% 
    summarise(lower_bound = max(min_ps), upper_bound = min(max_ps))
  
  # Apply trimming to ensure common support
  imp_data_trimmed <- imp_data %>%
    filter(new_ps >= ps_bounds$lower_bound & new_ps <= ps_bounds$upper_bound)
  
  # Calculate the overall probability of receiving treatment
  p_treat <- mean(imp_data$treatgroup == 1)
  
  # Calculate stabilized weights instead of unstabilized
  imp_data_trimmed$weight <- ifelse(imp_data_trimmed$treatgroup == 1,
                                    p_treat / imp_data_trimmed$new_ps,
                                    (1 - p_treat) / (1 - imp_data_trimmed$new_ps))
  
  # Update main dataset with weights
  data <- data %>%
    left_join(imp_data_trimmed %>% select(patid, .imp, weight), by = c("patid", ".imp")) %>%
    mutate(iptw = ifelse(is.na(iptw) & !is.na(weight), weight, iptw)) %>%
    select(-weight)
  
  cat("Imputation", imp_num, "complete. Trimmed sample size:", nrow(imp_data_trimmed), "\n")
}

# OUTCOME ANALYSIS ============================================================

# OUTCOME ANALYSIS ============================================================

cat("\nStarting outcome analysis...\n")

# Initialize lists to store model fits from each imputation
# This is necessary for pooling with Rubin's rules
fit_hes_unadj <- list()
fit_hes_iptw  <- list()
fit_death_unadj <- list()
fit_death_iptw  <- list()

# Analyze each imputed dataset
for (imp_num in 1:10) {
  cat("Analyzing imputation", imp_num, "of 10\n")
  
  # Get data for current imputation
  imp_data <- data[data$.imp == imp_num, ]
  imp_data_trimmed <- imp_data[!is.na(imp_data$iptw), ]
  
  # Set up survey designs
  svy_unadj <- svydesign(ids = ~1, data = imp_data)
  svy_iptw <- svydesign(ids = ~1, weights = ~iptw, data = imp_data_trimmed)
  
  # --- Model 1: COVID Hospitalization (Unadjusted) ---
  tryCatch({
    fit_hes_unadj[[imp_num]] <- svycoxph(
      Surv(timeinstudy2, covid_hes_present) ~ treat,
      design = svy_unadj
    )
  }, error = function(e) { cat("Error in HES unadj model, imp", imp_num, ":", e$message, "\n") })
  
  # --- Model 2: COVID Hospitalization (IPTW) ---
  tryCatch({
    fit_hes_iptw[[imp_num]] <- svycoxph(
      Surv(timeinstudy2, covid_hes_present) ~ treat,
      design = svy_iptw
    )
  }, error = function(e) { cat("Error in HES IPTW model, imp", imp_num, ":", e$message, "\n") })
  
  # --- Model 3: COVID Death (Unadjusted) ---
  tryCatch({
    fit_death_unadj[[imp_num]] <- svycoxph(
      Surv(timeinstudy3, covid_death_present) ~ treat,
      design = svy_unadj
    )
  }, error = function(e) { cat("Error in Death unadj model, imp", imp_num, ":", e$message, "\n") })
  
  # --- Model 4: COVID Death (IPTW) ---
  tryCatch({
    fit_death_iptw[[imp_num]] <- svycoxph(
      Surv(timeinstudy3, covid_death_present) ~ treat,
      design = svy_iptw
    )
  }, error = function(e) { cat("Error in Death IPTW model, imp", imp_num, ":", e$message, "\n") })
}

# POOLING RESULTS USING RUBIN'S RULES =========================================

cat("\nPooling results using Rubin's rules...\n")

# Convert the lists of fits to 'mira' objects, then pool.
# This ensures the pool() function recognizes the model objects correctly.
pooled_hes_unadj <- pool(as.mira(fit_hes_unadj))
pooled_hes_iptw  <- pool(as.mira(fit_hes_iptw))
pooled_death_unadj <- pool(as.mira(fit_death_unadj))
pooled_death_iptw  <- pool(as.mira(fit_death_iptw))

# Summarize the pooled results and exponentiate to get Hazard Ratios
summary_hes_unadj <- summary(pooled_hes_unadj, conf.int = TRUE, exponentiate = TRUE)
summary_hes_iptw  <- summary(pooled_hes_iptw,  conf.int = TRUE, exponentiate = TRUE)
summary_death_unadj <- summary(pooled_death_unadj, conf.int = TRUE, exponentiate = TRUE)
summary_death_iptw  <- summary(pooled_death_iptw,  conf.int = TRUE, exponentiate = TRUE)

# DIAGNOSTICS AND RESULTS ====================================================

# Balance assessment (remains the same)
cat("\nGenerating diagnostics...\n")
data_for_balance <- data %>% filter(!is.na(iptw))
balance_results <- bal.tab(x = ps_formula,
                           data = data_for_balance,
                           treat = "treat", 
                           weights = "iptw",
                           method = "weighting",
                           imp = ".imp",
                           stats = c("mean.diffs", "variance.ratios"))

cat("\nBalance assessment results:\n")
print(balance_results)


# Display the final pooled results in a clean table
cat("\n--- Final Pooled Results ---\n\n")

# Function to tidy up the summary output
clean_summary <- function(summary_obj, model_name) {
  df <- summary_obj[, c("term", "estimate", "2.5 %", "97.5 %", "p.value")]
  colnames(df) <- c("Term", "Hazard Ratio", "Lower 95% CI", "Upper 95% CI", "P-Value")
  df$Model <- model_name
  df <- df[df$Term == 'treatICS/LABA', ] # Assuming this is the term of interest
  return(df)
}

# Create a final results table
results_table <- bind_rows(
  clean_summary(summary_hes_unadj, "COVID HES - Unadjusted"),
  clean_summary(summary_hes_iptw,  "COVID HES - IPTW"),
  clean_summary(summary_death_unadj, "COVID Death - Unadjusted"),
  clean_summary(summary_death_iptw,  "COVID Death - IPTW")
)

cat("## Pooled Hazard Ratios for Treatment Effect ##\n")
print(results_table, row.names = FALSE)

end_time <- Sys.time()
cat("\nTotal analysis time:", format(end_time - start_time), "\n")