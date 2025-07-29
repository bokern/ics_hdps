# -----------------------------------------------------------------------------
# PROGRAM NAME:  05_effect_estimation
# PROJECT:      xzx26005
# AUTHOR:        Marleen Bokern
# DATE CREATED:  Aug 2024
# NOTES:        high-dimensional Propensity Score Analysis, estimating treatment effects for unweighted, HDPS-weighted, and predefined-covariate-weighted analyses

# -----------------------------------------------------------------------------
# Flag to enable debugging; if TRUE, script will only be run on 1,000 patients
# selected at random
debugMode = F

# Load libraries -------------------------------------------------------------

packages <- c(
  "tidyverse",
  "arrow",
  "dbplyr",
  "lubridate",
  "scales",
  "tableone",
  "Matching",
  "survival",
  "gridExtra",
  "survminer",
  "survey",
  "MetBrewer"
)

installed_packages <- packages %in% rownames(installed.packages())
if (any(installed_packages == FALSE)) {
  install.packages(packages[!installed_packages])
}

invisible(lapply(packages, function(pkg) {
  suppressMessages(library(pkg, character.only = TRUE, verbose = FALSE))
}))

set.seed(123)
setwd(Datadir_copd)
palette <- met.brewer("Cassatt2")

# Function to calculate risk differences with confidence intervals
calculate_risk_difference <- function(data, outcome_var, treatment_var, weights = NULL, alpha = 0.05) {
  
  # Convert to factors if needed
  data[[treatment_var]] <- as.factor(data[[treatment_var]])
  
  if (is.null(weights)) {
    # Unweighted analysis
    risk_table <- data %>%
      group_by(!!sym(treatment_var)) %>%
      summarise(
        n = n(),
        events = sum(!!sym(outcome_var), na.rm = TRUE),
        risk = mean(!!sym(outcome_var), na.rm = TRUE),
        .groups = 'drop'
      )
    
    # Calculate standard errors for unweighted risks
    risk_table$se <- sqrt(risk_table$risk * (1 - risk_table$risk) / risk_table$n)
    
  } else {
    # Weighted analysis
    risk_table <- data %>%
      group_by(!!sym(treatment_var)) %>%
      summarise(
        n = n(),
        weighted_n = sum(!!sym(weights), na.rm = TRUE),
        events = sum(!!sym(outcome_var), na.rm = TRUE),
        weighted_events = sum(!!sym(outcome_var) * !!sym(weights), na.rm = TRUE),
        risk = sum(!!sym(outcome_var) * !!sym(weights), na.rm = TRUE) / sum(!!sym(weights), na.rm = TRUE),
        .groups = 'drop'
      )
    
    # Calculate standard errors for weighted risks using effective sample size
    risk_table$se <- sqrt(risk_table$risk * (1 - risk_table$risk) / risk_table$weighted_n)
  }
  
  # Extract risks for exposed (1) and unexposed (0) groups
  risk_unexposed <- risk_table$risk[risk_table[[treatment_var]] == "0"]
  risk_exposed <- risk_table$risk[risk_table[[treatment_var]] == "1"]
  
  # Extract standard errors
  se_unexposed <- risk_table$se[risk_table[[treatment_var]] == "0"]
  se_exposed <- risk_table$se[risk_table[[treatment_var]] == "1"]
  
  # Calculate risk difference
  risk_difference <- risk_exposed - risk_unexposed
  
  # Calculate standard error of risk difference
  se_rd <- sqrt(se_unexposed^2 + se_exposed^2)
  
  # Calculate confidence interval
  z_alpha <- qnorm(1 - alpha/2)
  ci_lower <- risk_difference - z_alpha * se_rd
  ci_upper <- risk_difference + z_alpha * se_rd
  
  # Calculate p-value (two-sided test)
  z_stat <- risk_difference / se_rd
  p_value <- 2 * (1 - pnorm(abs(z_stat)))
  
  # Return results
  return(list(
    risk_unexposed = risk_unexposed,
    risk_exposed = risk_exposed,
    risk_difference = risk_difference,
    se_rd = se_rd,
    ci_lower = ci_lower,
    ci_upper = ci_upper,
    p_value = p_value
  ))
}

# Load data -------------------------------------------------------------------
all_schoenfeld_residuals <- list()

topVars <- c("100", "250", "500", "750", "1000")
outcomes <- c("covid_hes_present", "covid_death_present")
cohort_exts <- c("", "_no_triple")

for (cohort_ext in cohort_exts) { # Iterate through each cohort
  for (outcome in outcomes) {   # Iterate through each outcome
    
# cohort_ext <- ""
# outcome <- "covid_hes_present"

# Initialize risk difference results dataframe
risk_diff_results <- data.frame(
  Covariates = character(),
  Risk_Unexposed = numeric(),
  Risk_Exposed = numeric(),
  Risk_Difference = numeric(),
  SE_RD = numeric(),
  CI_Lower = numeric(),
  CI_Upper = numeric(),
  P_Value = numeric(),
  OutcomeEvents = numeric(),
  stringsAsFactors = FALSE
)

    print(paste("Analyzing outcome:", outcome))
    
    timeinstudy <- switch(outcome,
                          "covid_hes_present" = "timeinstudy2",
                          "covid_death_present" = "timeinstudy3")
    
    outcome_label <- switch(
      outcome,
      covid_hes_present = "COVID-19 Hospitalisation",
      covid_death_present = "COVID-19 Death"
    )
    
    # Load in patient summary
    pat_summary_data <- read_parquet("copd_wave1_60d.parquet") %>%
      mutate_at(c("patid"), as.character) %>%
      #remove if treatgroup is NA
      filter(!is.na(treatgroup))
    
    follow_up_times <- pat_summary_data %>% dplyr::select("patid", timeinstudy)
    
    if (cohort_ext == "_no_triple") {
      pat_summary_data <- pat_summary_data %>%
        filter(baseline_triple == 0)
    }
    
    #exclude if imd is "Missing"
    pat_summary_data <- pat_summary_data %>%
      filter(imd != "Missing")
    
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
      "treat",
      "baseline_triple",
      "covid_hes_present",
      "covid_death_present",
      "any_death_present",
      "any_covid_hes_present",
      "covid_present",
      "obese4cat",
      "exacerbations",
      "obese4cat",
      "bmi",
      "smokstatus"
    )
    
    # Unweighted analysis and analysis using predefined variables -------------
    PredefinedVars <- pat_summary_data %>% dplyr::select(-any_of(exclude_vars))
    
    file_path_predefined <- file.path(
      HDPS_folder,
      "outputs",
      paste0(
        ifelse(debugMode, "debug_", ""),
        "HDPS_predefined",
        "_",
        outcome,
        cohort_ext,
        "_trimmed.parquet"
      )
    )
    
    #read in parquet file
    hdpsPredefinedVars <- read_parquet(file_path_predefined)
    
    #merge in outcome
    hdpsPredefinedVars <- hdpsPredefinedVars %>%
      left_join(pat_summary_data %>% dplyr::select(patid, .data[[outcome]]), by = "patid")
    
    hdpsPredefinedVars$treatgroup <- relevel(hdpsPredefinedVars$treatgroup, ref = "0")
    
    # Estimate Treatment Effect (unweighted)
    surv_obj <- Surv(time = hdpsPredefinedVars[[timeinstudy]], event = hdpsPredefinedVars[[outcome]])
    # Unweighted Cox regression
    PredefinedUnweighted <- coxph(surv_obj ~ treatgroup, data = hdpsPredefinedVars)
    
    #plot schoenfeld residuals
    schoenfeld_residuals <- cox.zph(PredefinedUnweighted)
    
    all_schoenfeld_residuals[[paste("Unweighted", cohort_ext, outcome)]] <- data.frame(
      model = "Unweighted",
      cohort = cohort_ext,
      outcome = outcome,
      covariate = rownames(schoenfeld_residuals$table),
      p_value = schoenfeld_residuals$table[, "p"],
      stringsAsFactors = FALSE
    )

    # Kaplan-Meier Plot (Unweighted)
    km_fit_unweighted <- survfit(surv_obj ~ treatgroup, data = hdpsPredefinedVars)
    
    
    # Unweighted Logistic Regression
    logistic_unweighted <- glm(as.formula(paste(outcome, "~ treatgroup")),
                               data = hdpsPredefinedVars,
                               family = "binomial")
    
    # Calculate the residuals
    residuals <- residuals(logistic_unweighted)
    fitted_values <- fitted.values(logistic_unweighted)

    
    hdps_weighted_predefined <- coxph(
      surv_obj ~ treatgroup,
      data = hdpsPredefinedVars,
      weights = hdpsPredefinedVars$iptw_weight_predefined
    )
    
    # Kaplan-Meier Plot
    km_fit_weighted_predefined <- survfit(
      surv_obj ~ treatgroup,
      data = hdpsPredefinedVars,
      weights = hdpsPredefinedVars$iptw_weight_predefined
    )
    
    # Weighted Logistic Regression (Predefined)
    logistic_weighted_predefined <- glm(
      as.formula(paste(outcome, "~ treatgroup")),
      data = hdpsPredefinedVars,
      family = binomial(link = "logit"),
      weights = iptw_weight_predefined
    )
    
    # Calculate the residuals
    residuals <- residuals(logistic_weighted_predefined)
    fitted_values <- fitted.values(logistic_weighted_predefined)
    

    #get number of events
    num_events <- hdpsPredefinedVars %>%
      filter(!!sym(outcome) == 1) %>%
      nrow()
    #get number of weighted events
    num_weighted_events <- hdpsPredefinedVars %>%
      filter(!!sym(outcome) == 1) %>%
      summarise(weighted_events = sum(iptw_weight_predefined)) %>%
      pull()
    
    all_results <- data.frame(
      Covariates = character(),
      Results = character(),
      OutcomeEvents = numeric(),
      stringsAsFactors = FALSE
    )
    
    logistic_results <- data.frame(
      Covariates = character(),
      Results = character(),
      OutcomeEvents = numeric(),
      stringsAsFactors = FALSE)
    
    # Add the initial "Unweighted" and "Predefined" results as before
    all_results <- rbind(all_results,
                         data.frame(
                           Covariates = "Unweighted",
                           Results = ShowRegTable(PredefinedUnweighted, printToggle = FALSE),
                           OutcomeEvents = num_events))
    
    all_results <- rbind(all_results,
                         data.frame(
                           Covariates = "Prespecified covariates",
                           Results = ShowRegTable(hdps_weighted_predefined, printToggle = FALSE),
                           OutcomeEvents = num_weighted_events))
    
    logistic_results <- rbind(logistic_results,
                              data.frame(
                                Covariates = "Unweighted",
                                Results = ShowRegTable(logistic_unweighted, printToggle = FALSE),
                                OutcomeEvents = num_events))
    
    logistic_results <- rbind(logistic_results,
                              data.frame(
                                Covariates = "Prespecified covariates",
                                Results = ShowRegTable(logistic_weighted_predefined, printToggle = FALSE),
                                OutcomeEvents = num_weighted_events))
    
    # Calculate risk differences for unweighted analysis
    rd_unweighted <- calculate_risk_difference(
      data = hdpsPredefinedVars,
      outcome_var = outcome,
      treatment_var = "treatgroup"
    )
    
    risk_diff_results <- rbind(risk_diff_results, data.frame(
      Covariates = "Unweighted",
      Risk_Unexposed = rd_unweighted$risk_unexposed,
      Risk_Exposed = rd_unweighted$risk_exposed,
      Risk_Difference = rd_unweighted$risk_difference,
      SE_RD = rd_unweighted$se_rd,
      CI_Lower = rd_unweighted$ci_lower,
      CI_Upper = rd_unweighted$ci_upper,
      P_Value = rd_unweighted$p_value,
      OutcomeEvents = num_events
    ))
    
    # Calculate risk differences for predefined covariates weighted analysis
    rd_predefined <- calculate_risk_difference(
      data = hdpsPredefinedVars,
      outcome_var = outcome,
      treatment_var = "treatgroup",
      weights = "iptw_weight_predefined"
    )
    
    risk_diff_results <- rbind(risk_diff_results, data.frame(
      Covariates = "Prespecified covariates",
      Risk_Unexposed = rd_predefined$risk_unexposed,
      Risk_Exposed = rd_predefined$risk_exposed,
      Risk_Difference = rd_predefined$risk_difference,
      SE_RD = rd_predefined$se_rd,
      CI_Lower = rd_predefined$ci_lower,
      CI_Upper = rd_predefined$ci_upper,
      P_Value = rd_predefined$p_value,
      OutcomeEvents = num_weighted_events
    ))
    
    for (k in topVars) {
      
      print(paste("Analyzing top", k, "variables"))
      
      # Read in parquet file
      hdpsData <- read_parquet(
        paste0(
          HDPS_folder,
          "/outputs/HDPS_",
          k,
          "_",
          outcome,
          cohort_ext,
          "_trimmed.parquet"
        )
      )
      
      #move patid, pscore and iptw to the front
      hdpsData <- hdpsData %>%
        dplyr::select(patid, pscore, iptw_weight, everything())
      
      # Estimate Treatment Effect -----------------------------------------------
      #merge follow-up times in if not present
      if (!(timeinstudy %in% names(hdpsData))) {
        hdpsData <- hdpsData %>%
          left_join(follow_up_times, by = c("patid" = "patid"))
      }
      
      if (!("timeinstudy" %in% names(hdpsPredefinedVars))) {
        hdpsPredefinedVars <- hdpsPredefinedVars %>%
          left_join(follow_up_times, by = c("patid" = "patid"))
      }
      
      if (!(outcome %in% names(hdpsPredefinedVars))) {
        hdpsPredefinedVars <- hdpsPredefinedVars %>%
          left_join(hdpsData %>% dplyr::select("patid", outcome),
                    by = c("patid" = "patid"))
      }
      
      hdpsData[[timeinstudy]] <- as.numeric(hdpsData[[timeinstudy]])
      
      hdpsData$treatgroup <- relevel(hdpsData$treatgroup, ref = "0")
      
      # Estimate Treatment Effect for predefined covariates
      surv_obj <- Surv(time = hdpsData[[timeinstudy]], event = hdpsData[[outcome]])
      
      # Weighted Analysis (HDPS)
      hdpsWeighted <- coxph(
        surv_obj ~ hdpsData$treatgroup,
        data = hdpsData,
        weights = iptw_weight,
        robust = TRUE
      )
      
      #get weighted event count
      event_count <- hdpsData %>%
        filter(!!sym(outcome) == 1) %>%
        summarise(n = sum(iptw_weight)) %>%
        pull(n)
      
      # Add results to the all_results dataframe
      all_results <- rbind(all_results,
                           data.frame(
                             Covariates = as.character(k),
                             Results = ShowRegTable(hdpsWeighted, printToggle = FALSE),
                             OutcomeEvents = event_count
                           ))
      
      # Weighted Logistic Regression (HDPS)
      logistic_weighted <- glm(
        as.formula(paste(outcome, "~ treatgroup")),
        data = hdpsData,
        family = binomial(link = "logit"),
        weights = iptw_weight
      )
      
      # Calculate the residuals
      residuals <- residuals(logistic_weighted)
      fitted_values <- fitted.values(logistic_weighted)
      
      
      # Add results to the logistic_results dataframe
      logistic_results <- rbind(
        logistic_results,
        data.frame(
          Covariates = as.character(k),
          Results = ShowRegTable(logistic_weighted, printToggle = FALSE),
          OutcomeEvents = event_count
        )
      )
      
      # remove rows with the rowname Intercept from the logistic regression results
      logistic_results <- logistic_results[!grepl("Intercept", rownames(logistic_results)), ]
      
      # Calculate risk differences for HDPS weighted analysis
      rd_hdps <- calculate_risk_difference(
        data = hdpsData,
        outcome_var = outcome,
        treatment_var = "treatgroup",
        weights = "iptw_weight"
      )
      
      risk_diff_results <- rbind(risk_diff_results, data.frame(
        Covariates = as.character(k),
        Risk_Unexposed = rd_hdps$risk_unexposed,
        Risk_Exposed = rd_hdps$risk_exposed,
        Risk_Difference = rd_hdps$risk_difference,
        SE_RD = rd_hdps$se_rd,
        CI_Lower = rd_hdps$ci_lower,
        CI_Upper = rd_hdps$ci_upper,
        P_Value = rd_hdps$p_value,
        OutcomeEvents = event_count
      ))
      
    }
    
    # Clean up
    rm(list = ls()[!ls() %in% c(
      "Datadir_copd",
      "HDPS_folder",
      "hdpsResults",
      "debugMode",
      "k",
      "exclude_triple",
      "outcome",
      "topVars",
      "hdpsCovariates",
      "hdpsPredefinedVars",
      "tabWeighted_predefined",
      "hdpsWeighted_predefined",
      "cohort_ext",
      "hdpsCohort",
      "timeinstudy",
      "palette",
      "exposed",
      "outcomes",
      "follow_up_times",
      "all_results",
      "logistic_results",
      "smd_combined",
      "HDPS_github",
      "all_schoenfeld_residuals",
      "risk_diff_results",
      "calculate_risk_difference"
    )])
    
    # Output data -------------------------------------------------------------
    # combine all results and write to a single file
    combined_results <- do.call(rbind, all_results)
    logistic_combined_results <- do.call(rbind, logistic_results)
    
    # After processing all k values, write the results to a CSV file
    write_csv(
      all_results,
      paste0(
        HDPS_folder,
        "/outputs/Results/HDPS_Cox_Results_",
        outcome,
        cohort_ext,
        ".csv"
      )
    )
    
    # Write logistic regression results to a separate CSV file
    write_csv(
      logistic_results,
      paste0(
        HDPS_folder,
        "/outputs/Results/HDPS_Logistic_Results_",
        outcome,
        cohort_ext,
        ".csv"
      )
    )
    
    # Write risk difference results to CSV file
    write_csv(
      risk_diff_results,
      paste0(
        HDPS_folder,
        "/outputs/Results/HDPS_Risk_Differences_",
        outcome,
        cohort_ext,
        ".csv"
      )
    )
  }
}

# Save the Schoenfeld residuals to a CSV file
all_schoenfeld_residuals <- do.call(rbind, all_schoenfeld_residuals)
write_csv(
  all_schoenfeld_residuals,
  paste0(HDPS_folder, "/outputs/Schoenfeld_Residuals.csv")
)
