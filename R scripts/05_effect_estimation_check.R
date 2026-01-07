# -----------------------------------------------------------------------------
# PROGRAM NAME:  05_effect_estimation
# PROJECT:      xzx26005
# AUTHOR:        Marleen Bokern
# DATE CREATED:  Aug 2024
# NOTES:        high-dimensional Propensity Score Analysis for ACNU:
#                                  1. Summarise baseline covariates
#                                  2. Estimate propensity score
#                                  3. Propensity score matching and evaluation
#                                  4. Estimate the treatment effect

# -----------------------------------------------------------------------------
# Flag to enable debugging; if TRUE, script will only be run on 1,000 patients
# selected at random
debugMode = F

exclude_triple = T

if (exclude_triple == TRUE) {
  cohort_ext <- "_no_triple"
} else {
  cohort_ext <- ""
}
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

# Load data -------------------------------------------------------------------
all_schoenfeld_residuals <- list()

topVars <- seq(1, 250, 1)
outcome <- "covid_hes_present"
cohort_ext <- "_no_triple"

print(paste("Analyzing outcome:", outcome))

timeinstudy <- switch(outcome,
                      "covid_hes_present" = "timeinstudy2",
                      "covid_death_present" = "timeinstudy3")

# Load in patient summary
pat_summary_data <- read_parquet("copd_wave1_60d.parquet") %>%
  mutate_at(c("patid"), as.character) %>%
  #remove if treatgroup is NA
  filter(!is.na(treatgroup))

follow_up_times <- pat_summary_data %>% dplyr::select("patid", timeinstudy)


if (exclude_triple == TRUE) {
  pat_summary_data <- pat_summary_data %>%
    filter(baseline_triple == 0)
}

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
    "_trimmed_check.parquet"
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
# Unweighted Analysis
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

png(
  filename = paste0(
    HDPS_folder,
    "/outputs/Results/Schoenfeld_",
    outcome,
    "_unweighted",
    cohort_ext,
    "_check.png"
  )
)
plot(schoenfeld_residuals)
dev.off()

# Kaplan-Meier Plot (Unweighted)
km_fit_unweighted <- survfit(surv_obj ~ treatgroup, data = hdpsPredefinedVars)

# Unweighted Logistic Regression
logistic_unweighted <- glm(as.formula(paste(outcome, "~ treatgroup")),
                           data = hdpsPredefinedVars,
                           family = binomial(link = "logit"))

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

#plot schoenfeld residuals
schoenfeld_residuals <- cox.zph(hdps_weighted_predefined)

all_schoenfeld_residuals[[paste("Predefined", cohort_ext, outcome)]] <- data.frame(
  model = "predefined",
  cohort = cohort_ext,
  outcome = outcome,
  covariate = rownames(schoenfeld_residuals$table),
  p_value = schoenfeld_residuals$table[, "p"],
  stringsAsFactors = FALSE
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

all_results <- data.frame(
  Covariates = character(),
  Results = character(),
  stringsAsFactors = FALSE)

logistic_results <- data.frame(
  Covariates = character(),
  Results = character(),
  stringsAsFactors = FALSE)

# Add the initial "Unweighted" and "Predefined" results as before
all_results <- rbind(all_results,
                     data.frame(
                       Covariates = "Unweighted",
                       Results = ShowRegTable(PredefinedUnweighted, printToggle = FALSE)))
all_results <- rbind(all_results,
                     data.frame(
                       Covariates = "Predefined",
                       Results = ShowRegTable(hdps_weighted_predefined, printToggle = FALSE)))

logistic_results <- rbind(logistic_results,
                          data.frame(
                            Covariates = "Unweighted",
                            Results = ShowRegTable(logistic_unweighted, printToggle = FALSE)))
logistic_results <- rbind(logistic_results,
                          data.frame(
                            Covariates = "Predefined",
                            Results = ShowRegTable(logistic_weighted_predefined, printToggle = FALSE)))

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
      "_trimmed_check.parquet"
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
  
  # Kaplan-Meier Plot (Unweighted)
  km_fit_weighted_HDPS <- survfit(surv_obj ~ treatgroup,
                                  data = hdpsData,
                                  weights = iptw_weight)
  
  # Plot the Kaplan-Meier curve
  ggsurvplot <- ggsurvplot(
    km_fit_weighted_HDPS,
    data = hdpsData,
    conf.int = T,
    censor = F,
    ylim = c(0.97, 1),
    xlab = "Time in days",
    risk.table = "absolute",
    risk.table.title = "Number at risk",
    cumevents = TRUE,
    fontsize = 7,
    tables.height = 0.15,
    legend.labs = c("LABA/LAMA", "ICS"),
    legend.title = "",
    palette = c(palette[4], palette[9]),
    xlim = c(0, 183)
  )
  ggsurvplot$cumevents$layers[[1]]$data$cum.n.event <- round(ggsurvplot$cumevents$layers[[1]]$data$cum.n.event, 2)
  ggsurvplot$plot <- ggsurvplot$plot + scale_x_continuous(breaks = c(0, 50, 100, 150, 183))
  ggsurvplot$table$theme$axis.text.y$colour <- "black"
  ggsurvplot$table$theme$axis.text.y$size <- 24
  ggsurvplot$table$theme$axis.text.x$size <- 24
  ggsurvplot$table$labels$x <- ""
  ggsurvplot$table$theme$plot.title$size <- 28
  ggsurvplot$cumevents$theme$axis.text.y$colour <- "black"
  ggsurvplot$cumevents$theme$axis.text.y$size <- 24
  ggsurvplot$cumevents$theme$axis.text.x$size <- 24
  ggsurvplot$cumevents$labels$x <- ""
  ggsurvplot$cumevents$theme$plot.title$size <- 28
  ggsurvplot$plot$theme$axis.title.x$size <- 28
  ggsurvplot$plot$theme$axis.title.y$size <- 28
  ggsurvplot$plot$theme$axis.text.x$size <- 28
  ggsurvplot$plot$theme$axis.text.y$size <- 28
  ggsurvplot$plot$theme$legend.text$size <- 28
  ggsurvplot$plot$theme$legend.key.size <- unit(4, "lines")
  
  # Save the main plot
  file_path <- file.path(
    HDPS_folder,
    "outputs",
    paste0("km_curve_", outcome, "_", k, "_HDPS", cohort_ext, "_check.png")
  )
  # Combine the main plot, risk table, and cumulative events
  combined_plot <- ggarrange(
    ggsurvplot$plot, 
    ggsurvplot$table, 
    ggsurvplot$cumevents,
    nrow = 3,
    heights = c(2, 0.6, 0.6)  # Adjust height ratios as needed
  )
  ggsave(
    file_path,
    combined_plot,
    width = 30,
    height = 30,
    units = "cm")
  
  
  #plot schoenfeld residuals
  schoenfeld_residuals <- cox.zph(hdpsWeighted)
  
  all_schoenfeld_residuals[[paste0("HDPS_", k, cohort_ext, outcome)]] <- data.frame(
    model = paste0("HDPS", k),
    cohort = cohort_ext,
    outcome = outcome,
    covariate = rownames(schoenfeld_residuals$table),
    p_value = schoenfeld_residuals$table[, "p"],
    stringsAsFactors = FALSE
  )
  
  png(
    filename = paste0(
      HDPS_folder,
      "/outputs/Results/Schoenfeld_",
      outcome,
      "_",
      k,
      "_HDPS",
      cohort_ext,
      "_check.png"
    )
  )
  plot(schoenfeld_residuals)
  dev.off()
  
  # Add results to the all_results dataframe (Cox model)
  all_results <- rbind(all_results,
                       data.frame(
                         Covariates = as.character(k),
                         Results = ShowRegTable(hdpsWeighted, printToggle = FALSE)
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
  
  file_path_resid <- file.path(
    HDPS_folder,
    "outputs",
    "Results",
    paste0("log_resid_", outcome, "_HDPS_", k, cohort_ext, "_jittered_check.png")
  )
  png(file_path_resid)
  
  # Create a plot of the jittered residuals versus the jittered fitted values
  plot(
    jitter(fitted_values, factor = 0.3),
    jitter(residuals, factor = 10),
    main = paste(
      "Log residuals for",
      outcome,
      "(weighted using",
      k,
      "covariates)"
    ),
    cex.main = 0.9,
    xlab = "Fitted values",
    ylab = "Residuals",
    pch = 16,
    # Use solid circles for points
    cex = 0.5  # Reduce point size for better visibility
  )
  abline(h = 0, lty = 2)
  
  dev.off()
  
  # Add results to the logistic_results dataframe
  logistic_results <- rbind(
    logistic_results,
    data.frame(
      Covariates = as.character(k),
      Results = ShowRegTable(logistic_weighted, printToggle = FALSE)
    )
  )
  
  # remove rows with the rowname Intercept from the logistic regression results
  logistic_results <- logistic_results[!grepl("Intercept", rownames(logistic_results)), ]
  
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
  "all_schoenfeld_residuals"
)])

# Output data -------------------------------------------------------------
# After the loop for k, combine all results and write to a single file
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
    "_check.csv"
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
    "_check.csv"
  )
)

# Save the Schoenfeld residuals to a CSV file
all_schoenfeld_residuals <- do.call(rbind, all_schoenfeld_residuals)
write_csv(
  all_schoenfeld_residuals,
  paste0(HDPS_folder, "/outputs/Schoenfeld_Residuals_check.csv")
)
