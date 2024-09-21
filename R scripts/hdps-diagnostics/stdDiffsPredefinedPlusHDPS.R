# -----------------------------------------------------------------------------
# PROGRAM NAME:  006_stdDiffsPredefinedPlusHDPS.R
# PROJECT:      xzx26005
# AUTHOR:        John Tazare
# DATE CREATED:   10 Sep 2020
# NOTES:        high-dimensional Propensity Score Analysis for ACNU:
#                                  1. Summarise baseline covariates
#                                  2. Estimate propensity score
#                                  3. Propensity score matching and evaluation
#                                  4. Estimate the treatment effect
#
# REQUIRES: acnuCohort.csv & acnuFollowUp.csv
# -----------------------------------------------------------------------------

packages <- c("tidyverse", "ggplot2", "MetBrewer")
installed_packages <- packages %in% rownames(installed.packages())
if (any(installed_packages == FALSE)) {
  install.packages(packages[!installed_packages])
}

invisible(lapply(packages, function(pkg) {
  suppressMessages(library(pkg, character.only = TRUE, verbose = FALSE))
}))

set.seed(123)
setwd(HDPS_folder)

palette <- met.brewer("Cassatt2")

topVars <- c("100", "250", "500", "750", "1000")
outcomes <- c("covid_hes_present", "covid_death_present")
cohort_exts <- c("", "_no_triple")

for (cohort_ext in cohort_exts) {
  for (outcome in outcomes) {
    for (k in topVars) {

      print(paste0("Processing HDPS_", k, "_SMD_", outcome, cohort_ext ))
      data_path <- file.path(HDPS_folder, paste0("/outputs/HDPS_", k, "_SMD_", outcome, cohort_ext, "_trimmed.csv"))
      # Read the data
      data <- read_csv(data_path)
      
      # remove the first column
      data <- data[, -1]
      
      # Ensure smd columns are numeric
      data <- data %>%
        mutate(across(c(smd_unweighted, smd_weighted), as.numeric))
      
      #if variable called n is present remove it
      data <- data %>% filter(!str_detect(Variable, "^n$"))
      
      # Define covariate type based on variable name
      data <- data %>%
        mutate(covariate_type = if_else(str_detect(Variable, "^d\\d+_"), "HDPS", "Predefined")) %>% 
        group_by(covariate_type) %>%
        arrange(desc(covariate_type), desc(smd_unweighted)) %>%
        ungroup()
      
      #if smd_unweighted and smd_weighted are na remove the row
      data <- data %>% filter(!is.na(smd_unweighted) & !is.na(smd_weighted))
      
      # Get the number of predefined and HDPS covariates
      n_predefined <- sum(data$covariate_type == "Predefined")
      n_hdps <- sum(data$covariate_type == "HDPS")
      
      # Modify the data preparation steps
      data <- data %>%
        mutate(covariate_type = if_else(str_detect(Variable, "^d\\d+_"), "HDPS", "Predefined")) %>% 
        group_by(covariate_type) %>%
        arrange(desc(covariate_type), desc(smd_unweighted)) %>%
        ungroup()
      
      # Create custom labels for predefined covariates (in reverse order)
      n_predefined <- sum(data$covariate_type == "Predefined")
      
      #if Variable includes the string "past_asthma", change label to "Past asthma"
      data <- data %>%
        mutate(var_label = case_when(
          str_detect(Variable, "past_asthma") ~ "Past asthma",
          str_detect(Variable, "exacerb") ~ "Exacerbation",
          str_detect(Variable, "pneumo_vac") ~ "Pneumococcal vaccination",
          str_detect(Variable, "smokstatus") ~ "Smoking status",
          str_detect(Variable, "bmicat") ~ "BMI",
          str_detect(Variable, "eth") ~ "Ethnicity",
          str_detect(Variable, "age_index") ~ "Age",
          str_detect(Variable, "gender") ~ "Gender",
          str_detect(Variable, "allcancers") ~ "Cancer",
          str_detect(Variable, "cvd") ~ "CVD",
          str_detect(Variable, "diabetes") ~ "Diabetes",
          str_detect(Variable, "imd") ~ "IMD",
          str_detect(Variable, "flu_vac") ~ "Flu vaccination",
          str_detect(Variable, "kidney") ~ "Kidney disease",
          str_detect(Variable, "hypertension") ~ "Hypertension",
          str_detect(Variable, "immunosuppression") ~ "Immunosuppression",
          TRUE ~ ""
        ))
      
      # Update the adjusted rank calculation
      data <- data %>%
        mutate(adjusted_rank = case_when(
          covariate_type == "Predefined" ~ -row_number() * 1.5,  # Increase spacing
          covariate_type == "HDPS" ~ row_number() - n_predefined  # Adjust HDPS starting point
        ))
      
      # # Update custom breaks and labels
      # min_value <- min(data$adjusted_rank)  # This should be -24
      # max_value <- max(data$adjusted_rank)  # This should be 100
      # 
      # custom_breaks <- c(
      #   seq(min_value, -25, by = 25),  # Negative breaks: -24, -25
      #   seq(0, max_value, by = 25)     # Positive breaks: 0, 25, 50, 75, 100
      # )
      
      # Ensure var_label is correctly referenced here; you might need to adjust how you're creating custom_labels
      custom_labels <- unique(data$var_label) # Assuming you want unique labels from var_label
      
      # Create the main plot without x-axis labels for HDPS covariates
      p <- ggplot(data) +
        geom_point(aes(x = adjusted_rank, y = smd_unweighted, color = "Unweighted"), 
                   shape = 16, alpha = 1, size = 2) +
        geom_point(aes(x = adjusted_rank, y = smd_weighted, color = "HDPS-Weighted"), 
                   shape = 18, alpha = 1, size = 2) +
        geom_point(aes(x = adjusted_rank, y = smd_weighted_predefined, color = "Predefined Covariates"), 
                   shape = 17, alpha = 1, size = 2) +
        scale_x_continuous(
          expand = c(0.02, 0.02)
        ) +
        scale_y_reverse(limits = c(max(c(data$smd_unweighted, data$smd_weighted), na.rm = TRUE) * 1.1, 0), 
                        breaks = seq(0, 1, by = 0.1)) +
        labs(
          x = "Covariates",
          y = "Absolute standardised mean difference",
          color = "Weighting"
        ) +
        geom_hline(yintercept = c(0, 0.1), linetype = "dashed", color = "black", size = 0.5) +
        geom_vline(xintercept = 0, linetype = "dashed", color = "black", size = 0.5) +
        annotate("text", x = min(data$adjusted_rank)/2, y = max(c(data$smd_unweighted, data$smd_weighted), na.rm = TRUE) * 1.15, 
                 label = "Predefined Covariates", size = 3.5) +
        annotate("text", x = max(data$adjusted_rank)/2, y = max(c(data$smd_unweighted, data$smd_weighted), na.rm = TRUE) * 1.15, 
                 label = "HDPS covariates", size = 3.5) +
        theme_minimal() +
        theme(
          legend.position = c(0.85, 0.2),
          legend.background = element_rect(fill = "white", color = NA),
          legend.title = element_text(size = 10),
          legend.text = element_text(size = 8),
          axis.ticks.x = element_blank(), # Remove x-axis ticks
          panel.grid.minor = element_blank()
        ) +
        scale_color_manual(values = c("Unweighted" = palette[4], "HDPS-Weighted" = palette[7], "Predefined Covariates" = palette[9])) +
        guides(color = guide_legend(override.aes = list(size = 3)))
      
      # Print the plot
      print(p)
      
      # Save
      file_path <- file.path(
        HDPS_folder,
        "outputs",
        "diagnostics",
        paste0("SMD_", k, outcome, cohort_ext, ".png")
      )
      
      ggsave(file_path, p, width = 10, height = 8)
      
    }
  }
}


