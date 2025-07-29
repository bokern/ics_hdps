# -----------------------------------------------------------------------------
# PROGRAM NAME:  04_hdps_forest_plots
# PROJECT:      xzx26005
# AUTHOR:        John Tazare
# DATE CREATED:   10 Sep 2020
# NOTES:        high-dimensional Propensity Score Analysis for ACNU:
#                                  1. Summarise baseline covariates
#                                  2. Estimate propensity score
#                                  3. Propensity score matching and evaluation
#                                  4. Estimate the treatment effect
# -----------------------------------------------------------------------------

exclude_triple = T

if (exclude_triple == TRUE) {
  cohort_ext <- "_no_triple"
} else {
  cohort_ext <- ""
}

# Load libraries -------------------------------------------------------------

packages <- c("tidyverse", "arrow", "ggplot2", "patchwork")

installed_packages <- packages %in% rownames(installed.packages())
if (any(installed_packages == FALSE)) {
  install.packages(packages[!installed_packages])
}

invisible(lapply(packages, function(pkg) {
  suppressMessages(library(pkg, character.only = TRUE, verbose = FALSE))
}))

setwd(HDPS_folder)
models <- c("Logistic", "Cox")
cohort_exts <- c("", "_no_triple")
outcomes <- c("covid_hes_present", "covid_death_present")

  for (outcome in outcomes) {
    for (model in models) {
      
      # Read in results file
      results_file_main <- paste0(
        "outputs/Results/HDPS_",
        model,
        "_Results_",
        outcome,
        ".csv"
      )
      results_main <- read_csv(results_file_main)
      
      results_file_no_triple <- paste0(
        "outputs/Results/HDPS_",
        model,
        "_Results_",
        outcome,
        "_no_triple",
        ".csv"
      )
      results_no_triple <- read_csv(results_file_no_triple)
      
      # Modify covariate names  
      modify_covariates <- function(df) {
        df %>%
          mutate(
            Covariates = case_when(
              Covariates == "10" ~ "Prespecified + 10 HDPS",
              Covariates == "100" ~ "Prespecified + 100 HDPS",
              Covariates == "250" ~ "Prespecified + 250 HDPS",
              Covariates == "500" ~ "Prespecified + 500 HDPS",
              Covariates == "750" ~ "Prespecified + 750 HDPS",
              Covariates == "1000" ~ "Prespecified + 1000 HDPS",
              TRUE ~ Covariates
            )
          )
      }
      
      results_main <- results_main %>% modify_covariates()
      results_no_triple <- results_no_triple %>% modify_covariates()
      
      # Parse the results for hospitalization
      results_main <- results_main %>%
        mutate(
          estimate = as.numeric(str_extract(
            Results.exp.coef...confint., "\\d+\\.\\d+"
          )),
          lower_CI = as.numeric(
            str_extract(Results.exp.coef...confint., "(?<=\\[)\\d+\\.\\d+(?=,)")
          ),
          upper_CI = as.numeric(
            str_extract(
              Results.exp.coef...confint.,
              "(?<=,\\s)\\d+\\.\\d+(?=\\])"
            )
          ),
          row = factor(rev(row_number())),
          #events = OutcomeEvents, rounded to 1 decimal place
          events = round(OutcomeEvents, 1)
        )
      
      # Parse the results for death
      results_no_triple <- results_no_triple %>%
        mutate(
          estimate = as.numeric(str_extract(
            Results.exp.coef...confint., "\\d+\\.\\d+"
          )),
          lower_CI = as.numeric(
            str_extract(Results.exp.coef...confint., "(?<=\\[)\\d+\\.\\d+(?=,)")
          ),
          upper_CI = as.numeric(
            str_extract(
              Results.exp.coef...confint.,
              "(?<=,\\s)\\d+\\.\\d+(?=\\])"
            )
          ),
          row = factor(rev(row_number())),
          events = round(OutcomeEvents, 1)
        )
      
      # Create the forest plot for hospitalization
      forest_plot_main <- ggplot(results_main, aes(y = row)) +
        geom_point(aes(x = estimate), shape = 16, size = 3) +
        geom_linerange(aes(xmin = lower_CI, xmax = upper_CI)) +
        geom_vline(xintercept = 1, linetype = "dashed") +
        scale_x_continuous(trans = 'log10', breaks = c(0.8, 1, 2, 4)) +
        coord_cartesian(xlim = c(0.7, 4)) +
        labs(x = "Effect estimate", y = "", title = "With Triple Therapy") +
        theme_classic() +
        theme(
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.line.y = element_blank(), # This removes the y-axis line
          axis.line.x = element_line(color = "black"),# Ensure x-axis line is visible
          panel.border = element_blank(),# This removes the panel border
          plot.margin = margin(
            l = 0,
            r = 0,
            t = 0,
            b = 0,
            unit = "pt"
          )  # Adjust margin to add space on the left
        )
      
      label_plot_main <- ggplot(results_main, aes(y = row)) +
        geom_text(
          aes(x = 0, label = Covariates),
          hjust = 0,
          size = 5.5,
          fontface = "bold"
        ) +
        geom_text(aes(
          x = 1.2,  # Changed from 2.1 to 1
          label = paste(events)
        ),
        hjust = 0,
        size = 5.2) +
        geom_text(aes(
          x = 2.1,  
          label = sprintf("%.2f (%.2f - %.2f)", estimate, lower_CI, upper_CI)
        ),
        hjust = 1,  
        size = 5.2) +
        scale_x_continuous(limits = c(0, 2.1)) +  # Force x-axis to focus on our text area
        theme_void() +
        theme(plot.margin = margin(r = 0))
      
      # Create the forest plot for death
      forest_plot_no_triple <- ggplot(results_no_triple, aes(y = row)) +
        geom_point(aes(x = estimate), shape = 16, size = 3) +
        geom_linerange(aes(xmin = lower_CI, xmax = upper_CI)) +
        geom_vline(xintercept = 1, linetype = "dashed") +
        scale_x_continuous(trans = 'log10', breaks = c(0.8, 1, 2, 4)) +
        coord_cartesian(xlim = c(0.7, 4)) +
        labs(x = "Effect estimate", y = "", title = "Without Triple Therapy") +
        theme_classic() +
        theme(
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.line.y = element_blank(),
          # This removes the y-axis line
          axis.line.x = element_line(color = "black"),
          # Ensure x-axis line is visible
          panel.border = element_blank(),
          plot.margin = margin(
            l = 0,
            r = 0.1,
            t = 0,
            b = 0,
            unit = "pt"
          )  # Adjust margin to add space on the left
        )
      
      label_plot_no_triple <- ggplot(results_no_triple, aes(y = row)) +
        geom_text(aes(
          x = 0,
          label = paste(events)
        ),
        hjust = 0,
        size = 5.2) +
        geom_text(aes(
          x = 0.1,
          label = sprintf("%.2f (%.2f - %.2f)", estimate, lower_CI, upper_CI)
        ),
        hjust = 0,
        size = 5.2) +
        scale_x_continuous(limits = c(0, 0.3)) +  # Force x-axis to focus on our text area
        theme_void() +
        theme(plot.margin = margin(l = 0, r = 0, t = 0, b = 0, unit = "pt"))
      
      # Combine plots with adjusted layout
      combined_plot <- label_plot_main + forest_plot_main +
        label_plot_no_triple + forest_plot_no_triple +
        plot_layout(widths = c(2.3, 1.2, 1.3, 1.2))
      
      # Display the plot
      combined_plot
      
      ggsave(
        paste0(
          "outputs/Results/",
          model,
          "_forest_plots_",
          outcome,
          ".png"
        ),
        plot = combined_plot,
        width = 16,
        height = 7
      )
    }
  }