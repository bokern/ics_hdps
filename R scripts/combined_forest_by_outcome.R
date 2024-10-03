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

# outcome <- "covid_hes_present"
# model <- "Logistic"
#     
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
      "_no_triple.csv"
    )
    results_no_triple <- read_csv(results_file_no_triple)
    
    modify_covariates <- function(df) {
      df %>%
        mutate(
          Covariates = case_when(
            Covariates == "10" ~ "Predefined + 10 HDPS",
            Covariates == "100" ~ "Predefined + 100 HDPS",
            Covariates == "250" ~ "Predefined + 250 HDPS",
            Covariates == "500" ~ "Predefined + 500 HDPS",
            Covariates == "750" ~ "Predefined + 750 HDPS",
            Covariates == "1000" ~ "Predefined + 1000 HDPS",
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
        row = factor(rev(row_number()))
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
        row = factor(rev(row_number()))
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
        axis.line.y = element_blank(),
        # This removes the y-axis line
        axis.line.x = element_line(color = "black"),
        # Ensure x-axis line is visible
        panel.border = element_blank(),
        # This removes the panel border
        plot.margin = margin(
          l = 2,
          r = 0,
          t = 0,
          b = 0,
          unit = "pt"
        )  # Adjust margin to add space on the left
      )
    
    label_plot_main <- ggplot(results_main, aes(y = row)) +
      geom_text(
        aes(x = 1, label = Covariates),
        hjust = 0,
        size = 5,
        fontface = "bold"
      ) +
      geom_text(aes(
        x = 3.4,
        label = sprintf("%.2f (%.2f - %.2f)", estimate, lower_CI, upper_CI)
      ),
      hjust = 1.2,
      size = 5) +
      theme_void() +
      theme(plot.margin = margin(r = 5)) # Increase right margin
    
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
          l = 10,
          r = 0,
          t = 0,
          b = 0,
          unit = "pt"
        )  # Adjust margin to add space on the left
      )
    
    # Create the label plot for death
    label_plot_no_triple <- ggplot(results_no_triple, aes(y = row)) +
      geom_text(aes(
        x = 3,
        label = sprintf("%.2f (%.2f - %.2f)", estimate, lower_CI, upper_CI)
      ),
      hjust = 0.35,
      size = 5) +
      theme_void() +
      theme(plot.margin = margin(
        l = 10,
        r = 0,
        t = 0,
        b = 0,
        unit = "pt"
      ))
    
    # Combine plots with adjusted layout
    combined_plot <- label_plot_main + plot_spacer() + forest_plot_main + plot_spacer() +
      label_plot_no_triple + plot_spacer() + forest_plot_no_triple +
      plot_layout(widths = c(2, 0, 1.2, 0, 0.8, 0.1, 1.1))
    
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
      width = 15,
      height = 7
    )
    
  }
  
}