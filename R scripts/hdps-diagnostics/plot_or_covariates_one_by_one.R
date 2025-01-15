packages <- c("tidyverse", "ggplot2", "MetBrewer")
installed_packages <- packages %in% rownames(installed.packages())
if (any(installed_packages == FALSE)) {
  install.packages(packages[!installed_packages])
}
invisible(lapply(packages, function(pkg) {
  suppressMessages(library(pkg, character.only = TRUE, verbose = FALSE))
}))
set.seed(123)
# Set color palette
palette <- met.brewer("Cassatt2")
setwd(HDPS_folder)
# Define the file paths and model types
model_info <- list(
  cox = list(
    path = "outputs/Checks/HDPS_Cox_Results_covid_hes_present_no_triple_check.csv",
    name = "Cox",
    y_label = "Hazard Ratio"  # Added y-axis label for Cox
  ),
  logistic = list(
    path = "outputs/Checks/HDPS_Logistic_Results_covid_hes_present_no_triple_check.csv",
    name = "Logistic",
    y_label = "Odds Ratio"    # Added y-axis label for Logistic
  )
)
# Function to extract numbers from the string
extract_numbers <- function(x) {
  numbers <- str_extract_all(x, "\\d+\\.\\d+")[[1]]
  return(as.numeric(numbers))
}
# Loop through each model type
for (model_type in names(model_info)) {
  # Read in data
  results_file_path <- file.path(model_info[[model_type]]$path)
  results_file <- read.csv(results_file_path)
  
  # Process the results
  results_file$OR <- sapply(results_file$Results.exp.coef...confint., function(x) extract_numbers(x)[1])
  results_file$lower_ci <- sapply(results_file$Results.exp.coef...confint., function(x) extract_numbers(x)[2])
  results_file$upper_ci <- sapply(results_file$Results.exp.coef...confint., function(x) extract_numbers(x)[3])
  
  # Remove first row and fix covariates
  results_file <- results_file[-1,]
  results_file$Covariates[1] <- 0
  results_file$Covariates <- as.numeric(results_file$Covariates)
  
  # Create plot
  p <- ggplot(results_file, aes(x = Covariates)) +
    geom_line(aes(y = OR, group = 1), color = palette[8], size = 1.2) +
    geom_point(aes(y = OR), size = 0.8, color = palette[8]) +
    geom_line(aes(y = lower_ci, group = 1), color = palette[7], size = 1.2) +
    geom_point(aes(y = lower_ci), size = 0.8, color = palette[7]) +
    geom_line(aes(y = upper_ci, group = 1), color = palette[7], size = 1.2) +
    geom_point(aes(y = upper_ci), size = 0.8, color = palette[7]) +
    geom_hline(yintercept = 1, linetype = "dashed", color = "gray50") +
    geom_ribbon(aes(ymin = lower_ci, ymax = upper_ci), fill = palette[6], alpha = 0.2) +
    scale_y_continuous(limits = c(0.5, 2.2), trans = "log10", breaks = seq(0.5, 2, 0.5)) +
    labs(x = "Number of Covariates",
         y = model_info[[model_type]]$y_label) +  # Use model-specific y-axis label
    theme_minimal() +
    theme(panel.grid.minor = element_blank(),
          panel.grid.major.y = element_blank(),
          axis.text = element_text(size = 10),
          axis.title = element_text(size = 12, face = "bold"),
          plot.title = element_text(size = 14, face = "bold", hjust = 0.5))
  
  # Save the plot
  file_path <- file.path(
    HDPS_folder,
    "outputs",
    "diagnostics",
    paste0("Results_covariates_250_check_", tolower(model_info[[model_type]]$name), ".png")
  )
  ggsave(file_path, p, width = 10, height = 5)
}