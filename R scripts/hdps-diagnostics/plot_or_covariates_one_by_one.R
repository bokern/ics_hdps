packages <- c("tidyverse", "ggplot2", "MetBrewer")
installed_packages <- packages %in% rownames(installed.packages())
if (any(installed_packages == FALSE)) {
  install.packages(packages[!installed_packages])
}

invisible(lapply(packages, function(pkg) {
  suppressMessages(library(pkg, character.only = TRUE, verbose = FALSE))
}))

set.seed(123)

palette <- met.brewer("Cassatt2")
setwd(HDPS_folder)
#read in csv dataset:
results_file_path <- file.path(paste0("outputs/Checks/HDPS_Cox_Results_covid_hes_present_no_triple_check.csv"))
results_file <- read.csv(results_file_path)

# Function to extract numbers from the string
extract_numbers <- function(x) {
  numbers <- str_extract_all(x, "\\d+\\.\\d+")[[1]]
  return(as.numeric(numbers))
}

# Apply the function to extract point estimate and CI bounds
results_file$OR <- sapply(results_file$Results.exp.coef...confint., function(x) extract_numbers(x)[1])
results_file$lower_ci <- sapply(results_file$Results.exp.coef...confint., function(x) extract_numbers(x)[2])
results_file$upper_ci <- sapply(results_file$Results.exp.coef...confint., function(x) extract_numbers(x)[3])

#remove the first row of the dataset
results_file <- results_file[-1,]

#rename the first entry of Covariates to 0
results_file$Covariates[1] <- 0
#make the Covariates column numeric
results_file$Covariates <- as.numeric(results_file$Covariates)

p <- ggplot(results_file, aes(x = Covariates)) +
  geom_line(aes(y = OR, group = 1), color = palette[8], size = 1.2) +
  geom_point(aes(y = OR), size = 0.8, color = palette[8]) +
  geom_line(aes(y = lower_ci, group = 1), color = palette[7], size = 1.2) +
  geom_point(aes(y = lower_ci), size = 0.8, color = palette[7]) +
  geom_line(aes(y = upper_ci, group = 1), color = palette[7], size = 1.2) +
  geom_point(aes(y = upper_ci), size = 0.8, color = palette[7]) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "gray50") +
  geom_ribbon(aes(ymin = lower_ci, ymax = upper_ci), fill = palette[6], alpha = 0.2) +  # Shaded area between lower and upper CI
  #log y axis
  scale_y_continuous(limits = c(0.5, 2.2),trans = "log10", breaks = seq(0.5, 2, 0.5)) +
  labs(x = "Number of Covariates",
       y = "Odds Ratio") +
  theme_minimal() +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major.y = element_blank(),
        axis.text = element_text(size = 10),
        axis.title = element_text(size = 12, face = "bold"),
        plot.title = element_text(size = 14, face = "bold", hjust = 0.5))


file_path <- file.path(
  HDPS_folder,
  "outputs",
  "diagnostics",
  paste0("Results_covariates_250_check.png")
)

ggsave(file_path, p, width = 10, height = 5)
