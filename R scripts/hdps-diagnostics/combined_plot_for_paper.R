packages <- c("tidyverse", "arrow", "MetBrewer", "ggplot2", "patchwork", "cowplot", "ggpubr")
installed_packages <- packages %in% rownames(installed.packages())
if (any(installed_packages == FALSE)) {
  install.packages(packages[!installed_packages])
}
invisible(lapply(packages, function(pkg) {
  suppressMessages(library(pkg, character.only = TRUE, verbose = FALSE))
}))

palette <- met.brewer("Cassatt2")


k <- "250"
outcome <- "covid_hes_present"
cohort_ext <- ""

setwd(HDPS_folder)
#go to subfolder outputs
setwd("outputs")

#switch statements to create labels for outcome and cohort_ext
outcome_label <- switch(
  outcome,
  covid_hes_present = "COVID-19 Hospitalisation",
  covid_death_present = "COVID-19 Death"
)

if (cohort_ext == "") {
  cohort_label <- " including triple therapy users"
} else if (cohort_ext == "_no_triple") {
  cohort_label <- " excluding triple therapy users"
} else {
  cohort_label <- ""  # default case if needed
}

#  browser()
#list all rds files with 250 covariates, covid_hes_present 
rds_files <- list.files(pattern = paste0(k, "_.*", outcome, ".*\\.rds$"))

#if cohort_ext is "_no_triple", keep rds_files with notriple in the name
if (cohort_ext == "_no_triple") {
  rds_files <- rds_files[grep("_no_triple", rds_files)]
}
#if cohort_ext is "", remove rds_files with notriple in the name
if (cohort_ext == "") {
  rds_files <- rds_files[!grepl("_no_triple", rds_files)]
}

#add predefinedUnweighted and PredefinedWeighted plots to the list of rds_files
rds_files <- c(rds_files, paste0("Predefined_overlapUnweighted_", outcome, cohort_ext, ".rds"))
rds_files <- c(rds_files, paste0("Predefined_overlapWeighted_", outcome, cohort_ext, "_trimmed.rds"))
# Initialize a list to store the plots by filename
plot_list <- list()

# Loop over each file, load the plot, save it, and store it in the environment
for (i in seq_along(rds_files)) {
  plot <- readRDS(rds_files[i])
  
  if (grepl("Weighted", rds_files[i])) {
    plot <- plot + 
      guides(color = guide_legend(title = NULL))  # Remove the legend title
  }
  
  # Save plot as PNG
  ggsave(
    gsub(".rds", "_reimported.png", rds_files[i]),
    plot = plot,
    width = 3,
    height = 2,
    dpi = 300
  )
  
  # Store plot in the list with a name for reference
  plot_list[[rds_files[i]]] <- plot
  
  # Assign each plot to a named variable in the global environment
  assign(gsub(".rds", "_plot", rds_files[i]), plot, envir = .GlobalEnv)
  
  # Garbage collection to manage memory
  gc()
  rm(plot)
}

# Remove plot and plot_list from the global environment
rm(plot_list)

# Arrange the plots using ggarrange
unweighted_plot_name <- paste0("covariates_HDPS_", k, "_overlapUnweighted_", outcome, cohort_ext, "_plot")
weighted_plot_name <- paste0("covariates_HDPS_", k, "_overlapWeighted_", outcome, cohort_ext, "_trimmed_plot")
unweighted_predefined_plot_name <- paste0("Predefined_overlapUnweighted_", outcome, cohort_ext, "_plot")
weighted_predefined_plot_name <- paste0("Predefined_overlapWeighted_", outcome, cohort_ext, "_trimmed_plot")

unweighted_plot <- get(unweighted_plot_name)
weighted_plot <- get(weighted_plot_name)
weighted_predefined_plot <- get(weighted_predefined_plot_name)
unweighted_predefined_plot <- get(unweighted_predefined_plot_name)

base_theme <- theme_minimal() +
  theme(
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 14),
    legend.text = element_text(size = 12),
    plot.title = element_blank(),
    plot.margin = margin(20, 20, 20, 20),
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA),
    panel.grid = element_blank(),
    axis.line = element_line(color = "black"),
    panel.border = element_blank())


unweighted_plot <- get(unweighted_plot_name) +
  base_theme +
  theme(legend.position = "none") +
  scale_x_continuous(limits = c(0, 1))

weighted_plot <- get(weighted_plot_name) +
  base_theme +
  theme(
    legend.position = "none",
    axis.title.y = element_blank()
  ) +
  scale_color_manual(values = c("ICS" = palette[9], "LABA/LAMA" = palette[4])) +
  scale_fill_manual(values = c("ICS" = palette[9], "LABA/LAMA" = palette[4])) +
  scale_x_continuous(limits = c(0, 1)) + 
  #change y axis breaks to whole numbers
  scale_y_continuous(breaks = seq(0, 3, 1))

unweighted_predefined_plot <- get(unweighted_predefined_plot_name) +
  base_theme +
  theme(
    legend.position = "none",
    axis.title.x = element_blank()
  ) +
  scale_x_continuous(limits = c(0, 1))

weighted_predefined_plot <- get(weighted_predefined_plot_name) +
  base_theme +
  theme(
    legend.position = "none",
    axis.title.y = element_blank(),
    axis.title.x = element_blank()
  ) +
  scale_color_manual(values = c("ICS" = palette[9], "LABA/LAMA" = palette[4])) +
  scale_fill_manual(values = c("ICS" = palette[9], "LABA/LAMA" = palette[4])) +
  scale_x_continuous(limits = c(0, 1))

legend_plot <- get(weighted_plot_name) +
  scale_color_manual(values = c("ICS" = palette[9], "LABA/LAMA" = palette[4])) +
  scale_fill_manual(values = c("ICS" = palette[9], "LABA/LAMA" = palette[4])) +
  theme(
    legend.text = element_text(size = 12),
    legend.box = "horizontal",
    legend.spacing.x = unit(0.6, "cm"),
    legend.key.size = unit(0.6, "cm")
  ) +
  guides(
    color = guide_legend(
      override.aes = list(
        fill = c(palette[9], palette[4]),
        size = 5,
        shape = 15
      ), 
      title = NULL,
      nrow = 1,
      keywidth = unit(0.8, "cm"),
      keyheight = unit(0.8, "cm")
    )
  )

shared_legend <- get_legend(legend_plot)

# Arrange the top four plots in a 2x2 grid with labels
top_plots <- ggarrange(
  unweighted_predefined_plot, weighted_predefined_plot,
  unweighted_plot, weighted_plot,
  ncol = 2, nrow = 2,
  widths = c(1, 1),
  heights = c(1, 1),
  labels = c("A", "B", "C", "D"),
  font.label = list(size = 16, face = "bold")
)

# Combine all elements
combined_plot <- ggarrange(
  top_plots,
  as_ggplot(shared_legend),
  ncol = 1,
  heights = c(1.5, 0.15)
)

# Save the plot with a white background
ggsave(
  file.path(HDPS_folder, "outputs", "diagnostics", 
            paste0("combined_ps_plot_", k, "_", outcome, cohort_ext, "_for_paper.png")),
  combined_plot,
  width = 8,
  height = 6,
  dpi = 600,
  bg = "white"
)


