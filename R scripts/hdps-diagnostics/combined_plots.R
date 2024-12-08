packages <- c("tidyverse", "arrow", "MetBrewer", "ggplot2", "ggpubr")
installed_packages <- packages %in% rownames(installed.packages())
if (any(installed_packages == FALSE)) {
  install.packages(packages[!installed_packages])
}
invisible(lapply(packages, function(pkg) {
  suppressMessages(library(pkg, character.only = TRUE, verbose = FALSE))
}))

palette <- met.brewer("Cassatt2")

topVars <- c("100", "250", "500", "750")
outcomes <- c("covid_hes_present", "covid_death_present")
cohort_exts <- c("", "_no_triple")

k <- "1000"
outcome <- "covid_death_present"
cohort_ext <- ""

for (outcome in outcomes){
  print(paste0("Processing ", outcome))
  for (cohort_ext in cohort_exts){
    print(paste0("Processing ", outcome, cohort_ext))
    for (k in topVars) {
      print(paste0("Processing HDPS_", k, "_SMD_", outcome, cohort_ext ))
      
      gc()
 
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
        
        # Assign each plot to a named variable in the global environment
        assign(gsub(".rds", "_plot", rds_files[i]), plot, envir = .GlobalEnv)
        
        # Garbage collection to manage memory
        gc()
       
      }
      
      rm(plot)

      # Prepare the plot names
      unweighted_plot_name <- paste0("covariates_HDPS_", k, "_overlapUnweighted_", outcome, cohort_ext, "_plot")
      weighted_plot_name <- paste0("covariates_HDPS_", k, "_overlapWeighted_", outcome, cohort_ext, "_trimmed_plot")
    
      # Get the plots from their dynamic names using get()
      unweighted_plot <- get(unweighted_plot_name)
      weighted_plot <- get(weighted_plot_name)

      #import the chapters_prevalence data set for 250 covariates, covid_hes_present 
      setwd(HDPS_folder)
      setwd("outputs/diagnostics")
      # Read the PNG file and store it in the global environment
      smd_path <- file.path(
        HDPS_folder,
        "outputs",
        "diagnostics",
        paste0("SMD_", k, outcome, cohort_ext, ".rds"))
      
      smd_plot <- readRDS(smd_path)
      
      # Get the x-axis limits directly from the data
      smd_plot_xmin <- -18
      smd_plot_xmax <- max(smd_plot$data$adjusted_rank)
      
      calculate_icon_size <- function(k) {
        if (k <= 100) {
          return(2)
        } else if (k <= 500) {
          return(1.7)
        } else if (k <= 750) {
          return(1.3)
        } else {
          return(0.1)
        }
      }
      
      icon_size <- calculate_icon_size(k)
      
      # Ensure the x-axis limits are consistent across all plots
      smd_plot <- smd_plot +
        scale_x_continuous(limits = c(min(smd_plot_xmin, unweighted_plot$scales$scales$x$limits[1]),
                                      max(smd_plot_xmax, unweighted_plot$scales$scales$x$limits[2])),
                           expand = expansion(mult = 0.005)) +
        guides(color = guide_legend(override.aes = list(size = 6))) +
        #double the icon size for all plots
        geom_point(aes(x = adjusted_rank, y = smd_unweighted, color = "Unweighted"), 
                   shape = 16, alpha = 1, size = icon_size*2) +
        # Only plot weighted SMD where not NA
        geom_point(data = . %>% filter(!is.na(smd_weighted)),
                   aes(x = adjusted_rank, y = smd_weighted, color = "HDPS-Weighted"), 
                   shape = 18, alpha = 1, size = icon_size*2) +
        geom_point(data = . %>% filter(!is.na(smd_weighted_predefined)),
                   aes(x = adjusted_rank, y = smd_weighted_predefined, color = "Predefined Covariates"), 
                   shape = 17, alpha = 1, size = icon_size*2) +
        theme(legend.text = element_text(size = 20),
              legend.title = element_text(size = 24),
              legend.key.size = unit(1.1, "cm"),
              plot.margin = margin(t = 1, b = 0, l = 1, r = 1, "cm"),
              axis.text = element_text(size = 18),
              axis.title = element_text(size = 20))
   
      gc()
      
      # First ensure the top plots have consistent sizing
      unweighted_plot <- unweighted_plot + 
        geom_density(size = 1.5) +
        scale_y_continuous(limits = c(0, 4), breaks = seq(0, 4, 1)) +
        theme(
          legend.position = "none", 
          plot.background = element_rect(fill = "white", color = NA),
          panel.background = element_rect(fill = "white", color = NA),
          plot.title = element_blank(),
          aspect.ratio = 0.6,
          axis.text = element_text(size = 18),
          axis.title = element_text(size = 20),
          line = element_line(linewidth = 1.8),
          plot.margin = margin(t = 0, r = 0, b = -10, l = 0)
        )
      
      weighted_plot <- weighted_plot +
        geom_density(size = 1.5) +
        scale_color_manual(values = c("ICS" = palette[9], "LABA/LAMA" = palette[4])) +
        scale_fill_manual(values = c("ICS" = palette[9], "LABA/LAMA" = palette[4])) +
        scale_x_continuous(limits = c(0, 1)) +
        scale_y_continuous(limits = c(0, 4), breaks = seq(0, 4, 1)) +
        theme(
          legend.position = "none",
          plot.background = element_rect(fill = "white", color = NA),
          panel.background = element_rect(fill = "white", color = NA),
          aspect.ratio = 0.6,
          plot.title = element_blank(),
          axis.text = element_text(size = 18),
          axis.title = element_text(size = 20),
          axis.title.y = element_blank(),
          line = element_line(linewidth = 1.8),
          plot.margin = margin(t = 0, r = 0, b = -10, l = 0)
        )
      
      legend_plot <- get(weighted_plot_name) +
        scale_color_manual(values = c("ICS" = palette[9], "LABA/LAMA" = palette[4])) +
        scale_fill_manual(values = c("ICS" = palette[9], "LABA/LAMA" = palette[4])) +
        theme(
          legend.text = element_text(size = 18),
          legend.box = "horizontal",
          legend.spacing.x = unit(0.6, "cm"),
          legend.key.size = unit(0.6, "cm"),
          plot.margin = margin(t = -20, r = 0, b = 0, l = 0)
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
            keywidth = unit(1, "cm"),
            keyheight = unit(1, "cm")
          )
        )
      
      shared_legend <- get_legend(legend_plot)
      
      # First combine the top two plots
      top_plots <- ggarrange(
        unweighted_plot, weighted_plot,
        ncol = 2,
        widths = c(1, 1),
        labels = c("A", "B"),
        font.label = list(size = 24, face = "bold", vjust = 1, hjust = 0.5),
        align = "hv"
      )
      
      # Combine top plots with legend
      plots_with_legend <- ggarrange(
        shared_legend,
        top_plots,
        nrow = 2,
        heights = c(0.2, 1)  # Adjust this ratio to control legend height
      )
      
      # Finally combine everything with the smd plot
      combined_plot <- ggarrange(
        plots_with_legend,
        smd_plot,
        nrow = 2,
        heights = c(1.1, 1.5),  # Adjust these ratios as needed
        labels = c("", "C"),
        font.label = list(size = 24, face = "bold", vjust = -200, hjust = 0.5)
      ) +
        theme(
          plot.margin = margin(t = 20, b = 0, l = 0, r = 0),
          plot.background = element_rect(fill = "white", color = NA),
          panel.background = element_rect(fill = "white", color = NA)
        )
      
      # Add title
      combined_plot_with_title <- annotate_figure(
        combined_plot,
        top = text_grob(
          paste0("HDPS Analysis: Top ", k, " Covariates for ", outcome_label, cohort_label), 
          face = "bold", 
          size = 18,
          hjust = 0.5,
          #move title down 
          vjust = 2
        )
      )
      
      # Save with adjusted dimensions
      ggsave(
        file.path(HDPS_folder, "outputs", "diagnostics",
                  paste0("combined_HDPS_plot_", k, "_", outcome, cohort_ext, ".png")),
        combined_plot_with_title,
        width = 12,
        height = 15,
        dpi = 600,
        bg = "white"
      )      
      
      #remove all the plots with covariates_HDPS in the name or HDPS_weight_dist
      rm(list = ls(pattern = "covariates_HDPS|HDPS_weight_dist|plot"))
      gc()
      
    }
    gc()
  }
  gc()
}

