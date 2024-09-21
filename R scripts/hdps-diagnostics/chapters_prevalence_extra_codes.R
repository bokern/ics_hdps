###########################################################
# R script:    001_extraCodesVisualization.R
#
# Author:      [Your Name] (modified from John Tazare's script)
#
# Date:        [Current Date]
#
# Description: Plot for visualizing extra codes included in each
#              iteration of the HDPS model.
###########################################################

# Load and install required packages
packages <- c("tidyverse", "arrow", "MetBrewer", "gridExtra")
installed_packages <- packages %in% rownames(installed.packages())
if (any(installed_packages == FALSE)) {
  install.packages(packages[!installed_packages])
}
invisible(lapply(packages, function(pkg) {
  suppressMessages(library(pkg, character.only = TRUE, verbose = FALSE))
}))


# Set topVars as a vector
topVars <- c("100", "250", "500", "750", "1000")

outcomes <- c("covid_hes_present", "covid_death_present")
output_exts <- c("", "_no_triple")

for (output_ext in output_exts) {
  for (outcome in outcomes) {
    for (i in 2:length(topVars)) {

      prev_k <- topVars[i-1]
      current_k <- topVars[i]
      
      print(paste0("Processing ", prev_k, " to ", current_k, " ", outcome, " ", output_ext))
      
      # Load data from outputs folder for both previous and current iterations
      prev_data <- read.csv(paste0(
        HDPS_folder,
        "/outputs/top_",
        prev_k,
        "_HDPS_",
        outcome,
        output_ext,
        ".csv"
      ))
      
      current_data <- read.csv(paste0(
        HDPS_folder,
        "/outputs/top_",
        current_k,
        "_HDPS_",
        outcome,
        output_ext,
        ".csv"
      ))
      
      # Function to assign chapters to ICD and BNF codes
      assign_chapter <- function(variable) {
        # For ICD codes
        if (grepl("^[A-Z]", variable)) {
          if (substr(variable, 1, 1) == "D" &&
              as.numeric(substr(variable, 2, 3)) <= 49)
            return("D48")  # Neoplasms
          if (substr(variable, 1, 1) == "D" &&
              as.numeric(substr(variable, 2, 3)) >= 50)
            return("D89")  # Diseases of the blood...
          if (substr(variable, 1, 1) == "H" &&
              as.numeric(substr(variable, 2, 3)) <= 59)
            return("H59")  # Diseases of the eye and adnexa
          if (substr(variable, 1, 1) == "H" &&
              as.numeric(substr(variable, 2, 3)) >= 60)
            return("H95")  # Diseases of the ear and mastoid process
          return(substr(variable, 1, 1))  # For other ICD codes, return the first character
        }
        
        # For BNF codes
        if (grepl("^[0-9]", variable)) {
          return(substr(variable, 1, 2))  # Return first two digits for BNF codes
        }
        
        # If code doesn't match expected patterns, return the first character
        return(substr(variable, 1, 1))
      }
      
      
      # Filter current_data to only include new codes
      new_codes <- current_data %>%
        filter(!(variable %in% prev_data$variable)) %>%
        arrange(rank) %>%
        head(as.numeric(current_k) - as.numeric(prev_k))
      
      # Process new_codes data
      new_codes$dim <- as.factor(as.numeric(substr(new_codes$variable, 2, 2)))
      new_codes$chapter <- sapply(substr(new_codes$variable, 4, nchar(new_codes$variable)), assign_chapter)
      
      # Create summary dataset
      summary_df <- new_codes %>%
        group_by(dim, chapter) %>%
        summarise(tot = n(), .groups = 'drop')
      
      #if there are more than 50 bars, only keep the 50 with the highest tot
      if(nrow(summary_df) > 30){
        summary_df <- summary_df %>%
          arrange(desc(tot)) %>%
          slice(1:30)
      }
      
      # Add chapter indicators and descriptions
      chapter_ind <- c(
        "A" = 1,
        "B" = 1,
        "C" = 2,
        "D48" = 2,
        "D89" = 3,
        "E" = 4,
        "F" = 5,
        "G" = 6,
        "H59" = 7,
        "H95" = 8,
        "I" = 9,
        "J" = 10,
        "K" = 11,
        "L" = 12,
        "M" = 13,
        "N" = 14,
        "O" = 15,
        "P" = 16,
        "Q" = 17,
        "R" = 18,
        "S" = 19,
        "T" = 19,
        "U" = 20,
        "V" = 21,
        "W" = 21,
        "X" = 21,
        "Y" = 21,
        "Z" = 22,
        "01" = 23,
        "02" = 24,
        "03" = 25,
        "04" = 26,
        "05" = 27,
        "06" = 28,
        "07" = 29,
        "08" = 30,
        "09" = 31,
        "10" = 32,
        "11" = 33,
        "12" = 34,
        "13" = 35,
        "14" = 36,
        "15" = 37,
        "18" = 38,
        "19" = 39,
        "20" = 40,
        "21" = 41,
        "22" = 42,
        "23" = 43
      )
      
      chapter_descriptions <- c(
        "1" = "Certain infectious and parasitic diseases",
        "2" = "Neoplasms",
        "3" = "Blood and blood-forming organs diseases",
        "4" = "Endocrine, nutritional and metabolic diseases",
        "5" = "Mental and behavioural disorders",
        "6" = "Diseases of the nervous system",
        "7" = "Diseases of the eye and adnexa",
        "8" = "Diseases of the ear and mastoid process",
        "9" = "Diseases of the circulatory system",
        "10" = "Diseases of the respiratory system",
        "11" = "Diseases of the digestive system",
        "12" = "Diseases of the skin and subcutaneous tissue",
        "13" = "Diseases of the musculoskeletal system and connective tissue",
        "14" = "Diseases of the genitourinary system",
        "15" = "Pregnancy, childbirth and the puerperium",
        "16" = "Certain conditions originating in the perinatal period",
        "17" = "Congenital malformations, deformations and chromosomal abnormalities",
        "18" = "Symptoms, signs and abnormal findings",
        "19" = "Injury, poisoning and other external causes",
        "20" = "Codes for special purposes",
        "21" = "Health status factors and contact with health services",
        "22" = "Codes for special purposes",
        "23" = "Gastrointestinal system",
        "24" = "Cardiovascular system",
        "25" = "Respiratory system",
        "26" = "Central nervous system",
        "27" = "Infections",
        "28" = "Endocrine system",
        "29" = "Obstetrics, gynaecology and urinary tract disorders",
        "30" = "Malignant disease and immunosuppression",
        "31" = "Nutrition and blood",
        "32" = "Musculoskeletal and joint diseases",
        "33" = "Eye",
        "34" = "Ear, nose and oropharynx",
        "35" = "Skin",
        "36" = "Immunological products and vaccines",
        "37" = "Anaesthesia",
        "38" = "Preparations used in diagnosis",
        "39" = "Other drugs and preparations",
        "40" = "Dressings",
        "41" = "Appliances",
        "42" = "Incontinence appliances",
        "43" = "Stoma appliances"
      )
      
      summary_df$chapter_ind <- chapter_ind[summary_df$chapter]
      summary_df$description <- chapter_descriptions[as.character(summary_df$chapter_ind)]
      
      #remove chapter
      summary_df <- summary_df %>% dplyr::select(-chapter)
      # combine rows with the same chapter indicator, sum the total and keep the description
      summary_df <- summary_df %>%
        group_by(dim, chapter_ind, description) %>%
        summarise(tot = sum(tot), .groups = 'drop')
      
      # Create whitespace to separate each dimension by adding empty bars
      empty_bar <- 1.5
      to_add <- data.frame(
        dim = rep(levels(summary_df$dim), each = empty_bar),
        chapter_ind = NA,
        description = NA,
        tot = 0
      )
      
      # Combine the original data with the empty bars
      summary_df <- rbind(summary_df, to_add)
      summary_df <- summary_df %>% arrange(dim, tot)
      
      # Add id column for plotting
      summary_df$id <- seq(1, nrow(summary_df))
      
      # Get the name, angles and position of dimension chapter
      number_of_bar <- nrow(summary_df)
      angle <- 90 - 360 * (summary_df$id - 0.5) / number_of_bar
      summary_df$hjust <- ifelse(angle < -90, 1, 0)
      summary_df$angle <- ifelse(angle < -90, angle + 180, angle)
      summary_df$prev <- (summary_df$tot /as.numeric(k)) *100
      
      max_prev <- max(summary_df$prev, na.rm = TRUE)
      summary_df <- summary_df %>%
        #        group_by(dim) %>%
        mutate(normalized_prev = tot / sum(tot) * 100) #%>%
      #        ungroup()
      
      # Set a fixed maximum radius for the plot
      fixed_max_radius <- 8
      
      # Scale the normalized prevalence to fit within the fixed radius
      summary_df <- summary_df %>%
        mutate(scaled_prev = (normalized_prev*0.5) / max((normalized_prev), na.rm = TRUE) * fixed_max_radius)
      
      # Calculate label positions
      label_data <- summary_df %>%
        filter(!is.na(description)) %>%
        mutate(
          angle = 90 - 360 * (id - 0.5) / number_of_bar,
          hjust = ifelse(angle < -90, 1, 0),
          angle = ifelse(angle < -90, angle + 180, angle),
          label_y = scaled_prev * 1.05,  # Position labels just outside the fixed radius
          label_x = id
        )
      
      # Calculate percentage label positions
      label_percentages <- summary_df %>%
        mutate(
          angle = 90 - 360 * (id - 0.5) / number_of_bar,
          hjust = ifelse(angle < -90, 1, 0),
          angle = ifelse(angle < -90, angle + 180, angle),
          label_y = pmax(scaled_prev * 0.5, 0.5),  # Position percentages inside bars, with a minimum height
          label_x = id
        ) %>% 
        filter(normalized_prev != 0)
      
      # Modify the legend_data to include a legend title and labels
      legend_data <- data.frame(
        dim = c("1", "2", "3"),
        label = c("Clinical", "Prescriptions", "Hospitalisations"),
        x = 0,
        y = c(-5, -6, -7)  # Adjusted y values to move the legend downwards
      )
      
      # Define the colour palette outside the ggplot call
      colours <- c("1" = palette[1], "2" = palette[3], "3" = palette[5])
      # Modify the plotting code to reflect the changes in data
      # Update the title to show the current comparison
      plot_title <- paste("Extra Codes from", prev_k, "to", current_k, "Variables")
      
      # Update the ggplot code (example, adjust as needed)
      p <- ggplot(summary_df, aes(x = as.factor(id), y = scaled_prev, fill = dim)) +
        scale_fill_manual(values = colours,
                          labels = c("1" = "Clinical", "2" = "Prescriptions", "3" = "Hospital"),
                          guide = guide_legend(keywidth = 2, keyheight = 2, size = 24)) +
        geom_bar(stat = "identity", alpha = 1, width = 0.85) +
        ylim(-fixed_max_radius * 0.5, fixed_max_radius * 1.2) +  # Adjust the upper limit to accommodate labels
        theme_minimal() +
        theme(
          legend.position = c(0.51, 0.48),
          legend.justification = c(0.5, 0.2),
          legend.title = element_blank(),
          legend.text = element_text(size = 14),
          axis.text = element_blank(),
          axis.title = element_blank(),
          panel.grid = element_blank(),
          plot.margin = unit(rep(-1, 4), "cm")  # Adjust as needed
        ) +
        coord_polar(start = 0) +
        geom_text(
          data = label_data,
          aes(x = label_x, y = label_y, label = description, hjust = hjust),
          color = "black",
          alpha = 0.8,
          size = 4.5,  # Reduced size for better fit
          angle = label_data$angle,
          inherit.aes = FALSE
        ) +
        geom_text(
          data = label_percentages,
          aes(x = label_x, y = -fixed_max_radius *0.12 , label = sprintf("%.1f%%", prev), hjust = hjust),
          color = "black", 
          alpha = 0.8,
          size = 4.5,  # Reduced size for better fit
          angle = label_percentages$angle,
          inherit.aes = FALSE
        )
      
      print(p)
      
      # Update the file naming convention
      file_path <- file.path(
        HDPS_folder,
        "outputs",
        "diagnostics",
        paste0("extraCodesPlot_", prev_k, "to", current_k, "_", outcome, output_ext, ".png")
      )
      ggsave(file_path, p, width = 16, height = 16)
    }
  }
}