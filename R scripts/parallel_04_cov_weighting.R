# -----------------------------------------------------------------------------
# PROGRAM NAME:  502_HDPS_ACNU_Analysis
# PROJECT:      xzx26005
# AUTHOR:        John Tazare
# DATE CREATED:   10 Sep 2020
# NOTES:        high-dimensional Propensity Score Analysis for ACNU:
#                                  1. Summarise baseline covariates
#                                  2. Estimate propensity score
#                                  3. Propensity score matching and evaluation
#                                  4. Estimate the treatment effect

# -----------------------------------------------------------------------------
# Flag to enable debugging; if TRUE, script will only be run on 1,000 patients
# selected at random
debugMode <- F


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
  "MetBrewer",
  "furrr",
  "future"
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

# Assign priority of select to dplyr
select <- dplyr::select

topVars <- c("10", "100", "250", "500", "750", "1000")
outcomes <- c("covid_hes_present", "covid_death_present")
cohort_exts <- c("", "_no_triple")

cov_weighting <- function(cohort_ext, outcome, topVar) {
  
  for (cohort_ext in cohort_exts) {
    # Iterate through each outcome
    for (outcome in outcomes) {
      
      print(paste("Analyzing outcome:", outcome))
      
      timeinstudy <- switch(outcome,
                            "covid_hes_present" = "timeinstudy2",
                            "covid_death_present" = "timeinstudy3")
      
      # Load in patient summary
      pat_summary_data <- read_parquet("copd_wave1_60d.parquet") %>%
        mutate_at(c("patid"), as.character) %>%
        #remove if treatgroup is NA
        filter(!is.na(treatgroup))
      
      if (cohort_ext == "_no_triple") {
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
      
      hdpsPredefinedVars <- pat_summary_data %>% dplyr::select(-all_of(exclude_vars))
      
      # Remove treatgroup from the dataframe
      hdpsCovariates <- hdpsPredefinedVars %>% dplyr::select(-treatgroup)
      
      follow_up_times <- pat_summary_data %>% dplyr::select("patid", timeinstudy)
      
      hdpsCohort <- pat_summary_data %>% dplyr::select(all_of(c("patid", "treatgroup")), all_of(outcome))
      
      #Create pscores and a model based on the predefined variables only
      
      # Predefined covariates modeling
      predefined_vars <- setdiff(names(hdpsCovariates), "patid")
      
      psModelFunction_predefined <- as.formula(paste("treatgroup ~", paste(predefined_vars, collapse = " + "), sep = " "))
      
      psModel_predefined <- glm(psModelFunction_predefined,
                                family  = binomial(link = "logit"),
                                data    = hdpsPredefinedVars)
      
      hdpsPredefinedVars$pscore_predefined <- predict(psModel_predefined, type = "response")
      
      # browser()
      #save min and max propensity scores in each treatment group, save limits of area of common support
      ps_trim <- hdpsPredefinedVars %>% 
        dplyr::select(treatgroup, pscore_predefined) %>% 
        group_by(treatgroup) %>% 
        summarise(min = min(pscore_predefined), max = max(pscore_predefined)) %>% 
        ungroup() %>% 
        summarise(min = max(min), max = min(max))
      
      #exclude people outside of common support
      hdpsPredefinedVars <- hdpsPredefinedVars %>% 
        filter(pscore_predefined >= ps_trim$min & pscore_predefined <= ps_trim$max)
      
      # IPTW for predefined covariates
      hdpsPredefinedVars <- hdpsPredefinedVars %>%
        mutate(iptw_weight_predefined = case_when(
          treatgroup == 1 ~ 1 / pscore_predefined,
          treatgroup == 0 ~ 1 / (1 - pscore_predefined)
        ))
      
      # Merge follow-up times into hdpsPredefinedVars
      hdpsPredefinedVars[[timeinstudy]] <- follow_up_times[[timeinstudy]][match(hdpsPredefinedVars$patid, follow_up_times$patid)]
      
      # Save the hdpsPredefinedVars to parquet
      file_path_predefined <- file.path(
        HDPS_folder,
        "outputs",
        paste0(
          ifelse(debugMode, "debug_", ""),
          "HDPS_predefined",
          outcome,
          cohort_ext,
          "_trimmed.parquet"
        )
      )
      
      write_parquet(hdpsPredefinedVars, file_path_predefined)
      
      # Combine predefined and HDPS covariates if hdpsData exists
      all_vars_table <- setdiff(
        names(hdpsPredefinedVars),
        c("patid", "treatgroup", "iptw_weight_predefined", "pscore_predefined", outcome))
      
      # Vector of categorical variables that need transformation
      catVars <- all_vars_table[!all_vars_table %in% c("age_index", "pscore_predefined")]
      
      
      # Create survey design for predefined variables
      svy_design_predefined <- svydesign(
        ids = ~ 1,
        weights = ~ iptw_weight_predefined,
        data = hdpsPredefinedVars)
      
      tabWeighted_predefined <- svyCreateTableOne(
        vars = all_vars_table,
        strata = "treatgroup",
        data = svy_design_predefined,
        test = FALSE,
        factorVars = catVars)
      
      # Export TableOne for predefined covariates
      write.csv(print(
        tabWeighted_predefined,
        smd = TRUE,
        quote = FALSE,
        noSpaces = TRUE,
        printToggle = FALSE),
        file = paste0(
          HDPS_folder,
          "/outputs/Predefined_table1Weighted_",
          outcome,
          cohort_ext,
          "_trimmed.csv"
        )
      )
      
      # Plot overlap with IPTW for predefined covariates
      iptw_plot_predefined <- hdpsPredefinedVars %>%
        mutate(trtlabel = ifelse(treatgroup == 1, 'ICS', 'LABA/LAMA')) %>%
        ggplot(aes(x = pscore_predefined, color = trtlabel)) +
        geom_density(alpha = 0.5, linewidth = 1) +
        labs(
          x = 'Probability of receiving ICS (Predefined)',
          y = 'Density',
          title = 'Overlap plot with IPTW (Predefined)'
        ) +
        scale_color_manual(values = c("ICS" = palette[9], "LABA/LAMA" = palette[4])) +
        scale_x_continuous(breaks = seq(0, 1, 0.1)) +
        theme_minimal() +
        theme(legend.title = element_blank()) +
        guides(color = guide_legend(override.aes = list(fill = c(
          palette[9], palette[4]
        ))))
      
      ggsave(
        paste0(
          HDPS_folder,
          "/outputs/Predefined_overlapWeighted_",
          outcome,
          cohort_ext,
          "_trimmed.png"),
        iptw_plot_predefined, height = 4, width = 6)
      
      for (k in topVars) {
        # browser()
        print(paste("Assessing", k, "covariates"))
        
        # Load hdpsData
        file_path_top <- file.path(
          HDPS_folder,
          "outputs",
          paste0(
            ifelse(debugMode, "debug_", ""),
            "covariates_HDPS_",
            k,
            outcome,
            cohort_ext,
            ".parquet"
          )
        )
        
        hdpsData <- read_parquet(file_path_top)
        hdpsData <- hdpsData %>% left_join(hdpsCohort, by = "patid")
        
        # Combine predefined and HDPS covariates if hdpsData exists
        all_covariates <- hdpsData %>%
          dplyr::select(-any_of(c("iptw_weight", "pscore"))) %>%  # Remove HDPS-specific columns if they exist
          left_join(
            hdpsPredefinedVars %>% dplyr::select(patid, pscore_predefined, iptw_weight_predefined),
            by = "patid")
        
        # Define all variables for the table
        all_vars_table <- setdiff(
          names(all_covariates),
          c("patid",
            "treatgroup",
            "iptw_weight_predefined",
            "pscore_predefined",
            outcome
          )
        )
        
        # Summarise data before weighting  -----------------------------------------
        
        cols_to_exclude <- c("patid", outcome, "treatgroup")
        
        #vars are all the variables used as covariates
        vars <- setdiff(names(hdpsData), cols_to_exclude)
        
        # Vector of categorical variables that need transformation
        catVars <- vars[!vars %in% c("age_index")]
        
        if (!"treatgroup" %in% names(hdpsData)) {
          # Select only 'patid' and 'treatgroup' from hdpsCohort
          hdpsCohort_treat <- hdpsCohort %>% dplyr::select(patid, treatgroup)
          # Perform the left_join using the subsetted hdpsCohort
          hdpsData <- hdpsData %>% left_join(hdpsCohort_treat, by = "patid")
        }
        
        # Create a TableOne object
        tabUnweighted <- CreateTableOne(
          vars = vars,
          strata = "treatgroup",
          data = hdpsData,
          test = FALSE,
          factorVars = catVars)
        
        # Export TableOne to CSV
        write.csv(
          print(
            tabUnweighted,
            smd = TRUE,
            quote = FALSE,
            noSpaces = TRUE,
            printToggle = FALSE
          ),
          file = paste0(
            HDPS_folder,
            "/outputs/HDPS_main_",
            k,
            "_table1Unweighted_",
            outcome,
            cohort_ext,
            ".csv"
          )
        )
        
        # Fit propensity score model  ---------------------------------------------
        # Specify model
        
        psModelFunction <- as.formula(paste("treatgroup ~", paste(vars, collapse = " + "), sep = " "))
        
        psModel <- glm(psModelFunction,
                       family  = binomial(link = "logit"),
                       data    = hdpsData)
        
        hdpsData$pscore <- predict(psModel, type = "response")
        
        #save min and max propensity scores in each treatment group, save limits of area of common support
        ps_trim <- hdpsData %>% 
          dplyr::select(treatgroup, pscore) %>% 
          group_by(treatgroup) %>% 
          summarise(min = min(pscore), max = max(pscore)) %>% 
          ungroup() %>% 
          summarise(min = max(min), max = min(max))
        
        #exclude people outside of common support
        hdpsData <- hdpsData %>% 
          filter(pscore >= ps_trim$min & pscore <= ps_trim$max)
        
        # IPTW for predefined covariates
        hdpsData <- hdpsData %>%
          mutate(iptw_weight = case_when(treatgroup == 1 ~ 1 / pscore, treatgroup == 0 ~ 1 / (1 - pscore)))
        
        #plot the propensity scores
        plot <- ggplot(hdpsData, aes(x = pscore, color = factor(treatgroup))) +
          geom_density(linewidth = 1) +
          labs(
            x = "Propensity Score",
            y = "Density",
            title = "Propensity score distribution (unweighted)",
            color = NULL
          ) +
          scale_x_continuous(limits = c(0, 1)) +
          scale_color_manual(
            values = c(palette[4], palette[9]),
            labels = c("LABA/LAMA", "ICS/LABA")
          ) +
          theme_minimal() +
          guides(color = guide_legend(override.aes = list(fill = c(
            palette[4], palette[9]
          ))))
        
        ggsave(
          paste0(
            HDPS_folder,
            "/outputs/covariates_HDPS_",
            k,
            "_overlapUnweighted_",
            outcome,
            cohort_ext,
            ".png"
          )
        )
        
        # Propensity Score Weighting  ----------------------------------------------
        
        #IPTW
        hdpsData <- hdpsData %>%
          mutate(iptw_weight = case_when(treatgroup == 1 ~ 1 / pscore, treatgroup == 0 ~ 1 / (1 - pscore)))
        
        #calculate weighted propensity score
        hdpsData <- hdpsData %>%
          mutate(weighted_ps = iptw_weight * pscore) %>% 
          dplyr::select(patid, pscore, iptw_weight, weighted_ps, everything())
        
        # plot the weights by treatment group
        plot <- hdpsData %>%
          mutate(trtlabel = ifelse(treatgroup == 1, 'ICS', 'LABA/LAMA')) %>%
          ggplot(aes(x = iptw_weight, color = trtlabel)) +  # Change color mapping to trtlabel
          geom_density(alpha = 0.5, linewidth = 1) +
          scale_x_continuous(breaks = seq(0, 10, 1)) +
          xlim(0, 2) +
          labs(x = 'Inverse probability of treatment weight', y = 'Density', title = 'Weight distribution by treatment group') +
          scale_color_manual(values = c("ICS" = palette[9], "LABA/LAMA" = palette[4])) +  # Correctly referencing palaette colors
          theme_minimal() +
          theme(legend.title = element_blank()) +
          guides(color = guide_legend(override.aes = list(fill = c(palette[9], palette[4]))))
        
        file_path <- file.path(HDPS_folder,
                               "outputs",
                               paste0("HDPS_weight_dist_", k, outcome, cohort_ext, "_trimmed.png"))
        
        ggsave(file_path, plot, width = 8, height = 4)
        
        # Create a survey design object
        svy_design <- svydesign(ids = ~ 1,
                                weights = ~ iptw_weight,
                                data = hdpsData)
        
        # Create the weighted table
        tabWeighted <- svyCreateTableOne(
          vars = vars,
          strata = "treatgroup",
          data = svy_design,
          test = FALSE,
          factorVars = catVars)
        
        # Export TableOne to CSV
        write.csv(
          print(tabWeighted,
                smd = TRUE,
                quote = FALSE,
                noSpaces = TRUE,
                printToggle = FALSE),
          file = paste0(HDPS_folder,
                        "/outputs/covariates_HDPS_",
                        k,
                        "_table1Weighted_",
                        outcome,
                        cohort_ext,
                        "_trimmed.csv"))
        
        # Plot Overlap with IPTW
        iptw_plot <- hdpsData %>%
          mutate(trtlabel = ifelse(treatgroup == 1, 'ICS', 'LABA/LAMA')) %>%
          ggplot(aes(x = pscore, color = trtlabel)) + # Switched to color aesthetic
          geom_density(alpha = 0.5, linewidth = 1) +
          labs(x = 'Probability of receiving ICS', # Updated labels
               y = 'Density', title = 'Overlap plot with IPTW') +
          scale_color_manual(values = c("ICS" = palette[9], "LABA/LAMA" = palette[4])) + # Custom color scheme
          scale_x_continuous(breaks = seq(0, 1, 0.1)) +
          theme_minimal() +
          theme(legend.title = element_blank()) + # Remove legend title
          guides(color = guide_legend(override.aes = list(fill = c(
            palette[9], palette[4]
          ))))
        
        ggsave(paste0(HDPS_folder,
                      "/outputs/covariates_HDPS_",
                      k,
                      "_overlapWeighted_",
                      outcome,
                      cohort_ext,
                      "_trimmed.png"), height = 4, width = 6)
        
        
        #merge in iptw_weight_predefined from hdpsPredefinedVars
        hdpsData <- hdpsData %>%
          left_join(hdpsPredefinedVars %>%
                      dplyr::select(patid, iptw_weight_predefined),
                    by = "patid")
        
        #save patid and iptw_weight to a new object
        hdpsData_iptw_hdps <- hdpsData %>%
          dplyr::select(patid, iptw_weight_predefined)
        #remove iptw_weight_predefined from hdpsData
        hdpsData <- hdpsData %>%
          dplyr::select(-iptw_weight_predefined)
        
        # Merge the iptw_weight_predefined into hdpsData
        hdpsData <- hdpsData %>%
          left_join(hdpsPredefinedVars %>%
                      select(patid, iptw_weight_predefined),
                    by = "patid")
        
        #rmeove people with missing weighs
        hdpsData <- hdpsData %>%
          filter(!is.na(iptw_weight_predefined))
        # Create a survey design object using predefined weights but including all variables
        svy_design_predefined <- svydesign(ids = ~ 1,
                                           weights = ~ iptw_weight_predefined,
                                           data = hdpsData)
        
        
        #count empty weights
        empty_weights <- sum(is.na(svy_design_predefined$weights))
        
        # Create the weighted table using all variables
        tabWeighted_predefined <- svyCreateTableOne(
          vars = vars,  # This now includes all variables, both predefined and HDPS
          strata = "treatgroup",
          data = svy_design_predefined,
          test = FALSE,
          factorVars = catVars)
        
        # Export TableOne to CSV
        write.csv(
          print(tabWeighted_predefined,
                smd = TRUE,
                quote = FALSE,
                noSpaces = TRUE,
                printToggle = FALSE),
          file = paste0(HDPS_folder,
                        "/outputs/covariates_HDPS_",
                        k,
                        "_table1Weighted_predefined_",
                        outcome,
                        cohort_ext,
                        "_trimmed.csv"))
        
        extract_smd <- function(tableone_obj, label) {
          smd_data <- as.data.frame(print(tableone_obj, smd = TRUE, printToggle = FALSE))
          smd_data$Variable <- rownames(smd_data)
          smd_data$Iteration <- label
          smd_data <- smd_data %>%
            dplyr::select(Iteration, Variable, SMD)
          return(smd_data)
        }
        
        smd_combined <- extract_smd(tabUnweighted, paste0("Unweighted_", k, "_", outcome)) %>%
          rename(smd_unweighted = SMD) %>%
          left_join(
            extract_smd(tabWeighted, paste0("Weighted_", k, "_", outcome)) %>%
              rename(smd_weighted = SMD),
            by = c("Variable")
          ) %>%
          left_join(
            extract_smd(
              tabWeighted_predefined,
              paste0("Predefined_Weighted_", k, "_", outcome)
            ) %>%
              rename(smd_weighted_predefined = SMD),
            by = c("Variable")
          ) %>%
          mutate_at(vars(contains("smd")), ~ ifelse(. == "<0.001", 0.0001,.)) %>%
          mutate(
            analysis = paste0(k, "_", outcome),
            smd_unweighted = as.numeric(smd_unweighted),
            smd_weighted = as.numeric(smd_weighted),
            smd_weighted_predefined = as.numeric(smd_weighted_predefined)
          ) %>%
          select(analysis,
                 Variable,
                 smd_unweighted,
                 smd_weighted,
                 smd_weighted_predefined)
        
        # Save SMDs to a CSV file
        write_csv(smd_combined,
                  paste0(HDPS_folder, "/outputs/HDPS_", k, "_SMD_", outcome, cohort_ext, "_trimmed.csv"))
        
        
        #save the HDPS data to a parquet file
        write_parquet(hdpsData, paste0(HDPS_folder, "/outputs/HDPS_", k, "_", outcome, cohort_ext, "_trimmed.parquet"))
        
      }
      
    }
  }
  
}

gc()
# Set up parallel processing
numCores <- availableCores()
plan(multisession, workers = (numCores-1))

# Create a tibble with all combinations of parameters
params <- expand_grid(
  cohort_ext = cohort_exts,
  outcome = outcomes,
  topVar = topVars
)

# Use future_pmap_dfr to parallelize the function
results <- future_pmap_dfr(
  params,
  ~cov_weighting(..1, ..2, ..3)
)
