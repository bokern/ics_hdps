# -----------------------------------------------------------------------------
# PROGRAM NAME:  table_SMDs_predefined_covariates.R
# AUTHOR:        Marleen Bokern
# DATE CREATED:   Oct 2024
# NOTES:       Generates tables of Standardized Mean Differences (SMDs) for predefined covariates. It compares unweighted and weighted (using IPTW) SMDs for different outcomes and cohorts. The SMDs are calculated for the predefined covariates before weighting, and after weighitng using the predefined covariates, and different iterations of the HDPS
#                
# -----------------------------------------------------------------------------


packages <- c(
  "tidyverse",
  "arrow",
  "dbplyr",
  "lubridate",
  "scales",
  "tableone",
  "Matching",
  "survival",
  "MetBrewer",
  "gridExtra",
  "survminer",
  "survey",
  "furrr",
  "future",
  "parallel",
  "WeightIt"
)

installed_packages <- packages %in% rownames(installed.packages())
if (any(installed_packages == FALSE)) {
  install.packages(packages[!installed_packages])
}

invisible(lapply(packages, function(pkg) {
  suppressMessages(library(pkg, character.only = TRUE, verbose = FALSE))
}))

outcomes <- c("covid_hes_present", "covid_death_present")
cohort_exts <- c("", "_no_triple")
topVars <- c("100", "250", "500", "750", "1000")

for (outcome in outcomes) {
  for (cohort_ext in cohort_exts) {
    timeinstudy <- switch(outcome,
                          "covid_hes_present" = "timeinstudy2",
                          "covid_death_present" = "timeinstudy3")
    
    #read in this file   file_path_predefined <- file.path(
    data_path <- file.path(
      HDPS_folder,
      "outputs",
      paste0(
        "HDPS_predefined_",
        outcome,
        cohort_ext,
        "_trimmed.parquet"
      )
    )
    
    
    # Load the data
    hdpsData <- arrow::read_parquet(data_path)
    
    
    # Vector of categorical variables that need transformation
    catVars <- setdiff(
      names(hdpsData),
      c(
        "patid",
        "age_index",
        "treatgroup",
        "iptw_weight_predefined",
        "pscore_predefined",
        timeinstudy
      )
    )
    
    
    #data for unweighted table
    # Combine predefined and HDPS covariates if hdpsData exists
    all_vars_table <- setdiff(
      names(hdpsData),
      c(
        "patid",
        "treatgroup",
        "iptw_weight_predefined",
        "pscore_predefined",
        timeinstudy
      )
    )
    
    tryCatch({
      tabUnweighted <- CreateTableOne(
        vars = all_vars_table,
        strata = "treatgroup",
        data = hdpsData,
        test = FALSE,
        smd = TRUE,
        factorVars = catVars
      )
    }, error = function(e) {
      print(paste("Error occurred:", e$message))
    })
    
    #create df from tableone object
    smd_data <- as.data.frame(print(tabUnweighted, smd = TRUE, printToggle = FALSE))
    smd_data <- smd_data %>%
      tibble::rownames_to_column(var = "covariate")
  
    
    svy_design_predefined <- svydesign(
      ids = ~ 1,
      weights = ~ iptw_weight_predefined,
      data = hdpsData
    )
    
    #generate smd table weighted using predefined covariates
    tabWeighted_predefined <- svyCreateTableOne(
      vars = all_vars_table,
      strata = "treatgroup",
      data = svy_design_predefined,
      test = FALSE,
      smd = TRUE,
      factorVars = catVars
    )
    
    smd_data_predefined <- as.data.frame(print(
      tabWeighted_predefined,
      smd = TRUE,
      printToggle = FALSE
    ))
    smd_data_predefined <- smd_data_predefined %>%
      tibble::rownames_to_column(var = "covariate") %>%
      rename(smd_predefined = SMD)
    
    
    # take the smd column from tabWeighted_predefined, rename it to smd_predefined, and add it to smd_data
    smd_data <- smd_data %>%
      left_join(smd_data_predefined %>%
                  dplyr::select(covariate, smd_predefined),
                by = "covariate")
    
    
    
    for (k in topVars) {
      #read in hdpsData with weights
      data_path <- file.path(
        HDPS_folder,
        "outputs",
        paste0("HDPS_", k, "_", outcome, cohort_ext, "_trimmed.parquet")
      )
      hdpsData <- arrow::read_parquet(data_path)
      
      #remove all vriables that start with d1_, d2_ or d3_
      hdpsData <- hdpsData %>%
        dplyr::select(-starts_with("d1_"),
                      -starts_with("d2_"),
                      -starts_with("d3_")) %>%
        dplyr::select(-iptw_weight_predefined, -pscore)
      
      # Vector of categorical variables that need transformation
      catVars <- setdiff(
        names(hdpsData),
        c(
          "patid",
          "age_index",
          "treatgroup",
          "iptw_weight",
          "pscore",
          "baseline_triple",
          outcome
        )
      )
      
      
      #data for unweighted table
      # Combine predefined and HDPS covariates if hdpsData exists
      all_vars_table <- setdiff(
        names(hdpsData),
        c(
          "patid",
          "treatgroup",
          "iptw_weight",
          "baseline_triple",
          "pscore",
          timeinstudy,
          outcome
        )
      )
      
      svy_design_hdps <- svydesign(ids = ~ 1,
                                   weights = ~ iptw_weight,
                                   data = hdpsData)
      
      #generate smd table weighted using predefined covariates
      tabWeighted_hdps <- svyCreateTableOne(
        vars = all_vars_table,
        strata = "treatgroup",
        data = svy_design_hdps,
        test = FALSE,
        smd = TRUE,
        factorVars = catVars
      )
      
      smd_var_name <- paste0("smd_hdps_", k)
      
      smd_data_hdps <- as.data.frame(print(
        tabWeighted_hdps,
        smd = TRUE,
        printToggle = FALSE
      ))
      smd_data_hdps <- smd_data_hdps %>%
        tibble::rownames_to_column(var = "covariate") %>%
        rename(!!smd_var_name := SMD)
      
      # take the smd column from tabWeighted_predefined, rename it to smd_predefined, and add it to smd_data
      smd_data <- smd_data %>%
        left_join(smd_data_hdps %>%
                    dplyr::select(covariate, !!smd_var_name),
                  by = "covariate")
      
    }
    
    # Save overall SMDs to a CSV file
    write_csv(
      smd_data,
      paste0(
        HDPS_folder,
        "/outputs/SMDs_predefined_covariates",
        outcome, 
        "_",
        cohort_ext,
        ".csv"
      )
    )
    
  }
}