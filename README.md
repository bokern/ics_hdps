# ics_hdps
 Code for HDPS analysis of ICS and COVID outcomes. 
 Data management code is available under https://github.com/bokern/ics_covid/tree/main/do-files. Codelists are available under https://github.com/bokern/ics_covid/tree/main/codelists

The analysis pipeline consists of the main scripts:

1. Data Extraction (01_extract_codes.R)
- Extracts and processes relevant patient data during a 12-month lookback period
- Reads observation and prescription files for patients in the COPD cohort
- Filters data by date ranges (lookback period before 12 months prior to March 2020)
- Processes HES inpatient diagnosis data

2. Code Mapping (02_map_codes.R)
- Processes primary care observation data
- Maps SNOMED CT codes from the CPRD observation files to ICD-10 codes to avoid sparsity
- Merges CPRD browser data with SNOMED-to-ICD10 mappings
- Filters observations to the relevant study period
- Creates an "ever-mapped" dataset capturing unique ICD-10 codes per patient
- Exports both mapped and unmatched observations

  02_map_codes_drugs.R
- Processes primary care prescription data
- Maps CPRD prodcodes to BNF chapters
- Filters drug issues to the covariate assessment window.
- Manually maps the most frequent unmapped product codes
- Exports both mapped and unmatched observations

  02_map_codes_hes.R
- Truncates ICD-10 codes to 3 characters

3. Recurrence Assessment and HDPS Implementation (03_assess_recurrence_multioutcome.R)
- Calculates code prevalence and recurrence within each healthcare dimension
- Identifies eligible covariates by excluding predefined variables
- Ranks codes by their potential for confounding bias
- Applies a bias formula to prioritize codes that may confound treatment-outcome relationships
- Creates datasets with different numbers of top-ranked covariates (100, 250, 500, 750, 1000)
- Implements instrumental variable checks to exclude potential instruments

4.parallelized propensity score weighting (04_cov_weighting_parallel.R)

- Calculates and applies inverse probability of treatment weights (IPTW)
- Analyzes both predefined and high-dimensional covariates
- Utilizes parallel processing to improve performance
- Generates diagnostic plots and tables for assessing balance

5. Treatment effect estimation (05_effect_estimation.R)

- estimates treatment effects for COPD patients using different weighting methods: 
  - Unweighted analysis
  - Analysis using prespecified covariates
  - Analysis using high-dimensional propensity score (HDPS) with various covariate counts (100, 250, 500, 750, 1000)
- Input Data: cohort data from parquet files for COPD patients in Wave 1 (60 days)
- Analysis Methods: Cox Proportional Hazards Models for time-to-event analysis, Logistic Regression Models for binary outcome analysis, Kaplan-Meier Curves with risk tables and cumulative events
- Diagnostics: Schoenfeld residuals for Cox models, Residual plots for logistic regression models, Propensity score balance assessment
- Outputs: Detailed results tables for Cox and logistic regression models, Kaplan-Meier curve visualizations with risk tables, Diagnostic plots for model assessment, Summarized results in CSV format
Usage:
  - Set the working directory and debugging flag as needed. The script will loop through:
  - Different outcomes (COVID-19 hospitalization, COVID-19 death)
  - Different cohorts (all patients, non-triple therapy patients)
  - Different HDPS variable selection thresholds (100-1000 variables)



HDPS Diagnostic Visualizations
This folder contains R scripts for generating diagnostic visualizations for high-dimensional propensity score (HDPS) analyses. These scripts help assess the quality and performance of HDPS models by creating informative plots and tables.
Overview
The code in this repository generates various diagnostic plots for evaluating propensity score models, with a focus on comparing conventional predefined covariates with high-dimensional automated approaches. The visualizations help assess balance, prevalence, and treatment effects across different covariate sets.
Key Visualizations
1. Covariate Concepts Plot (chapters_prevalence.R)
Visualizes the prevalence of high-level clinical concepts within the top N ranked HDPS covariates:
- Creates circular bar plots showing the distribution of clinical codes across different data dimensions (Clinical, Prescriptions, Hospital)
- Groups codes by clinical chapters (e.g., respiratory system, cardiovascular system)
- Provides a visual understanding of which clinical domains are being captured by the HDPS algorithm

2. Standardized Mean Differences Plot (stdDiffsPredefinedPlusHDPS.R)
- Generates plots comparing standardized mean differences (SMDs) before and after weighting:
- Shows covariate balance between treatment groups
- Compares unweighted vs. weighted balance using different approaches
- Highlights the relative impact of predefined covariates vs. HDPS covariates
- Demonstrates how adding HDPS covariates improves balance across groups

3. Combined PS Overlap Plot (combined_plot_for_paper.R)
- Creates publication-ready plots that combine:
  - Propensity score distributions for treatment groups
  - Compares overlap between unweighted and weighted scenarios
  - Arranges multiple plots in a grid for easy comparison
  - Customizes visualization settings for publication quality

4. Treatment Effect by Covariate Count (plot_or_covariates_one_by_one.R)
  - Generates plots showing how treatment effect estimates change with increasing numbers of covariates (1-250 HDPS-covariates):
  - Plots hazard/odds ratios against number of covariates included (both hazard and odds ratios)
  - Shows confidence intervals as lines and ribbons
  - Helps identify when effect estimates stabilize with additional covariates

5. SMD Tables for Predefined Covariates (table_SMDs_predefined_predefined_covariates.R)
  - Creates detailed tables showing standardized mean differences:
  - Compares balance across predefined covariates only
  - Shows unweighted and weighted balance metrics
  - Allows comparison of balance achieved with different HDPS covariate sets
  - Outputs results to CSV files for further analysis

Usage
The scripts are designed to work with propensity score analysis outputs from HDPS models. Common configuration parameters include:

k: Number of covariates (e.g., 100, 250, 500, 750, 1000)
outcome: Outcome variable (e.g., "covid_hes_present", "covid_death_present")
cohort_ext: Cohort extension (e.g., "", "_no_triple")

Output
The visualizations are saved as:
- PNG files for direct viewing
- RDS files for further modification and analysis
- CSV files for tabular data

- All outputs are stored in an "outputs/diagnostics" directory structure.

