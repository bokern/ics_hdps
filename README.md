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
