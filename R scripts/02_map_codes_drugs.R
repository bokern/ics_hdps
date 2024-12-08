#############################
# Project: HDPS
# Programmer name: Marleen Bokern
# Date Started:  06/24
#############################
# Map prodcodes to BNF paragraphs
#############################
packages <- c(
  "tidyverse",
  "MetBrewer",
  "arrow",
  "readstata13",
  "data.table",
  "autoCovariateSelection"
)
installed_packages <- packages %in% rownames(installed.packages())
if (any(installed_packages == FALSE)) {
  install.packages(packages[!installed_packages])
}

invisible(lapply(packages, function(pkg) {
  suppressMessages(library(pkg, character.only = TRUE, verbose = FALSE))
}))

palette <- met.brewer("Cassatt2")

#unzip files
setwd(Copd_aurum_extract)

#code_browser as string variables
file_path <- file.path(Denominator,
                       "CPRD Aurum",
                       "Code browsers",
                       "2022_05" ,
                       "CPRDAurumProduct.txt")
code_browser <- read.delim(file_path,
                           sep = "\t",
                           header = TRUE,
                           colClasses = "character") %>% dplyr::select(ProdCodeId, dmdid, BNFChapter, Term.from.EMIS, DrugIssues)
# Convert DrugIssues to numeric and format without scientific notation
code_browser$DrugIssues <- format(as.numeric(code_browser$DrugIssues), scientific = FALSE)

# Read in BNF SNOMED mapping file (only needs to be done once, so moved outside the loop)
file_path <- file.path(HDPS_folder, "bnf_to_dmd.csv")
bnf_snomed <- read.csv(file_path, header = TRUE, colClasses = "character") %>% dplyr::select(-X, -vmp_previous, -vmp_previous_date, -type, -vtm, -nm, -vtm_nm)

#left join code_browser with bnf_snomed
code_map <- code_browser %>%
  left_join(bnf_snomed, by = c("dmdid" = "id"))

#create bnf_paragraph
code_map$bnf_paragraph <- case_when(
  !is.na(code_map$bnf_code) &
    nchar(code_map$bnf_code) >= 6 ~ substr(code_map$bnf_code, 1, 6),!is.na(code_map$BNFChapter) &
    nchar(code_map$BNFChapter) == 7 ~ str_pad(
      substr(code_map$BNFChapter, 1, 5),
      6,
      side = "left",
      pad = "0"
    ),!is.na(code_map$BNFChapter) &
    nchar(code_map$BNFChapter) == 8 ~ substr(code_map$BNFChapter, 1, 6),
  TRUE ~ NA_character_
)

#pick out codes with bnf_paragraph = NA and sort on DrugIssues
missing <- code_map %>%
  filter(is.na(bnf_paragraph)) %>%
  arrange(desc(DrugIssues)) %>%
  dplyr::select(ProdCodeId, dmdid, BNFChapter, Term.from.EMIS, DrugIssues)
#make drugissues numeric
code_map$DrugIssues <- as.numeric(code_map$DrugIssues)
#sort on bnf_paragraph and DrugIssues
code_map <- code_map %>%
  arrange(bnf_paragraph, desc(DrugIssues))
#count if bnf_paragraph = NA
nrow(code_map[is.na(code_map$bnf_paragraph), ])
nrow(code_map[is.na(code_map$bnf_paragraph) & (code_map$DrugIssues >= 500000), ])

#subset of code_map with bnf_paragraph = NA
na_code_map <- code_map[is.na(code_map$bnf_paragraph), ]

ggplot(na_code_map, aes(x = DrugIssues)) +
  geom_histogram(bins = 150, fill = palette[1], color = "black") +
  labs(
    title = "Histogram of DrugIssues for unmapped codes",
    x = "DrugIssues",
    y = "Count"
  ) +
  theme_minimal() +
  scale_y_log10() +
  #make x axis labels readable
  scale_x_continuous(labels = scales::comma)

# Map specific terms to BNF paragraphs when bnf_paragraph is empty
code_map$bnf_paragraph <- case_when(
  is.na(code_map$bnf_paragraph) &
    grepl("paracetamol.*mg.*", code_map$Term.from.EMIS, ignore.case = TRUE) ~ "040701",
  is.na(code_map$bnf_paragraph) &
    grepl(
      "mirtazapine.*tablets",
      code_map$Term.from.EMIS,
      ignore.case = TRUE
    ) ~ "040304",
  is.na(code_map$bnf_paragraph) &
    grepl("apixaban.*tablets.*", code_map$Term.from.EMIS, ignore.case = TRUE) ~ "020802",
  is.na(code_map$bnf_paragraph) &
    grepl("aspirin.*mg", code_map$Term.from.EMIS, ignore.case = TRUE) ~ "020900",
  is.na(code_map$bnf_paragraph) &
    grepl(
      "rosuvastatin.*tablets",
      code_map$Term.from.EMIS,
      ignore.case = TRUE
    ) ~ "021200",
  is.na(code_map$bnf_paragraph) &
    grepl(
      "candesartan.*tablets",
      code_map$Term.from.EMIS,
      ignore.case = TRUE
    ) ~ "020505",
  is.na(code_map$bnf_paragraph) &
    grepl("quinine bisulfate", code_map$Term.from.EMIS, ignore.case = TRUE) ~ "050401",
  is.na(code_map$bnf_paragraph) &
    grepl("epimax.*cream", code_map$Term.from.EMIS, ignore.case = TRUE) ~ "212200",
  is.na(code_map$bnf_paragraph) &
    grepl(
      "sitagliptine.*tablets",
      code_map$Term.from.EMIS,
      ignore.case = TRUE
    ) ~ "060102",
  is.na(code_map$bnf_paragraph) &
    grepl("ventolin", code_map$Term.from.EMIS, ignore.case = TRUE) ~ "030101",
  is.na(code_map$bnf_paragraph) &
    grepl(
      "amlodipin.*tablets.*mg",
      code_map$Term.from.EMIS,
      ignore.case = TRUE
    ) ~ "020602",
  is.na(code_map$bnf_paragraph) &
    grepl("gaviscon", code_map$Term.from.EMIS, ignore.case = TRUE) ~ "010102",
  is.na(code_map$bnf_paragraph) &
    grepl("erythromycin", code_map$Term.from.EMIS, ignore.case = TRUE) ~ "050105",
  is.na(code_map$bnf_paragraph) &
    grepl(
      "fluticasone.*nasal spray",
      code_map$Term.from.EMIS,
      ignore.case = TRUE
    ) ~ "120201",
  is.na(code_map$bnf_paragraph) &
    grepl("salbutamol aerosol", code_map$Term.from.EMIS, ignore.case = TRUE) ~ "030101",
  is.na(code_map$bnf_paragraph) &
    grepl("co-codamol", code_map$Term.from.EMIS, ignore.case = TRUE) ~ "040701",
  TRUE ~ code_map$bnf_paragraph
)

#ad hoc changes after examining most important covariates in hdps: correct mapping for prednisolone
code_map$bnf_paragraph <- case_when(
  grepl("1114241000033116", code_map$ProdCodeId, ignore.case = TRUE) ~ "060302",
  grepl("1131941000033119", code_map$ProdCodeId, ignore.case = TRUE) ~ "060302",
  TRUE ~ code_map$bnf_paragraph
)

#save codes that remain unmatched to drugs_counts_unmatched.csv
unmatched <- code_map %>%
  filter(is.na(bnf_paragraph)) %>%
  dplyr::select(ProdCodeId, dmdid, BNFChapter, Term.from.EMIS, DrugIssues) %>%
  #select the top 100 unmatched codes by DrugIssues
  arrange(desc(DrugIssues)) %>%
  head(100)

#save to file
write.csv(
  unmatched,
  file = file.path(HDPS_folder, "drugs_counts_unmatched.csv"),
  row.names = FALSE
)

code_map_exp <- code_map %>%  filter(!is.na(bnf_paragraph)) %>%
  dplyr::select(ProdCodeId, bnf_paragraph) %>%
  distinct()

# Ensure ProdCodeId is character and format without scientific notation
code_map_exp$ProdCodeId <- format(as.numeric(code_map_exp$ProdCodeId), scientific = FALSE)

# Convert bnf_paragraph to character explicitly
code_map_exp$bnf_paragraph <- as.character(code_map_exp$bnf_paragraph)

# Write to file, saving all variables as text
file_path <- file.path(HDPS_folder, "map_prodcode_bnf.txt")

#remove spaces
code_map_exp$ProdCodeId <- gsub(" ", "", code_map_exp$ProdCodeId)

# Write the data frame to txt, ensuring no scientific notation and all columns as text
write.table(
  code_map_exp,
  file = file_path,
  row.names = FALSE,
  quote = TRUE,
  fileEncoding = "UTF-8"
)

#remove ICS, LABA and LAMA from drugs_all
# Function to read the first part of each line from a text file
read_first_part_of_lines <- function(filepath) {
  lines <- readLines(filepath)
  parts <- strsplit(lines, "\\s+") # Split each line at spaces
  first_parts <- sapply(parts, `[`, 1)    # Take the first part of each split line
  # Remove any non-numeric characters and convert to numeric
  numeric_parts <- as.character(gsub("[^0-9]", "", first_parts))
  return(numeric_parts)
}

# Apply the function to each file
ics_single_codes <- read_first_part_of_lines(file.path(Github_codelists, "cl_ics_single.txt"))
laba_single_codes <- read_first_part_of_lines(file.path(Github_codelists, "cl_laba_single.txt"))
lama_single_codes <- read_first_part_of_lines(file.path(Github_codelists, "cl_lama_single.txt"))
laba_lama_codes <- read_first_part_of_lines(file.path(Github_codelists, "cl_laba_lama.txt"))
ics_laba_codes <- read_first_part_of_lines(file.path(Github_codelists, "cl_ics_laba.txt"))
triple_therapy_codes <- read_first_part_of_lines(file.path(Github_codelists, "cl_triple_therapy.txt"))

# Combine all codes into one vector
all_codes <- c(
  ics_single_codes,
  laba_single_codes,
  lama_single_codes,
  laba_lama_codes,
  ics_laba_codes,
  triple_therapy_codes
)

setwd(Copd_aurum_extract)
# List of drug files to process
drug_files <- list.files(Copd_aurum_extract,
                         pattern = "Drugs_for_hdps_.*\\.parquet",
                         full.names = TRUE)

# Initialize an empty dataframe to store all drugs_empty_bnf
all_drugs_empty_bnf <- data.frame()
drugs_all <- data.frame()

for (file in drug_files) {
  # Extract file number for naming
  file_num <- str_extract(basename(file), "\\d+")
  
  # Read and process each drug file
  drugs <- arrow::read_parquet(file) %>% dplyr::select(-pracid, -enterdate)
  drugs <- setDT(drugs)
  #remove ICS, LABA and LAMA from drugs_all
  drugs <- drugs[!(prodcodeid %in% all_codes)]
  # Merge with code_browser
  drugs <- drugs %>%
    left_join(code_map, by = c("prodcodeid" = "ProdCodeId")) %>%
    dplyr::select(patid, bnf_paragraph, issuedate) %>%
    #rename BNFchapter to code
    rename(code = bnf_paragraph)
  
  # Save cut of data with empty BNF Chapter and BNF Code
  drugs_empty_bnf <- drugs %>% filter(is.na(drugs$code))
  
  # Append drugs_empty_bnf to all_drugs_empty_bnf
  all_drugs_empty_bnf <- rbind(all_drugs_empty_bnf, drugs_empty_bnf)
  
  #remove rows with no bnf_paragraph from drugs_all
  drugs <- drugs %>% filter(!is.na(code))
  drugs_all <- rbind(drugs_all, drugs)
}

# After the loop, save the combined all_drugs_empty_bnf
write_parquet(all_drugs_empty_bnf, "all_drugs_empty_bnf.parquet")

#drop bnf_code and BNFChapter from drugs_all
drugs_all <- drugs_all %>% dplyr::select(c(patid, code, issuedate))

setwd(Datadir_copd)

write_parquet(drugs_all, "drugs_for_hdps_mapped.parquet")
