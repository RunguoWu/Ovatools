#################################################################
##      generate study population and index date               ##
#################################################################
# this file generates the list of study population
# each line represents a unique ID with an index date
# it is based on two files shared from Kirsten, who implemented her algorithm
# and code list to generate all TV/pelvic US records and all CA125 records
# I include the earliest diagnostic US records and CA125 records between 2013-04-01 and 2017-12-31
# the index date is whichever the earliest for each
# I exclude those included with a procedural US records within 1 year before diagnostic US (incl. the same day)
# and exclude those included with any TV/pelvic US records and CA125 1 year before indedx date

rm(list = ls())

source("//Qmcr/dfs/WIPH/Primary Care/OVATOOLS/r_project/scripts/file_path.R")
source(file.path(r_wd, "scripts/fn_data_clean.R"))

# Read in stata data prepared by Kirsten
# us_all <- read_dta(file.path(proc_data, "US_alldates.dta"))
# ca125_all <- read_dta(file.path(proc_data, "CA125_alldates.dta"))
# Save a RDS version
# saveRDS(us_all, file.path(work_data, "US_alldates.RDS"))
# saveRDS(ca125_all, file.path(work_data, "CA125_alldates.RDS"))

# Step 1:  ----------------------------------------------------------------

## implement inclusion and exclusion criteria ----
# read the US and CA125 orginal data
us_all <- readRDS(file.path(work_data, "US_alldates.RDS"))
ca125_all <- readRDS(file.path(work_data, "CA125_alldates.RDS"))

he_ova <- inc_excl(us_all = us_all, 
                   ca125_all = ca125_all, 
                   start_date="2013-04-01", 
                   end_date="2017-12-31", 
                   symptom=T)

## Add demographics ----
# DOB data from CPRD
# demo_cprd <- read_dta(file.path(proc_data, "demographics_agesex_all.dta"))
# saveRDS(demo_cprd, file.path(work_data, "demographics.RDS"))
demo_cprd <- readRDS(file.path(work_data, "demographics.RDS"))

# Townsend 2011, 1 = least deprived
townsend <- read_dta(file.path(raw_link_data, "e_aurum_townsend2011_10.dta"))
# HES patient characteristics
hes_pat <- read_dta(file.path(raw_link_data, "e_aurum_hes_patient_21_001655_dm.dta"))
# CPRD ethnicity
# ethn_cprd <- read_dta(file.path(proc_data, "ethnicity_cprd_alldates.dta"))
# saveRDS(ethn_cprd, file.path(work_data, "ethnicity_cprd_alldates.RDS"))
ethn_cprd <- readRDS(file.path(work_data, "ethnicity_cprd_alldates.RDS"))


he_ova <- add_demo(he_ova = he_ova, 
                   demo_cprd = demo_cprd, 
                   townsend = townsend, 
                   ethn_cprd = ethn_cprd, 
                   hes_pat = hes_pat)

## Add death records ----

# linked death records
death <- read_dta(file.path(raw_link_data, "e_aurum_death_patient_21_001655_dm.dta"))
he_ova <- add_death(he_ova = he_ova, death = death)

## save ----
# saveRDS(he_ova, file.path(work_data, "he_ova_step1_withSymptom.rds"))
saveRDS(he_ova, file.path(work_data, "he_ova_step1_withSymptom_prescr.rds"))
rm(demo_cprd, townsend, hes_pat, ethn_cprd, us_all, ca125_all)

# Step 2:  ----------------------------------------------------------------
## prepare HES cancer/gynae diagnosis records by admission----
# linked hospital admission records
hes_hosp <- read_dta(file.path(raw_link_data, "e_aurum_hes_hospital_21_001655_dm.dta"))
# linked diagnoses in hospital
hes_hosp_diag <- read_dta(file.path(raw_link_data, "e_aurum_hes_diagnosis_hosp_21_001655_dm.dta"))
# hes_hosp has 79 admission dates are missing
# use episode starte date to fill
hes_epis <- read_dta(file.path(raw_link_data, "e_aurum_hes_episodes_21_001655_dm.dta"))

# ICD codes
# benign gynae
benign_gynae_code <- read_dta(file.path(gen_wd, "Codelists/CPRD Aurum - benign gynae disease.dta"))
benign_gynae_code <- paste(unique(benign_gynae_code$icd10), collapse = "|^")
benign_gynae_code <- paste0("^", benign_gynae_code)
# ICD 10 code refer to Garth et al. 2020, and meeting discussion
icd_codes <- c(
  ova = "^C56|^C57.0|^C48.1|^C48.2",   
  panc = "^C25",
  lung = "^C34",
  uter = "^C54|^C55",
  loGI = "^C18|^C19|^C20|^C21",
  benign = benign_gynae_code, # include D39.1
  borderline = "^D39.1"
)

hes_admin_diag <- gen_admin_diag(hes_hosp = hes_hosp, 
                                 hes_hosp_diag = hes_hosp_diag, 
                                 hes_epis = hes_epis, 
                                 icd_codes = icd_codes)

saveRDS(hes_admin_diag, file.path(work_data, "hes_admin_diag.rds"))
rm(hes_hosp, hes_hosp_diag, hes_epis)

# Step 3:  ----------------------------------------------------------------
# read in the study population dataset
he_ova <- readRDS(file.path(work_data, "he_ova_step1_withSymptom_prescr.rds"))

# only keep the ids in the study population
hes_admin_diag <- readRDS(file.path(r_wd, "working_data", "hes_admin_diag.rds"))

# CPRD cancer, from Kirsten
# allcancers_cprd <- read_dta(file.path(proc_data, "cprd_allcancers_clean.dta"))
# saveRDS(allcancers_cprd, file.path(work_data, "cprd_allcancers_clean.dta"))
allcancers_cprd <- readRDS(file.path(work_data, "cprd_allcancers_clean.dta"))

he_ova <- exclude_early_oc(he_ova=he_ova, 
                           hes_admin_diag=hes_admin_diag, 
                           allcancers_cprd=allcancers_cprd)

he_ova <- pat_diag(he_ova=he_ova, 
                   hes_admin_diag=hes_admin_diag)

he_ova <- add_history(he_ova=he_ova, 
                      hes_admin_diag=hes_admin_diag, 
                      allcancers_cprd=allcancers_cprd)

he_ova <- add_incident(he_ova)

# since apply US+symptoms 2023-10
# saveRDS(he_ova, file.path(r_wd, "working_data", "he_ova_step3_withSymptom.rds"))
# write_dta(he_ova, file.path(r_wd, "working_data", "he_ova_step3_withSymptom.dta"))

saveRDS(he_ova, file.path(r_wd, "working_data", "he_ova_step3_withSymptom_prescr.rds"))
write_dta(he_ova, file.path(r_wd, "working_data", "he_ova_step3_withSymptom_prescr.dta"))

rm(allcancers_cprd)


