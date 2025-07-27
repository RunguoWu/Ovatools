#################################################################
##             prepare cost data for analysis                  ##
#################################################################
rm(list = ls())

source("//Qmcr/dfs/WIPH/Primary Care/OVATOOLS/r_project/scripts/file_path.R")
source(file.path(r_wd, "scripts/fn_data_clean.R"))
source(file.path(r_wd, "scripts/cost/fn_hrg_prep.R"))

# HRG prep ----------------------------------------------------------------

hes_hrg <- read_dta(file.path(raw_link_data, "e_aurum_hes_hrg_21_001655_dm.dta"))
# clean hrg raw data:
# keep year>=2012 
# use at least 1 year before diagnosis for control
# use 2011 to avoid some admission in 2012 use 2011 HRG
# at admission/spell level
# remove duplicates
# impute missing with available episode hrg
hes_hrg_clean <- hrg_clean(hes_hrg, 2011)
saveRDS(hes_hrg_clean, file.path(work_data, "hes_hrg_clean.rds"))

# associate with costs
# all ids since 2013
# do not limited to study population
hes_hrg_cost <- hrg_cost(hes_hrg_clean)

# inflated to 21/22
inflator <- read_excel(file.path(work_data, "NHS inflation rates.xlsx"), sheet = 1)

hes_hrg_cost <- hes_hrg_cost %>%
  left_join(inflator) %>%
  mutate(Unit_cost_inf= Unit_cost * inflator)

saveRDS(hes_hrg_cost, file.path(work_data, "hes_hrg_cost.rds"))


# create study population for cost ----------------------------------------
# same method as create the study population for the HE study
# but to estimate cost, we do not need to stick to the clinical pathway
# instead, we should include cancer incidents as many as possible
# so do not limited to people with symptoms before US.

## Step 1:  ----

# implement inclusion and exclusion criteria 
# read the US and CA125 orginal data
us_all <- readRDS(file.path(work_data, "US_alldates.RDS"))
ca125_all <- readRDS(file.path(work_data, "CA125_alldates.RDS"))

he_ova <- inc_excl(us_all = us_all, 
                   ca125_all = ca125_all, 
                   start_date="2013-04-01", 
                   end_date="2017-12-31", 
                   symptom=F) # no need symptoms

# Add demographics
# DOB data from CPRD
demo_cprd <- read_dta(file.path(proc_data, "demographics.dta"))
# Townsend 2011, 1 = least deprived
townsend <- read_dta(file.path(raw_link_data, "e_aurum_townsend2011_10.dta"))
# HES patient characteristics
hes_pat <- read_dta(file.path(raw_link_data, "e_aurum_hes_patient_21_001655_dm.dta"))
# CPRD ethnicity
ethn_cprd <- read_dta(file.path(proc_data, "ethnicity_cprd_alldates.dta"))

he_ova <- add_demo(he_ova = he_ova, 
                   demo_cprd = demo_cprd, 
                   townsend = townsend, 
                   ethn_cprd = ethn_cprd, 
                   hes_pat = hes_pat)

# Add death records

# linked death records
death <- read_dta(file.path(raw_link_data, "e_aurum_death_patient_21_001655_dm.dta"))
he_ova <- add_death(he_ova = he_ova, death = death)

rm(demo_cprd, townsend, hes_pat, ethn_cprd, us_all, ca125_all)

## Step 2:  ----
hes_admin_diag <- readRDS(file.path(work_data, "hes_admin_diag.rds"))

## Step 3:  ----

# CPRD cancer, from Kirsten
allcancers_cprd <- read_dta(file.path(proc_data, "allcancers.dta"))

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
saveRDS(he_ova, file.path(work_data, "he_ova_withoutSymptom.rds"))

# datasets for analyses ---------------------------------------------------
# for estimating costs associated with cancer and benign gynae diseases
# we do not need to limited to all inclusion/exclusion criteria
# so we use the preiouvs he_ova without implementing symptoms
# use the more comprehensive population, without necessarily linking to symptoms
# he_ova <- readRDS(file.path(work_data, "he_ova_withoutSymptom.rds"))

###########################################################
# make the analytical population consistent - 2024-08-01
he_ova <- readRDS(file.path(work_data, "he_ova_step4_withSymptom_prescr_NCRASupdated.RDS"))
###########################################################

# read in the hosptial admission records data for all
hes_admin_diag <- readRDS(file.path(work_data, "hes_admin_diag.rds"))

hes_hrg_cost <- readRDS(file.path(work_data, "hes_hrg_cost.rds"))

# hes_hrg_cost: cleaned HES HRG data with cost
# hes_admin_diag: HES admission dates
# he_ova: working dataset for HE analysis of ovarian and other cancers

## data for all ----
# this dataset include all study population, diagnosed with cancer or not
# data are split by year from index dates
# incidents are labelled as in same years, one years ago, .. over 4 years
cost_all <- all_hrg(hes_hrg_cost=hes_hrg_cost, 
                    hes_admin_diag=hes_admin_diag, 
                    he_ova=he_ova)

saveRDS(cost_all, file.path(work_data, "cost_all.rds"))

## dataset by cancer types----
cost_ova <- cancer_hrg(hes_hrg_cost=hes_hrg_cost, 
           hes_admin_diag=hes_admin_diag,                        
           he_ova=he_ova, 
           cancer_type="ova", 
           one_year_only = FALSE) 
saveRDS(cost_ova, file.path(work_data, "cost_ova.rds"))

cost_lung <- cancer_hrg(hes_hrg_cost=hes_hrg_cost, 
                       hes_admin_diag=hes_admin_diag, 
                       he_ova=he_ova, 
                       cancer_type="lung", 
                       one_year_only = FALSE)
saveRDS(cost_lung, file.path(work_data, "cost_lung.rds"))

cost_panc <- cancer_hrg(hes_hrg_cost=hes_hrg_cost, 
                       hes_admin_diag=hes_admin_diag, 
                       he_ova=he_ova, 
                       cancer_type="panc", 
                       one_year_only = FALSE)
saveRDS(cost_panc, file.path(work_data, "cost_panc.rds"))

cost_loGI <- cancer_hrg(hes_hrg_cost=hes_hrg_cost, 
                       hes_admin_diag=hes_admin_diag, 
                       he_ova=he_ova, 
                       cancer_type="loGI", 
                       one_year_only = FALSE)
saveRDS(cost_loGI, file.path(work_data, "cost_loGI.rds"))

cost_uter <- cancer_hrg(hes_hrg_cost=hes_hrg_cost, 
                        hes_admin_diag=hes_admin_diag, 
                        he_ova=he_ova, 
                        cancer_type="uter", 
                        one_year_only = FALSE)
saveRDS(cost_uter, file.path(work_data, "cost_uter.rds"))

cost_benign <- cancer_hrg(hes_hrg_cost=hes_hrg_cost, 
                        hes_admin_diag=hes_admin_diag, 
                        he_ova=he_ova, 
                        cancer_type="benign", 
                        one_year_only = TRUE)
# better just use those without cancer diagnosis within 1 year as a cohort
saveRDS(cost_benign, file.path(work_data, "cost_benign.rds"))

cost_other <- cancer_hrg(hes_hrg_cost=hes_hrg_cost, 
                        hes_admin_diag=hes_admin_diag, 
                        he_ova=he_ova, 
                        cancer_type="other", 
                        one_year_only = FALSE)
saveRDS(cost_other, file.path(work_data, "cost_other.rds"))

