library(tidyverse)
library(haven)
library(foreign)
library(parallel)
library(doSNOW)
library(dtplyr)
library(data.table)
library(survival)


# Step 1: create the study population -------------------------------------

## US + symptoms ----

# symptoms <- read_dta(file.path(proc_data, "observation_ovarysymptom.dta"))
# 
# symptoms <- symptoms %>%
#   lazy_dt() %>%
#   mutate(symptom_obs_date= as.Date(obsdate, "%d/%m/%Y"),
#          symptom_ent_date= as.Date(enterdate, "%d/%m/%Y")
#   ) %>%
#   mutate(symptom_obs_date=if_else(is.na(symptom_obs_date),symptom_ent_date, symptom_obs_date)) %>%
#   rename(epatid = e_patid) %>%
#   filter(symptom_obs_date >= "2013-01-01" & symptom_obs_date <= "2017-12-31") %>%
#   filter(symptom !=3) %>% # filter out "suspected cancer", because it is not ovary specific
#   select(epatid, symptom_obs_date) %>%
#   as.data.frame()
# 
# symptoms_wide <- symptoms %>%
#   group_by(epatid) %>%
#   mutate(symptom_obs_date_n = paste0("symptom_obs_date_", 1:n())) %>%
#   spread(symptom_obs_date_n, symptom_obs_date)
# 
# saveRDS(symptoms_wide, file.path(r_wd, "working_data", "symptoms_wide.RDS"))

# US + symptoms + prescription ----
# update from Kirsten
# 2024-04-24
# symptoms <- read_dta(file.path(proc_data, "ovary_sympt&prescr_alldates.dta"))
# 
# symptoms <- symptoms %>%
#   lazy_dt() %>%
#   mutate(symptom_obs_date = date) %>%
#   filter(symptom_obs_date >= "2013-01-01" & symptom_obs_date <= "2017-12-31") %>%
#   filter(ovarysymptom !=3) %>% # filter out "suspected cancer", because it is not ovary specific
#   select(epatid, symptom_obs_date) %>%
#   as.data.frame()
# 
# # there are a few ids with up to 2533 obs
# # too large after converting to wide form
# # split them into two
# symp1 <- symptoms %>% 
#   group_by(epatid) %>% 
#   mutate(n= n()) %>% 
#   filter(n<200) %>% 
#   select(-n)
#   
# symp2 <- symptoms %>% 
#   group_by(epatid) %>% 
#   mutate(n= n()) %>% 
#   filter(n>=200) %>% 
#   select(-n)
# 
# symptoms_wide1 <- symp1 %>%
#   group_by(epatid) %>%
#   mutate(symptom_obs_date_n = paste0("symptom_obs_date_", 1:n())) %>%
#   spread(symptom_obs_date_n, symptom_obs_date)
# 
# symptoms_wide2 <- symp2 %>%
#   group_by(epatid) %>%
#   mutate(symptom_obs_date_n = paste0("symptom_obs_date_", 1:n())) %>%
#   spread(symptom_obs_date_n, symptom_obs_date)
# 
# saveRDS(symptoms_wide1, file.path(r_wd, "working_data", "symptoms_prescr_wide1.RDS"))
# saveRDS(symptoms_wide2, file.path(r_wd, "working_data", "symptoms_prescr_wide2.RDS"))

# only include US with ovary symptoms reported with 3 months before
us_symptom_screen <- function(dt){
  
  # symptoms_wide <- readRDS(file.path(r_wd, "working_data", "symptoms_wide.RDS"))
  # use the updated symptom + prescription
  symptoms_wide1 <- readRDS(file.path(r_wd, "working_data", "symptoms_prescr_wide1.RDS"))
  symptoms_wide2 <- readRDS(file.path(r_wd, "working_data", "symptoms_prescr_wide2.RDS"))
  
  # spread long to wide for symptom records
  # wide1
  rt1 <- dt %>% 
    filter(!is.na(us_date)) %>% 
    inner_join(symptoms_wide1)%>% 
    lazy_dt() %>% 
    # calculate the gap between US date and symptom date
    # only keep ids with at least one gap within 3 months before US date
    filter(if_any(starts_with("symptom_obs_"), ~ as.numeric(us_date - .x)>0 & 
                    as.numeric(us_date - .x)<93)) %>% 
    select(-starts_with("symptom_obs_")) %>%
    as.data.frame()
  
  rm(symptoms_wide1)
  
  # wide2
  rt2 <- dt %>% 
    filter(!is.na(us_date)) %>% 
    inner_join(symptoms_wide2)%>% 
    lazy_dt() %>% 
    # calculate the gap between US date and symptom date
    # only keep ids with at least one gap within 3 months before US date
    filter(if_any(starts_with("symptom_obs_"), ~ as.numeric(us_date - .x)>0 & 
                    as.numeric(us_date - .x)<93)) %>% 
    select(-starts_with("symptom_obs_")) %>%
    as.data.frame()
  
  rm(symptoms_wide2)
  
  rt <- bind_rows(rt1, rt2)
  
  return(rt)
}

## inclusion criteria ----
inclusion <- function(us_all, ca125_all,
                      start_date="2013-04-01", 
                      end_date="2017-12-31", 
                      symptom=T){

  # 2013-04-01 - 2017-12-31
  # UScat==1 : diagnostic 
  # all records
  dt1 <- us_all %>% 
    filter(us_date >= start_date & us_date <= end_date & UScat==1) %>% 
    select(-UScat)
  
  dt2 <- ca125_all %>% 
    select(-studystart, -studyend, -studyperiod) %>% 
    filter(obsdate >= start_date & obsdate <= end_date) %>% 
    rename(ca125_date=obsdate)
  
  # eliminate procedure on the same day
  # or procedure within 2013-04-01 - 2017-12-31 and 1 year before diagnostic US/CA125
  # this is a wide form list of procedural US
  dt11 <- us_all %>% 
    filter(us_date >= start_date & us_date <= end_date & UScat==6) %>% 
    rename(us_date_proc = us_date) %>% 
    select(epatid, us_date_proc) %>% 
    group_by(epatid) %>% 
    mutate(us_date_proc_n = paste0("us_date_proc_", 1:n())) %>% 
    spread(us_date_proc_n, us_date_proc)
  
  # tested, add as.numeric for date or not, the results are the same
  dt1 <- dt1 %>% 
    left_join(dt11) %>% 
    mutate(excl = case_when(
      if_any(starts_with("us_date_proc_"), ~ us_date - .x <= 365.25 & us_date - .x >= 0) ~ 1,
      .default = 0
    )) %>% 
    filter(excl==0) %>% 
    select(-starts_with("us_date_proc_"), -excl) %>% 
    distinct()  
  
  dt2 <- dt2 %>% # also 
    left_join(dt11) %>% 
    mutate(excl = case_when(
      if_any(starts_with("us_date_proc_"), ~ ca125_date - .x <= 365.25 & ca125_date - .x >= 0) ~ 1,
      .default = 0
    )) %>% 
    filter(excl==0) %>% 
    select(-starts_with("us_date_proc_"), -excl) %>% 
    distinct()  
  
  # use ovary symptom to exclude patients' US records
  # after PMG meeting 2023-10
  if (symptom) {
    dt1 <- us_symptom_screen(dt1)
  }

  # generate a dataset only include the earliest diagnostic US/CA125 records
  # US
  dt1 <- dt1 %>% 
    group_by(epatid) %>% 
    mutate(us_date_earliest = min(us_date)) %>% 
    filter(us_date == us_date_earliest) %>% 
    distinct() %>% 
    select(-us_date_earliest)
  # CA125
  dt2 <- dt2 %>% 
    group_by(epatid) %>% 
    mutate(ca125_date_earliest = min(ca125_date)) %>% 
    filter(ca125_date == ca125_date_earliest) %>% 
    distinct() %>% 
    select(-ca125_date_earliest)
  
  # combine dt1 and dt2
  # keep all records
  he_ova <- merge(dt1, dt2, all = T)
  
  # generate index date as the earlier one if have both US and CA125
  # Gartha and Kirsten decided CA125=0 is not reliable
  he_ova <- he_ova %>% 
    mutate(ca125_date = if_else(CA125==0, NA, ca125_date)) %>% 
    mutate(CA125 = if_else(CA125==0, NA, CA125)) %>% 
    filter(!is.na(ca125_date) | !is.na(us_date)) %>% 
    rowwise() %>% 
    mutate(index_date = min(us_date, ca125_date, na.rm = T)) 
  
  return(he_ova)
}


## exclusion criteria----
# enter the return of the inclusion

exclusion <- function(us_all, ca125_all, he_ova, 
                      start_date="2013-04-01"){

  # generate datasets for exclusion
  # 1. US in 2012-04-01 to 2013-03-31
  # 2. CA125 in 2012-04-01 to 2013-03-31
  t1 <- as.Date(start_date) - 365
  t2 <- as.Date(start_date) - 1
  
  # only need the latest records
  ex1 <- us_all %>% 
    lazy_dt() %>% 
    filter(us_date >= t1 & us_date <= t2) %>% 
    group_by(epatid) %>% 
    mutate(us_date_latest = max(us_date)) %>%
    filter(us_date == us_date_latest) %>% 
    distinct() %>% 
    select(-doppler, -us_date_latest, -UScat) %>% 
    rename(us_date_ex=us_date) %>% 
    as.data.frame()
  
  ex2 <- ca125_all %>% 
    lazy_dt() %>% 
    filter(obsdate >= t1 & obsdate <= t2) %>%
    group_by(epatid) %>% 
    mutate(ca125_date_latest = max(obsdate)) %>%
    filter(obsdate == ca125_date_latest) %>% 
    distinct() %>% 
    select(-ca125_date_latest, -studystart, -studyend, -studyperiod, -CA125, -units) %>% 
    rename(ca125_date_ex=obsdate)%>% 
    as.data.frame()
  
  # combine ex1 and ex2
  # keep all records
  ex <- merge(ex1, ex2, all = T)
  
  # implement exclusion criteria 
  ex_list <- he_ova %>% 
    left_join(ex) %>% 
    mutate(exclusion= case_when(as.numeric(index_date - ca125_date_ex) <= 365.25 ~ 1, 
                                as.numeric(index_date - us_date_ex) <= 365.25 ~ 1,
                                .default = 0
    )) %>% 
    filter(exclusion==1) %>% 
    select(epatid) %>% 
    distinct()
  # exclude these ids
  he_ova <- he_ova %>% 
    filter(!epatid %in% ex_list$epatid)
  
  return(he_ova)
}

## combine inclusion and exclusion ----

inc_excl <- function(us_all, ca125_all, 
                     start_date="2013-04-01", 
                     end_date="2017-12-31", 
                     symptom=T){
  
  he_ova <- inclusion(us_all = us_all, 
                      ca125_all = ca125_all, 
                      start_date = start_date, 
                      end_date = end_date, 
                      symptom=symptom)
  
  he_ova <- exclusion(us_all = us_all, 
                      ca125_all = ca125_all, 
                      he_ova = he_ova, 
                      start_date = start_date)
  
  return(he_ova)
  
}


## Add demographics ----

add_demo <- function(he_ova, demo_cprd, townsend, ethn_cprd, hes_pat){

  # Townsend 2011, 1 = least deprived
  townsend <- townsend %>% 
    lazy_dt() %>% 
    rename(epatid = e_patid) %>% 
    # impute missing townsend score by the majority score in the practice
    group_by(e_pracid) %>% 
    add_count(e2011_townsend_10) %>% 
    mutate(townsend_10_imp = if_else(e2011_townsend_10=="", e2011_townsend_10[which.max(n)], 
                                     e2011_townsend_10)) %>% 
    ungroup() %>% 
    select(-e2011_townsend_10, -n) %>% 
    as.data.frame()
  
  # HES patient characteristics
  hes_pat <- hes_pat %>% 
    lazy_dt() %>% 
    select(e_patid, gen_ethnicity) %>%  # e_pracid is uncomplete from hes_pat, so use that from townsend
    rename(epatid = e_patid) %>% 
    as.data.frame()
  
  # merge into he_ova
  he_ova <- he_ova %>% 
    left_join(demo_cprd) %>% 
    select(-female, -yob) %>% 
    left_join(townsend) %>% 
    left_join(hes_pat) 
  
  he_ova <- he_ova %>% 
    lazy_dt() %>% 
    mutate(town5 = case_when(
      townsend_10_imp %in% c("1", "2") ~ "Q1", # least deprived
      townsend_10_imp %in% c("3", "4") ~ "Q2",
      townsend_10_imp %in% c("5", "6") ~ "Q3",
      townsend_10_imp %in% c("7", "8") ~ "Q4",
      townsend_10_imp %in% c("9", "10") ~ "Q5"
    )) %>% 
    mutate(town5 = factor(town5, levels= c("Q1","Q2","Q3","Q4","Q5"))) %>% 
    as.data.frame()
  
  he_ova <- he_ova %>% 
    mutate(ethn = case_when(
      gen_ethnicity %in% c("Bl_Afric", "Bl_Carib", "Bl_Other") ~ "Black",
      gen_ethnicity %in% c("Indian", "Pakistani", "Bangladesi") ~ "South Asian",
      gen_ethnicity %in% c("Chinese", "Oth_Asian") ~ "Other Asian",
      gen_ethnicity == "White" ~ "White",
      gen_ethnicity %in% c("Mixed", "Other") ~ "Mixed or others",
      gen_ethnicity %in% c("", "Unknown") | is.na(gen_ethnicity) ~ NA
    )) %>% 
    mutate(ethn = factor(ethn, levels = c("White", "Black", "South Asian", "Other Asian", "Mixed or others"))) %>% 
    as.data.frame()
  
  # CPRD ethnicity
  ethn_cprd <- ethn_cprd %>% 
    lazy_dt() %>% 
    mutate(ethn_cprd = case_when(ethnicity==1 ~ "Asian or British Asian",
                                 ethnicity==2 ~ "Black or Black British",
                                 ethnicity==3 ~ "Mixed",
                                 ethnicity==4 ~ "White",
                                 ethnicity==5 ~ "Other"
    )) %>% 
    # because it has multiple records for one id, and some id has more than one ethn 
    # keep the ethn that appears in most times
    group_by(epatid) %>% 
    add_count(ethn_cprd) %>% 
    mutate(ethn_cprd_max = ethn_cprd[which.max(n)]) %>% 
    ungroup() %>%
    select(epatid, ethn_cprd_max) %>% 
    rename(ethn_cprd=ethn_cprd_max) %>% 
    distinct() %>% 
    as.data.frame()
  
  
  he_ova <- he_ova %>% 
    left_join(ethn_cprd) %>% 
    mutate(ethn_comb = case_when(ethn=="White" | ethn_cprd=="White" ~ "White",
                                 ethn %in% c("South Asian", "Other Asian") | ethn_cprd == "Asian or British Asian" ~ "Asian",
                                 ethn=="Black" | ethn_cprd=="Black or Black British" ~ "Black",
                                 ethn=="Mixed or others" | ethn_cprd %in% c("Mixed", "Other") ~ "Mixed or others",
                                 .default = NA
    )) %>% 
    mutate(ethn_comb = factor(ethn_comb, levels = c("White", "Asian", "Black", "Mixed or others"))) 
  
  
  he_ova <- he_ova %>% 
    mutate(ethn_comb_imp = if_else(is.na(ethn_comb), "White", ethn_comb)) %>% 
    mutate(ethn_comb_imp = factor(ethn_comb_imp, levels = c("White", "Asian", "Black", "Mixed or others"))) 
  
  # generate age at index date
  he_ova$age <- round(as.numeric(he_ova$index_date - he_ova$dobdate)/365.25)
  # keep aged >=18 at index date
  he_ova <- subset(he_ova, age>=18)
  
  # create age group
  he_ova <- he_ova %>% 
    mutate(age_group10 = case_when(
      age<20 ~ "< 20",
      age>=20 & age<30 ~ "20-29",
      age>=30 & age<40 ~ "30-39",
      age>=40 & age<50 ~ "40-49",
      age>=50 & age<60 ~ "50-59",
      age>=60 & age<70 ~ "60-69",
      age>=70 & age<80 ~ "70-79",
      age>=80 & age<90 ~ "80-89",
      age>=90 ~ "90+",
    ))
  
  # according to the tabulate results,
  # for ova, use group <40, 40-49, 50-59, 60-69, 70-79, 80+
  # for uter and loGI, use group <50, 50-59, 60-69, 70-79, 80+
  # for lung and panc, use group <60, 60-69, 70-79, 80+ 
  he_ova <- he_ova %>% 
    mutate(age_group_40_80 = case_when(
      age<40 ~ "<40",
      age>=40 & age<50 ~ "40-49",
      age>=50 & age<60 ~ "50-59",
      age>=60 & age<70 ~ "60-69",
      age>=70 & age<80 ~ "70-79",
      age>=80 ~ ">=80",
    ),
    age_group_50_80 = case_when(
      age<50 ~ "<50",
      age>=50 & age<60 ~ "50-59",
      age>=60 & age<70 ~ "60-69",
      age>=70 & age<80 ~ "70-79",
      age>=80 ~ ">=80",
    ),
    age_group_60_80 = case_when(
      age<60 ~ "<60",
      age>=60 & age<70 ~ "60-69",
      age>=70 & age<80 ~ "70-79",
      age>=80 ~ ">=80",
    ),
    age_group_40_80=factor(age_group_40_80, levels=c("<40","40-49","50-59","60-69","70-79",">=80")),
    age_group_50_80=factor(age_group_50_80, levels=c("<50","50-59","60-69","70-79",">=80")),
    age_group_60_80=factor(age_group_60_80, levels=c("<60","60-69","70-79",">=80"))
    )
  
  return(he_ova)
  
}

## Add death ----
add_death <- function(he_ova, death){
  
  # linked death records
  death <- read_dta(file.path(raw_link_data, "e_aurum_death_patient_21_001655_dm.dta"))
  
  # code if a cancer type is the cause of death
  death$death_cancer_ova <- 0
  death$death_cancer_panc <- 0
  death$death_cancer_lung <- 0
  death$death_cancer_uter <- 0
  death$death_cancer_loGI <- 0
  death$death_cancer <- 0
  
  for (i in c("", 1:15)) {
    text <- paste0("cause", i) 
    death$death_cancer_ova[grepl("^C56|^C57.0|^C48.1|^C48.2", death[[text]])] <- 1 # exclude D39.1 as borderline tumour
    death$death_cancer_panc[grepl("^C25", death[[text]])] <- 1
    death$death_cancer_lung[grepl("^C34", death[[text]])] <- 1
    death$death_cancer_uter[grepl("^C54|^C55", death[[text]])] <- 1
    death$death_cancer_loGI[grepl("^C18|^C19|^C20|^C21", death[[text]])] <- 1
    death$death_cancer[grepl("^C", death[[text]]) & !grepl("^C44", death[[text]])] <- 1
  }
  
  death <- death %>% 
    lazy_dt() %>% 
    mutate(dor = as.Date(dor, "%d/%m/%Y"), 
           dod = as.Date(dod, "%d/%m/%Y")) %>% 
    mutate(death_date = if_else(is.na(dod), dor, dod)) %>% 
    rename(epatid = e_patid) %>% 
    select(epatid, death_date, 
           death_cancer_ova, death_cancer_panc, death_cancer_lung, 
           death_cancer_uter, death_cancer_loGI, death_cancer) %>% 
    as.data.frame()
  
  he_ova <- he_ova %>% 
    left_join(death)
  
  # 5 died before index date and 5 died on the same day
  # check!
  # there is no error in data preparation
  # xx <- he_ova[as.numeric(he_ova$death_date - he_ova$index_date)<0 & !is.na(he_ova$death_date), 
  #        c("epatid", "us_date", "ca125_date", "index_date", "death_date")]
  # xxx <- xx$epatid
  # us_all[us_all$epatid %in% xxx, ]
  # they do come from CPRD/HES DID
  # remove the 5 who died before index date
  he_ova <- he_ova %>% 
    filter(is.na(death_date) | as.numeric(death_date - index_date>=0))
  
  return(he_ova)
}


# Step 2: admission & diagnosis -------------------------------------------
gen_admin_diag <- function(hes_hosp, hes_hosp_diag, hes_epis, icd_codes){
  
  # cancer types
  ova <- icd_codes["ova"]
  panc <- icd_codes["panc"]
  lung <- icd_codes["lung"]
  uter <- icd_codes["uter"]
  loGI <- icd_codes["loGI"]
  benign <- icd_codes["benign"]
  borderline <- icd_codes["borderline"]

  # generate a set of cancer types
  # by unique patid and spell id
  # long to wide
  hes_hosp_cancer <- hes_hosp_diag %>%
    lazy_dt() %>% 
    # ICDs here are with a single decimal or no decimal
    mutate(cancer_ova = if_else(grepl(ova, icd), 1, 0),
           cancer_panc = if_else(grepl(panc, icd), 1, 0),
           cancer_lung = if_else(grepl(lung, icd), 1, 0),
           cancer_uter = if_else(grepl(uter, icd), 1, 0), # D39.0 removed
           cancer_loGI = if_else(grepl(loGI, icd), 1, 0), #^D37.3|^D37.4|^D37.5 removed
           cancer = if_else(grepl("^C", icd) & !grepl("^C44", icd), 1, 0), #^D37|^D38|^D39|^D4 removed
           benign_gynae = if_else(grepl(benign, icd), 1, 0),
           borderline = if_else(grepl(borderline, icd), 1, 0)
    ) %>%
    group_by(e_patid, spno) %>% # within an admission, as long as one record is 1
    summarise(cancer_ova = max(cancer_ova),
              cancer_panc = max(cancer_panc),
              cancer_lung = max(cancer_lung),
              cancer_uter = max(cancer_uter),
              cancer_loGI = max(cancer_loGI),
              cancer = max(cancer),
              benign_gynae = max(benign_gynae),
              borderline = max(borderline)
    ) %>% 
    as.data.frame()
  
  # merge
  # the former include all latter spno id
  hes_hosp <- hes_hosp %>%
    mutate(admidate=as.Date(admidate, "%d/%m/%Y"),
           discharged=as.Date(discharged, "%d/%m/%Y"),
           elecdate=as.Date(elecdate, "%d/%m/%Y")) %>%
    left_join(hes_hosp_cancer, by=c("e_patid", "spno")) %>%
    rename(epatid = e_patid)
  
  # impute missing admission dates
  hes_epis <- hes_epis %>%
    select(e_patid, spno, epikey, admidate, epistart, epiend) %>%
    rename(epatid=e_patid, admidate2=admidate)
  
  hosp_to_fill <- hes_hosp %>%
    filter(is.na(admidate)) %>%
    select(epatid, spno) %>%
    left_join(hes_epis) %>%
    mutate(epistart_imp = if_else(epistart=="", epiend, epistart)) %>%
    mutate(epistart_imp=as.Date(epistart_imp, "%d/%m/%Y")) %>%
    group_by(epatid, spno) %>%
    mutate(epistart_imp2 = min(epistart_imp)) %>%
    filter(epistart_imp==epistart_imp2) %>%
    select(epatid, spno, epistart_imp) %>%
    distinct()
  
  hes_hosp <- hes_hosp %>%
    left_join(hosp_to_fill) %>%
    mutate(admidate = if_else(is.na(admidate), epistart_imp, admidate)) %>%
    select(-epistart_imp)
  
  # hes_hosp has 6189 discharge dates are missing
  # use episode end date to fill 
  hosp_to_fill2 <- hes_hosp %>%
    filter(is.na(discharged)) %>%
    select(epatid, spno) %>%
    left_join(hes_epis) %>%
    mutate(epiend_imp = if_else(epiend=="", epistart, epiend)) %>%
    mutate(epiend_imp=as.Date(epiend_imp, "%d/%m/%Y")) %>%
    group_by(epatid, spno) %>%
    mutate(epiend_imp2 = max(epiend_imp)) %>%
    filter(epiend_imp==epiend_imp2) %>%
    select(epatid, spno, epiend_imp) %>%
    distinct()
  
  hes_hosp <- hes_hosp %>%
    left_join(hosp_to_fill2) %>%
    mutate(discharged = if_else(is.na(discharged), epiend_imp, discharged)) %>%
    select(-epiend_imp)
  
  # code those discharged date earlier than admission date as the same day
  # because the median difference between admission and discharge is 0
  # 53 are discharge < admission
  hes_hosp$discharged[as.numeric(hes_hosp$discharged-hes_hosp$admidate)<0] <- 
    hes_hosp$admidate[as.numeric(hes_hosp$discharged-hes_hosp$admidate)<0]
  
  # summary(as.numeric(hes_hosp$discharged-hes_hosp$admidate))
  # there are a few super long stay, due to very early admission dates
  # the issue could be addressed after select a later starting point, such as 2013
  
  return(hes_hosp)
}


# Step 3: add incident events and code history ----------------------------

## exclude those with ovarian cancer + borderline before index dates----
# use both HES and CPRD data
exclude_early_oc <- function(he_ova, hes_admin_diag, allcancers_cprd){
  
  hes_hosp2 <- hes_admin_diag %>%
    filter(epatid %in% he_ova$epatid)
  
  # HES records
  ova_date_excl <- hes_hosp2 %>% 
    lazy_dt() %>% 
    filter(cancer_ova==1 | borderline==1) %>% 
    group_by(epatid) %>% 
    mutate(ova_date = min(admidate)) %>%
    filter(admidate == ova_date) %>% 
    select("epatid", "ova_date") %>% 
    distinct() %>% 
    as.data.frame()

  ova_date_excl_cprd <- allcancers_cprd %>% 
    filter(epatid %in% he_ova$epatid) %>% 
    filter(cancerdate<"2024-01-01") %>% 
    mutate(cancer_ova=case_when(site==18 ~ 1,# include borderline ovary tumour
                                .default = 0)) %>% 
    filter(cancer_ova==1) %>% 
    group_by(epatid) %>% 
    mutate(ova_date_cprd = min(cancerdate)) %>%
    filter(cancerdate == ova_date_cprd) %>% 
    select("epatid", "ova_date_cprd") %>% 
    distinct()
  
  he_ova <- he_ova %>% 
    left_join(ova_date_excl) %>% 
    left_join(ova_date_excl_cprd) %>% 
    mutate(prior_oc=case_when(ova_date <= index_date ~ 1, 
                              ova_date_cprd <= index_date ~ 1, 
                              .default=0)) %>% 
    filter(prior_oc == 0) %>% 
    select(-ova_date, -ova_date_cprd, -prior_oc) 
  
  return(he_ova)
}

## earliest diagnosis dates by cancer types----
# cancer_type coulde be: ova, lung, panc, loGI, uter and cancer (for all)
# also can be benign
earliest_date <- function(hes_hosp3, cancer_type){
  
  if (cancer_type=="cancer" | cancer_type=="borderline")
    cancer_name <- cancer_type else
      if (cancer_type=="benign")
        cancer_name <- "benign_gynae" else
          cancer_name <- paste0("cancer_", cancer_type)
        
        c_date <- paste0(cancer_type, "_date")
        
        rt <- hes_hosp3 %>% 
          filter(get(cancer_name)==1) %>% 
          group_by(epatid) %>% 
          mutate(!!c_date := min(admidate)) %>%
          filter(admidate == get(c_date)) 
        
        if (cancer_type=="cancer"){
          rt <- rt %>% 
            select("epatid", c_date) %>% 
            distinct()
        } else {
          rt <- rt %>% 
            select("epatid","spno", c_date) %>% 
            distinct() 
        }
        
        return(rt)
}

# merge into a dataset of all admission records for
# those in the study population
# no one have ovarian cancer records before index date
# each one has a earliest date for all types of cancer 
# if one has one type of specified cancers, a earliest date is recorded 
# hes_hosp4 <- hes_hosp3 %>% 
#   left_join(cancer_date, by=c("epatid")) %>% 
#   # a few admissions with the earliest cancer diagnosis on the same day
#   # here we keep all of them
#   # so one person can have more than one admission with the same earliest diagnosis date for a cancer, e.g. ovarian
#   left_join(ova_date, by=c("epatid", "spno")) %>% 
#   left_join(panc_date, by=c("epatid", "spno")) %>% 
#   left_join(lung_date, by=c("epatid", "spno")) %>% 
#   left_join(uter_date, by=c("epatid", "spno")) %>% 
#   left_join(loGI_date, by=c("epatid", "spno")) %>% 
#   left_join(benign, by=c("epatid", "spno"))
# 
# saveRDS(hes_hosp4, file.path(r_wd, "working_data", "hes_hosp_cancer_date.rds"))

## Patient level cancer and date records ----- 

pat_diag <- function(he_ova, hes_admin_diag){
  
  # hospital admission records, 
  # with missing admission dates imputed
  # only include those in the study population
  # only include those with the admission date later than the index date
  # only include those with cancer or benign gynae
  hes_hosp3 <- hes_admin_diag %>% 
    ## remove admission dates before index dates----
  left_join(he_ova[, c("epatid", "index_date")]) %>% 
    filter(epatid %in% he_ova$epatid) %>% 
    filter(admidate >= index_date) %>% 
    filter(cancer==1 | benign_gynae==1 | cancer_ova==1) 
  
  # earliest diagnosis dates by cancer types
  ova_date <- earliest_date(hes_hosp3, "ova")
  panc_date <- earliest_date(hes_hosp3, "panc")
  lung_date <- earliest_date(hes_hosp3, "lung")
  uter_date <- earliest_date(hes_hosp3, "uter")
  loGI_date <- earliest_date(hes_hosp3, "loGI")
  cancer_date <- earliest_date(hes_hosp3, "cancer")
  benign <- earliest_date(hes_hosp3, "benign")
  borderline <- earliest_date(hes_hosp3, "borderline")
  
  # Patient level cancer and date records
  he_ova <- he_ova %>% 
    left_join(cancer_date, by=c("epatid")) %>% 
    # a few admissions with the earliest cancer diagnosis on the same day
    # only keep one record here
    left_join(ova_date[!duplicated(ova_date$epatid), c("epatid", "ova_date")], by="epatid") %>% 
    left_join(panc_date[!duplicated(panc_date$epatid), c("epatid", "panc_date")], by="epatid") %>% 
    left_join(lung_date[!duplicated(lung_date$epatid), c("epatid", "lung_date")], by="epatid") %>% 
    left_join(uter_date[!duplicated(uter_date$epatid), c("epatid", "uter_date")], by="epatid") %>% 
    left_join(loGI_date[!duplicated(loGI_date$epatid), c("epatid", "loGI_date")], by="epatid") %>% 
    left_join(benign[!duplicated(benign$epatid), c("epatid", "benign_date")], by=c("epatid")) %>% 
    left_join(borderline[!duplicated(borderline$epatid), c("epatid", "borderline_date")], by=c("epatid")) %>% 
    mutate(cancer=if_else(is.na(cancer_date), 0, 1),
           cancer_ova=if_else(is.na(ova_date), 0, 1),
           cancer_panc=if_else(is.na(panc_date), 0, 1),
           cancer_lung=if_else(is.na(lung_date), 0, 1),
           cancer_uter=if_else(is.na(uter_date), 0, 1),
           cancer_loGI=if_else(is.na(loGI_date), 0, 1),
           benign_gynae=if_else(is.na(benign_date), 0, 1),
           borderline=if_else(is.na(borderline_date), 0, 1)
    )
  
  return(he_ova)
}

## baseline cancer records----
# identify patient with non-ovarian cancer or benign gynae diagnosis before index date
add_history <- function(he_ova, hes_admin_diag, allcancers_cprd){
  
  baseline_disease <- hes_admin_diag %>% 
    # join index date
    left_join(he_ova[, c("epatid", "index_date")]) %>% 
    filter(epatid %in% he_ova$epatid) %>% 
    # before index date
    filter(admidate < index_date) %>% 
    # baseline non-ovarian cancer history
    # baseline benign gynae
    # also add specific cancer types 
    mutate(b_cancer_nonOC = if_else(cancer==1 & cancer_ova!=1, 1, 0),
           b_benign = if_else(benign_gynae==1, 1, 0), 
           b_cancer_panc = if_else(cancer_panc==1, 1, 0), 
           b_cancer_lung = if_else(cancer_lung==1, 1, 0),
           b_cancer_uter = if_else(cancer_uter==1, 1, 0),
           b_cancer_loGI = if_else(cancer_loGI==1, 1, 0)
    ) %>% 
    select(epatid, b_cancer_nonOC, b_benign, b_cancer_panc, b_cancer_lung, 
           b_cancer_uter, b_cancer_loGI) %>% 
    group_by(epatid) %>% 
    mutate(b_cancer_nonOC=max(b_cancer_nonOC),  # any history = 1
           b_benign=max(b_benign), 
           b_cancer_panc=max(b_cancer_panc),
           b_cancer_lung=max(b_cancer_lung),
           b_cancer_uter=max(b_cancer_uter),
           b_cancer_loGI=max(b_cancer_loGI)
    ) %>% 
    ungroup() %>% 
    filter(if_any(starts_with("b_"), ~ . ==1)) %>% 
    distinct()
  
  baseline_disease_cprd <- allcancers_cprd %>% 
    filter(epatid %in% he_ova$epatid) %>% 
    filter(cancerdate<"2024-01-01") %>% # filter out one records = 2999
    left_join(he_ova[, c("epatid", "index_date")]) %>% 
    # before index date
    filter(cancerdate < index_date) %>% 
    mutate(b_cancer_nonOC_cprd=case_when(site!=18 ~ 1,# include borderline ovary tumour
                                         .default = 0),
           b_cancer_panc_cprd = if_else(site==19, 1, 0), 
           b_cancer_lung_cprd = if_else(site==11, 1, 0),
           b_cancer_uter_cprd = if_else(site==23, 1, 0),
           b_cancer_loGI_cprd = if_else(site==6, 1, 0)# colorectal
    ) %>% 
    select(epatid, b_cancer_nonOC_cprd, b_cancer_panc_cprd, b_cancer_lung_cprd, 
           b_cancer_uter_cprd, b_cancer_loGI_cprd) %>% 
    group_by(epatid) %>% 
    mutate(b_cancer_nonOC_cprd=max(b_cancer_nonOC_cprd),
           b_cancer_panc_cprd=max(b_cancer_panc_cprd),
           b_cancer_lung_cprd=max(b_cancer_lung_cprd),
           b_cancer_uter_cprd=max(b_cancer_uter_cprd),
           b_cancer_loGI_cprd=max(b_cancer_loGI_cprd)
    ) %>% 
    ungroup() %>% 
    filter(if_any(starts_with("b_"), ~ . ==1)) %>% 
    distinct()
  
  # merge back to the patient-level cancer record
  he_ova <- he_ova %>% 
    left_join(baseline_disease) %>% 
    left_join(baseline_disease_cprd) %>% 
    mutate(b_cancer_nonOva = case_when(b_cancer_nonOC==1 ~ 1,
                                       b_cancer_nonOC_cprd==1 ~ 1,
                                       .default = 0),
           b_benign_gynae = if_else(is.na(b_benign), 0, b_benign),
           b_panc = case_when(b_cancer_panc==1 ~ 1,
                              b_cancer_panc_cprd==1 ~ 1,
                              .default = 0),
           b_lung = case_when(b_cancer_lung==1 ~ 1,
                              b_cancer_lung_cprd==1 ~ 1,
                              .default = 0),
           b_uter = case_when(b_cancer_uter==1 ~ 1,
                              b_cancer_uter_cprd==1 ~ 1,
                              .default = 0),
           b_loGI = case_when(b_cancer_loGI==1 ~ 1,
                              b_cancer_loGI_cprd==1 ~ 1,
                              .default = 0)
    ) %>% 
    select(-b_cancer_nonOC, -b_cancer_nonOC_cprd, -b_benign, -b_cancer_panc,
           -b_cancer_lung, -b_cancer_uter, -b_cancer_loGI, -b_cancer_panc_cprd,
           -b_cancer_lung_cprd, -b_cancer_uter_cprd, -b_cancer_loGI_cprd)
  
  return(he_ova)
}


## cancer within 1 year after index date ----
add_incident <- function(he_ova){
  
  # 1-year incident
  he_ova <- he_ova %>% 
    mutate(ova_1yr_date = if_else(!is.na(ova_date) & 
                                    as.numeric(ova_date - index_date) <=365.25, 
                                  ova_date, NA),
           cancer_ova_1yr = if_else(!is.na(ova_1yr_date), 1, 0),
           
           panc_1yr_date = if_else(!is.na(panc_date) & 
                                     as.numeric(panc_date - index_date) <=365.25 & 
                                     b_panc==0, panc_date, NA),
           cancer_panc_1yr = if_else(!is.na(panc_1yr_date), 1, 0),
           
           lung_1yr_date = if_else(!is.na(lung_date) & 
                                     as.numeric(lung_date - index_date) <=365.25 &
                                     b_lung==0, lung_date, NA),
           cancer_lung_1yr = if_else(!is.na(lung_1yr_date), 1, 0),
           
           uter_1yr_date = if_else(!is.na(uter_date) & 
                                     as.numeric(uter_date - index_date) <=365.25 &
                                     b_uter==0, uter_date, NA),
           cancer_uter_1yr = if_else(!is.na(uter_1yr_date), 1, 0),
           
           loGI_1yr_date = if_else(!is.na(loGI_date) & 
                                     as.numeric(loGI_date - index_date) <=365.25 &
                                     b_loGI==0, loGI_date, NA),
           cancer_loGI_1yr = if_else(!is.na(loGI_1yr_date), 1, 0),
           
           # any (non-ova) cancer in 1yr excludes any (non-ova) cancer history
           # however, x-type cancer in 1yr only excludes x-type cancer history
           # therefore, any cancer in 1yr does not include all x-type cancer in 1yr
           cancer_1yr_date = if_else(!is.na(cancer_date) & 
                                       as.numeric(cancer_date - index_date) <=365.25 &
                                       b_cancer_nonOva==0, cancer_date, NA),
           cancer_1yr = if_else(!is.na(cancer_1yr_date), 1, 0),
           
           benign_1yr_date = if_else(!is.na(benign_date) & 
                                       as.numeric(benign_date - index_date) <=365.25, 
                                     benign_date, NA),
           benign_1yr = if_else(!is.na(benign_1yr_date), 1, 0)
    )
  
  # all incident
  he_ova <- he_ova %>% 
    mutate(ova_icd_date = ova_date,
           cancer_ova_icd = if_else(!is.na(ova_icd_date), 1, 0),
           
           panc_icd_date = if_else(!is.na(panc_date) & 
                                     b_panc==0, panc_date, NA),
           cancer_panc_icd = if_else(!is.na(panc_icd_date), 1, 0),
           
           lung_icd_date = if_else(!is.na(lung_date) & 
                                     b_lung==0, lung_date, NA),
           cancer_lung_icd = if_else(!is.na(lung_icd_date), 1, 0),
           
           uter_icd_date = if_else(!is.na(uter_date) & 
                                     b_uter==0, uter_date, NA),
           cancer_uter_icd = if_else(!is.na(uter_icd_date), 1, 0),
           
           loGI_icd_date = if_else(!is.na(loGI_date) & 
                                     b_loGI==0, loGI_date, NA),
           cancer_loGI_icd = if_else(!is.na(loGI_icd_date), 1, 0),
           
           cancer_icd_date = if_else(!is.na(cancer_date) &
                                       b_cancer_nonOva==0, cancer_date, NA),
           cancer_icd = if_else(!is.na(cancer_icd_date), 1, 0),
           
           benign_icd_date = benign_date,
           benign_icd = if_else(!is.na(benign_icd_date), 1, 0)
    )
  
  he_ova <- he_ova %>% 
    mutate(death = if_else(is.na(death_date), 0, 1)) %>% 
    relocate(death, .before = death_date) %>% 
    mutate(censor_date = if_else(is.na(death_date), as.Date("2021-03-31"), death_date)) %>% 
    relocate(censor_date, .after = death_date) %>% 
    mutate(ova_fu_time = as.numeric(censor_date - ova_1yr_date)) %>% 
    mutate(ova_fu_time = if_else(ova_fu_time==0, ova_fu_time + 1, ova_fu_time)) %>% # avoid 0 fu time
    mutate(ova_diag_time = as.numeric(ova_1yr_date - index_date)) 
  
  return(he_ova)
}



# Table one ---------------------------------------------------------------

ova_tableone <- function(he_ova, adj_pop_tag){
  
  library(tableone)
  
  vars <- c("age", "CA125", "town5", "ethn4", "b_cancer_nonOva", "b_benign_gynae",
            "cancer_ova_1yr", "cancer_panc_1yr", "cancer_lung_1yr",
            "cancer_uter_1yr", "cancer_loGI_1yr", "cancer_other_1yr", "benign_1yr",
            "ova_diag_time", "fu_time"
  )
  
  cat <- c("town5", "ethn4", "b_cancer_nonOva", "b_benign_gynae",
           "cancer_ova_1yr", "cancer_panc_1yr", "cancer_lung_1yr",
           "cancer_uter_1yr", "cancer_loGI_1yr", "cancer_other_1yr", "benign_1yr")
  
  tab <- CreateTableOne(vars = vars, data = he_ova, factorVars = cat)

  rt <- print(tab, digits = 2)
  
  write.csv(rt, file.path(output, paste0("tb", adj_pop_tag, ".csv")))
  
  return(tab)
  
}

tableone2 <- function(d4s){
  
  library(tableone)
  
  d4s$age <- d4s$age_cent60*10+60
  
  vars <- c("age", "town5", "ethn3", "cancer_type2", "any_late", "fu_time", 
            "death_cancer", "death_other")
  cat <- c("town5", "ethn3", "cancer_type2", "any_late","death_cancer", "death_other")
  
  tab <- CreateTableOne(vars = vars, data = d4s, factorVars = cat)
  
  write.csv(print(tab), file.path(output, "summary_4risk_models.csv"))
  
  return(tab)
}


stage_tab <- function(he_ova, cancer_type){
  
  tab <- table(he_ova[[paste0("stage_", cancer_type, "_imp")]], useNA = "ifany")
  
  tot <- table(he_ova[[paste0("cancer_", cancer_type, "_1yr")]], useNA = "ifany")["1"]
  
  c <- tab["stage 1"] + tab["stage 2"]
  p <- round(c/tot, 3)
  s12 <- c(c, p)
  
  c <- tab["stage 3"] + tab["stage 4"]
  p <- round(c/tot, 3)
  s34 <- c(c, p)
  
  c <- tab[is.na(names(tab))]
  p <- round(c/tot, 3)
  mis <- c(c, p)
  
  rt <- rbind(c(cancer_type, ""), s12, s34, mis)
  
  return(rt)
}



# formula: transition probability at x year to probability at y year-----------------
# TP to rate to TP

prp <- function(tp_x, tp_x_year, tp_y_year=1){
  
  rate <- -log(1 - tp_x)/tp_x_year
  
  tp_y <- 1 - exp(-rate * tp_y_year)
  
  return(tp_y)
  
}


# Creat cf file -----------------------------------------------------------


cf_gen <- function(mod_exp, mod_wei, mod_gom){
  
  can <- list()
  
  for (i in c("exp", "wei", "gom")) {
    
    coef <- get(paste0("mod_", i))$coefficients
    
    if ("age_cent60" %in% names(coef)) cf_b <- coef[(which(names(coef)=="age_cent60")):length(coef)] else
      cf_b <- coef[(which(names(coef)=="town5_Q1")):length(coef)]
    
    if(i=="exp") {
      can[[i]][["cf_b"]] <- c("Intercept"=as.numeric(coef["rate"]), cf_b)
      can[[i]][["cf_t"]] <- c()
      can[[i]][["shape"]] <- NA
    }
    if(i=="wei"){
      can[[i]][["cf_b"]] <- c("Intercept"=as.numeric(coef["scale"]), cf_b)
      can[[i]][["cf_t"]] <- c()
      can[[i]][["shape"]] <- as.numeric(exp(coef["shape"]))
    }
    if(i=="gom"){
      can[[i]][["cf_b"]] <- c("Intercept"=as.numeric(coef["rate"]), cf_b)
      can[[i]][["cf_t"]] <- c()
      can[[i]][["shape"]] <- as.numeric(coef["shape"])
    }
  }
  
  return(can)
}


# Export functions --------------------------------------------------------

# Linear
# 95% CI
lm.csv <- function(fit, out_folder){ # input = output of coxph() 
  coef <- round(summary(fit)$coefficients[,1], digits = 3)
  CI <- paste(round(summary(fit)$coefficients[,1] - 1.96*summary(fit)$coefficients[,2],digits = 3), round(summary(fit)$coefficients[,1]+ 1.96*summary(fit)$coefficients[,2],digits = 3), sep = " , ")
  CI <- paste0("(", CI, ")")
  coef <- paste(coef, CI, sep = " ")
  coef <- data.frame(Var=rownames(summary(fit)$coefficients), coef=coef)
  write.csv(coef, paste0(out_folder,"/", deparse(substitute(fit)), ".csv"))
  return(coef)
}

# Cox
# 95% CI
coxph.csv <- function(fit, out_folder){ # input = output of coxph() 
  HR <- round(summary(fit)$conf.int[,1], digits = 2)
  CI <- paste(round(summary(fit)$conf.int[,3],digits = 2), round(summary(fit)$conf.int[,4],digits = 2), sep = "-")
  CI <- paste0("(", CI, ")")
  HR <- paste(HR, CI, sep = " ")
  HR <- data.frame(Var=rownames(summary(fit)$conf.int), HR=HR)
  write.csv(HR, paste0(out_folder,"/", deparse(substitute(fit)), ".csv"))
  return(HR)
}

# parametric model
flexsurvreg.csv <- function(fit, out_folder){
  coef <- fit$coefficients
  se <- sqrt(diag(fit$cov))
  up <- round(exp(coef + 1.96*se), digits = 2)
  low <- round(exp(coef - 1.96*se), digits = 2)
  CI <- paste0("(", low, " to ", up, ")")
  HR <- round(exp(coef), digits = 2)
  HR <- paste(HR, CI, sep = " ")
  names(HR) <- names(coef)
  write.csv(HR, paste0(out_folder,"/", deparse(substitute(fit)), ".csv"))
  return(HR)
}


# summarise functions -----------------------------------------------------

# summary of prevalence of ovarian cancer by age
tb_pct <- function(he_ova, cancer_type){
  tb <- table(he_ova$age_group10, he_ova[[paste0("cancer_", cancer_type, "_1yr")]], useNA = "ifany")
  tb_prop <- prop.table(tb, 1)
  rt <- cbind(tb, tb[,1]+tb[,2], round(tb_prop[,2]*100, 2))[, c(3, 2, 4)]
  colnames(rt) <- c("Total", paste0(cancer_type, " cancer in 1yr"), "Percentage %") 
  return(rt)
}


# PSA fns -----------------------------------------------------------------
# for qol
test_psa2 <- function(qol_list, cf){
  
  qol_det <- c(cf$qol$cf_b, cf$qol$cf_t)
  
  qol_psa <- c( )
  
  for (i in 1:1000) {
    qol <- c(qol_list[[i]]$cf_b, qol_list[[i]]$cf_t)
    qol_psa <- cbind(qol_psa, qol)
  }
  
  qol_prob <- rowQuantiles(qol_psa, probs = c(0.025, 0.975))
  
  to_plot <- data.frame(
    # exclude intercept, which is relatively too large to present with other coef.
    vars = names(qol_det)[-1], 
    mean = qol_det[-1],
    ci_l = qol_prob[,1][-1],
    ci_h = qol_prob[,2][-1]
  )
  
  p <- ggplot(to_plot, aes(x = vars, y = mean))+
    geom_bar(stat = "identity")+
    geom_errorbar(aes(ymin = ci_l, ymax = ci_h), width=0.2)+
    theme_minimal()
  
  return(p)
}

# test for cost
test_psa2 <- function(cost_coef_int_psa_all, cancer_type){
  
  cost_coef <- cost_coef_int_psa_all[[cancer_type]][-1001]
  
  det_p1 <- c(cost_coef_int_psa_all[[cancer_type]][[1001]]$p1$cf_b , cost_coef_int_psa_all[[cancer_type]][[1001]]$p1$cf_t)
  psa_p1 <- c( )
  for (i in 1:1000) {
    x <- c(cost_coef[[i]]$p1$cf_b , cost_coef[[i]]$p1$cf_t)
    psa_p1 <- cbind(psa_p1, x)
  }
  prob_p1 <- rowQuantiles(psa_p1, probs = c(0.025, 0.975))
  to_plot1 <- data.frame(
    # exclude intercept, which is relatively too large to present with other coef.
    vars = names(det_p1),
    mean = det_p1,
    ci_l = prob_p1[,1],
    ci_h = prob_p1[,2],
    p = "p1"
  )
  
  det_p2 <- c(cost_coef_int_psa_all[[cancer_type]][[1001]]$p2$cf_b , cost_coef_int_psa_all[[cancer_type]][[1001]]$p2$cf_t)
  psa_p2 <- c( )
  for (i in 1:1000) {
    x <- c(cost_coef[[i]]$p2$cf_b , cost_coef[[i]]$p2$cf_t)
    psa_p2 <- cbind(psa_p2, x)
  }
  prob_p2 <- rowQuantiles(psa_p2, probs = c(0.025, 0.975))
  to_plot2 <- data.frame(
    # exclude intercept, which is relatively too large to present with other coef.
    vars = names(det_p2),
    mean = det_p2,
    ci_l = prob_p2[,1],
    ci_h = prob_p2[,2],
    p = "p2"
  )
  
  to_plot <- rbind(to_plot1, to_plot2)
  
  
  p <- ggplot(to_plot, aes(x = vars, y = mean))+
    geom_bar(stat = "identity")+
    geom_errorbar(aes(ymin = ci_l, ymax = ci_h), width=0.2)+
    facet_wrap( ~ p, nrow = 2)+
    theme_minimal()
  
  return(p)
}

# test for accuracy data
test_psa3 <- function(accuracy){
  
  psa <- accuracy[-1001]
  
  det_p1 <- accuracy[[1001]]$accuracy_early[!grepl("cost", names(accuracy[[1001]]$accuracy_early))]
  
  psa_p1 <- c( )
  for (i in 1:1000) {
    x <- psa[[i]]$accuracy_early[!grepl("cost", names(psa[[i]]$accuracy_early))]
    psa_p1 <- cbind(psa_p1, x)
  }
  
  prob_p1 <- rowQuantiles(psa_p1, probs = c(0.025, 0.975))
  
  to_plot1 <- data.frame(
    # exclude intercept, which is relatively too large to present with other coef.
    vars = names(det_p1),
    mean = det_p1,
    ci_l = prob_p1[,1],
    ci_h = prob_p1[,2],
    p = "early"
  )
  
  det_p2 <- accuracy[[1001]]$accuracy_late[!grepl("cost", names(accuracy[[1001]]$accuracy_late))]
  
  psa_p2 <- c( )
  for (i in 1:1000) {
    x <- psa[[i]]$accuracy_late[!grepl("cost", names(psa[[i]]$accuracy_late))]
    psa_p2 <- cbind(psa_p2, x)
  }
  
  prob_p2 <- rowQuantiles(psa_p2, probs = c(0.025, 0.975))
  
  to_plot2 <- data.frame(
    # exclude intercept, which is relatively too large to present with other coef.
    vars = names(det_p2),
    mean = det_p2,
    ci_l = prob_p2[,1],
    ci_h = prob_p2[,2],
    p = "late"
  )
  
  to_plot <- rbind(to_plot1, to_plot2)
  
  
  p <- ggplot(to_plot, aes(x = vars, y = mean))+
    geom_bar(stat = "identity")+
    geom_errorbar(aes(ymin = ci_l, ymax = ci_h), width=0.2)+
    facet_wrap( ~ p, nrow = 2)+
    theme_minimal()
  
  return(p)
}

# test for stage_shift data
test_psa4 <- function(stage_shift){
  
  psa <- stage_shift[-1001]
  
  det_p1 <- do.call(c, stage_shift[[1001]])
  
  psa_p1 <- c( )
  for (i in 1:1000) {
    x <- do.call(c, psa[[i]])
    psa_p1 <- cbind(psa_p1, x)
  }
  
  prob_p1 <- rowQuantiles(psa_p1, probs = c(0.025, 0.975))
  
  to_plot <- data.frame(
    vars = names(det_p1),
    mean = det_p1,
    ci_l = prob_p1[,1],
    ci_h = prob_p1[,2]
  )
  
  p <- ggplot(to_plot, aes(x = vars, y = mean))+
    geom_bar(stat = "identity")+
    geom_errorbar(aes(ymin = ci_l, ymax = ci_h), width=0.2)+
    theme_minimal()
  
  return(p)
}

# test for FP data
test_psa5 <- function(FP_data){
  
  psa <- FP_data[-1001]
  
  det_p1 <- do.call(c, FP_data[[1001]])
  det_p1 <- det_p1[!grepl("cost", names(det_p1))]
  
  
  psa_p1 <- c( )
  for (i in 1:1000) {
    x <- do.call(c, psa[[i]])
    x <- x[!grepl("cost", names(x))]
    psa_p1 <- cbind(psa_p1, x)
  }
  
  prob_p1 <- rowQuantiles(psa_p1, probs = c(0.025, 0.975))
  
  to_plot <- data.frame(
    vars = names(det_p1),
    mean = det_p1,
    ci_l = prob_p1[,1],
    ci_h = prob_p1[,2]
  )
  
  p <- ggplot(to_plot, aes(x = vars, y = mean))+
    geom_bar(stat = "identity")+
    geom_errorbar(aes(ymin = ci_l, ymax = ci_h), width=0.2)+
    theme_minimal()
  
  return(p)
}



# define formula convert mean and 95% CI to beta distribution parameters
sample_beta <- function(mu, lb, ub){
  
  z <- qnorm(0.975)
  
  se <- (ub - lb)/(2*z)
  
  alpha <- ((1 - mu)/se^2 - 1/mu) * mu^2
  beta <- alpha * (1/mu - 1)
  
  rt <- rbeta(1000, alpha, beta)
  
  rt <- c(rt, unname(mu))
  
  return(rt)
}

sample_beta2 <- function(yes, total){
  
  alpha <- yes
  
  beta <- total - yes
  
  rt <- rbeta(1000, alpha, beta)
  
  rt <- c(rt, yes/total)
  
  return(rt)
}

sample_norm <- function(mu, lb, ub){
  
  z <- qnorm(0.975)
  
  se <- (ub - lb)/(2*z)
  
  rt <- rnorm(1000, mu, se)

  rt <- c(rt, unname(mu))
  
  return(rt) 
}





