#################################################################
##                  linke to cancer stage data                 ##
#################################################################
rm(list = ls())

source("//Qmcr/dfs/WIPH/Primary Care/OVATOOLS/r_project/scripts/file_path.R")

source(file.path(r_wd, "scripts/fn_data_clean.R"))

he_ova <- readRDS(file.path(work_data, "he_ova_step3_withSymptom_prescr.RDS"))

stage <- read_dta(file.path(gen_wd, "Data/Cancer registry files/Cancer_registry_clean.dta"))


# Merge -------------------------------------------------------------------

## Ovary -------------------------------------------------------------------

# from NCRAS
ova_stage <- stage %>% 
  select(epatid, cancerdate, site, stage, invasive_ovary, behaviour) %>% 
  # 18 coded for ovary
  filter(site==18) %>%   
  # only keep the earliest records
  group_by(epatid) %>% 
  mutate(cancerdate_min=min(cancerdate)) %>% 
  filter(cancerdate==cancerdate_min) %>% 
  distinct() %>% 
  select(-site, -cancerdate_min) %>% 
  rename(ova_date_ncras=cancerdate, 
         invasive_ova=invasive_ovary,
         stage_ova=stage,
         behaviour_ova=behaviour
         )

# merge into the current study population
he_ova <- he_ova %>% 
  # select(epatid, index_date, cancer_ova_icd, ova_icd_date, 
  #        cancer_ova_1yr, ova_1yr_date) %>% 
  left_join(ova_stage) %>% 
# remove IDs with ova_date_ncras < index date
# 237 removed
  filter(is.na(ova_date_ncras) | 
           !is.na(ova_date_ncras) & ova_date_ncras >= index_date) 

# # check
# table(he_ova$cancer_ova_1yr, he_ova$stage_ova, useNA = "ifany")
# table(he_ova$cancer_ova_icd, he_ova$stage_ova, useNA = "ifany")

## Lung --------------------------------------------------------------------
# from NCRAS
lung_stage <- stage %>% 
  select(epatid, cancerdate, site, stage, behaviour) %>% 
  # 11 coded for lung
  filter(site==11) %>%   
  # only keep the earliest records
  group_by(epatid) %>% 
  mutate(cancerdate_min=min(cancerdate)) %>% 
  filter(cancerdate==cancerdate_min) %>% 
  distinct() %>% 
  select(-site, -cancerdate_min) %>% 
  rename(lung_date_ncras=cancerdate, 
         stage_lung=stage,
         behaviour_lung=behaviour
  )

# merge into the current study population
he_ova <- he_ova %>% 
  # select(epatid, index_date, cancer_ova_icd, ova_icd_date, 
  #        cancer_ova_1yr, ova_1yr_date) %>% 
  left_join(lung_stage)

## Pancreatic --------------------------------------------------------------
# from NCRAS
panc_stage <- stage %>% 
  select(epatid, cancerdate, site, stage, behaviour) %>% 
  # 19 coded for panc 
  filter(site==19) %>%   
  # only keep the earliest records
  group_by(epatid) %>% 
  mutate(cancerdate_min=min(cancerdate)) %>% 
  filter(cancerdate==cancerdate_min) %>% 
  distinct() %>% 
  select(-site, -cancerdate_min) %>% 
  rename(panc_date_ncras=cancerdate, 
         stage_panc=stage,
         behaviour_panc=behaviour
  )

# merge into the current study population
he_ova <- he_ova %>% 
  # select(epatid, index_date, cancer_ova_icd, ova_icd_date, 
  #        cancer_ova_1yr, ova_1yr_date) %>% 
  left_join(panc_stage)


## Uterus ------------------------------------------------------------------
uter_stage <- stage %>% 
  select(epatid, cancerdate, site, stage, behaviour) %>% 
  # 23 coded for uter 
  filter(site==23) %>%   
  # only keep the earliest records
  group_by(epatid) %>% 
  mutate(cancerdate_min=min(cancerdate)) %>% 
  filter(cancerdate==cancerdate_min) %>% 
  distinct() %>% 
  select(-site, -cancerdate_min) %>% 
  rename(uter_date_ncras=cancerdate, 
         stage_uter=stage,
         behaviour_uter=behaviour
  )

# merge into the current study population
he_ova <- he_ova %>% 
  # select(epatid, index_date, cancer_ova_icd, ova_icd_date, 
  #        cancer_ova_1yr, ova_1yr_date) %>% 
  left_join(uter_stage)


## Lower GI ----------------------------------------------------------------
loGI_stage <- stage %>% 
  select(epatid, cancerdate, site, stage, behaviour) %>% 
  # 6 coded for loGI 
  filter(site==6) %>%   
  # only keep the earliest records
  group_by(epatid) %>% 
  mutate(cancerdate_min=min(cancerdate)) %>% 
  filter(cancerdate==cancerdate_min) %>% 
  distinct() %>% 
  select(-site, -cancerdate_min) %>% 
  rename(loGI_date_ncras=cancerdate, 
         stage_loGI=stage,
         behaviour_loGI=behaviour
  )

# merge into the current study population
he_ova <- he_ova %>% 
  # select(epatid, index_date, cancer_ova_icd, ova_icd_date, 
  #        cancer_ova_1yr, ova_1yr_date) %>% 
  left_join(loGI_stage)



# update cancer variables with new input from NCRAS -----------------------

## Ovary-----------------

# no ovarian cancer cases earlier than index date
# do not count in borderline i.e. invasive==0
he_ova <- he_ova %>% 
  mutate(
    # all incident
    cancer_ova_icd = case_when(
      !is.na(invasive_ova) & invasive_ova==0 ~ 0,
      cancer_ova_icd==0 & !is.na(ova_date_ncras) ~ 1,
      .default = cancer_ova_icd),
    ova_icd_date = case_when(
      cancer_ova_icd==0 ~ NA,
      is.na(ova_icd_date) & !is.na(ova_date_ncras) ~ ova_date_ncras,
      !is.na(ova_icd_date) & !is.na(ova_date_ncras) & 
        ova_icd_date > ova_date_ncras ~ ova_date_ncras,
      .default = ova_icd_date
    ),
    # 1 year
    ova_1yr_date = case_when(
      !is.na(ova_icd_date) & 
        as.numeric(ova_icd_date - index_date) <=365.25 ~ ova_icd_date,
      .default = NA
    ),
    cancer_ova_1yr = case_when(
      !is.na(ova_1yr_date) ~ 1,
      .default = 0
    )
  )

## Lung------------------

he_ova <- he_ova %>% 
  mutate(
    # baseline first
    b_lung = case_when(
      b_lung==0 & !is.na(lung_date_ncras) & lung_date_ncras < index_date ~ 1, 
      .default = b_lung
    ),
    # all incident
    cancer_lung_icd = case_when(
      b_lung==1 ~ 0,
      cancer_lung_icd==0 & !is.na(lung_date_ncras) & lung_date_ncras>=index_date ~ 1, 
      .default = cancer_lung_icd
    ),
    lung_icd_date = case_when(
      cancer_lung_icd==0 ~ NA,
      is.na(lung_icd_date) & !is.na(lung_date_ncras) & 
        lung_date_ncras >= index_date ~ lung_date_ncras,
      !is.na(lung_icd_date) & !is.na(lung_date_ncras) & 
        lung_date_ncras >= index_date &
        lung_icd_date > lung_date_ncras ~ lung_date_ncras,
      .default = lung_icd_date
    ),
    # 1 year
    lung_1yr_date = case_when(
      b_lung==1 ~ NA,
      !is.na(lung_icd_date) & 
        as.numeric(lung_icd_date - index_date) <=365.25 ~ lung_icd_date,
      .default = NA
    ),
    cancer_lung_1yr = case_when(
      !is.na(lung_1yr_date) ~ 1,
      .default = 0
    )
  )

## Panc ----------------

he_ova <- he_ova %>% 
  mutate(
    # baseline first
    b_panc = case_when(
      b_panc==0 & !is.na(panc_date_ncras) & panc_date_ncras < index_date ~ 1, 
      .default = b_panc
    ),
    # all incident
    cancer_panc_icd = case_when(
      b_panc==1 ~ 0,
      cancer_panc_icd==0 & !is.na(panc_date_ncras) & panc_date_ncras>=index_date ~ 1, 
      .default = cancer_panc_icd
    ),
    panc_icd_date = case_when(
      cancer_panc_icd==0 ~ NA,
      is.na(panc_icd_date) & !is.na(panc_date_ncras) & 
        panc_date_ncras >= index_date ~ panc_date_ncras,
      !is.na(panc_icd_date) & !is.na(panc_date_ncras) & 
        panc_date_ncras >= index_date &
        panc_icd_date > panc_date_ncras ~ panc_date_ncras,
      .default = panc_icd_date
    ),
    # 1 year
    panc_1yr_date = case_when(
      b_panc==1 ~ NA,
      !is.na(panc_icd_date) & 
        as.numeric(panc_icd_date - index_date) <=365.25 ~ panc_icd_date,
      .default = NA
    ),
    cancer_panc_1yr = case_when(
      !is.na(panc_1yr_date) ~ 1,
      .default = 0
    )
  )

## uter ---------------------
he_ova <- he_ova %>% 
  mutate(
    # baseline first
    b_uter = case_when(
      b_uter==0 & !is.na(uter_date_ncras) & uter_date_ncras < index_date ~ 1, 
      .default = b_uter
    ),
    # all incident
    cancer_uter_icd = case_when(
      b_uter==1 ~ 0,
      cancer_uter_icd==0 & !is.na(uter_date_ncras) & uter_date_ncras>=index_date ~ 1, 
      .default = cancer_uter_icd
    ),
    uter_icd_date = case_when(
      cancer_uter_icd==0 ~ NA,
      is.na(uter_icd_date) & !is.na(uter_date_ncras) & 
        uter_date_ncras >= index_date ~ uter_date_ncras,
      !is.na(uter_icd_date) & !is.na(uter_date_ncras) & 
        uter_date_ncras >= index_date &
        uter_icd_date > uter_date_ncras ~ uter_date_ncras,
      .default = uter_icd_date
    ),
    # 1 year
    uter_1yr_date = case_when(
      b_uter==1 ~ NA,
      !is.na(uter_icd_date) & 
        as.numeric(uter_icd_date - index_date) <=365.25 ~ uter_icd_date,
      .default = NA
    ),
    cancer_uter_1yr = case_when(
      !is.na(uter_1yr_date) ~ 1,
      .default = 0
    )
  )

## loGI ---------------------
he_ova <- he_ova %>% 
  mutate(
    # baseline first
    b_loGI = case_when(
      b_loGI==0 & !is.na(loGI_date_ncras) & loGI_date_ncras < index_date ~ 1, 
      .default = b_loGI
    ),
    # all incident
    cancer_loGI_icd = case_when(
      b_loGI==1 ~ 0,
      cancer_loGI_icd==0 & !is.na(loGI_date_ncras) & loGI_date_ncras>=index_date ~ 1, 
      .default = cancer_loGI_icd
    ),
    loGI_icd_date = case_when(
      cancer_loGI_icd==0 ~ NA,
      is.na(loGI_icd_date) & !is.na(loGI_date_ncras) & 
        loGI_date_ncras >= index_date ~ loGI_date_ncras,
      !is.na(loGI_icd_date) & !is.na(loGI_date_ncras) & 
        loGI_date_ncras >= index_date &
        loGI_icd_date > loGI_date_ncras ~ loGI_date_ncras,
      .default = loGI_icd_date
    ),
    # 1 year
    loGI_1yr_date = case_when(
      b_loGI==1 ~ NA,
      !is.na(loGI_icd_date) & 
        as.numeric(loGI_icd_date - index_date) <=365.25 ~ loGI_icd_date,
      .default = NA
    ),
    cancer_loGI_1yr = case_when(
      !is.na(loGI_1yr_date) ~ 1,
      .default = 0
    )
  )


# Update other cancers ----------------------------------------------------

all_stage <- stage %>% 
  select(epatid, cancerdate, site, stage, invasive_ovary, behaviour) %>% 
  # invasive==0 only for ovarian cancer 
  filter(invasive_ovary !=0 | is.na(invasive_ovary)) %>%   
  # only keep the earliest records
  group_by(epatid) %>% 
  mutate(cancerdate_min=min(cancerdate)) %>% 
  filter(cancerdate==cancerdate_min) %>% 
  distinct(epatid, .keep_all = T) %>% 
  select(-cancerdate_min) %>% 
  rename(all_date_ncras=cancerdate, 
         stage_all=stage,
         behaviour_all=behaviour
  )

# merge into the current study population
he_ova <- he_ova %>% 
  left_join(all_stage) %>% 
  # remove possible baseline ovarian cancer
  # actually, none, as they have been removed from he_ova when merging ova_stage
  filter(!(site==18 & !is.na(all_date_ncras) & all_date_ncras<index_date)) %>% 
  select(-site) %>%
  # baseline cancer
  mutate(
    b_cancer_nonOva = case_when(
      b_cancer_nonOva==0 & !is.na(all_date_ncras) & all_date_ncras < index_date ~ 1,
      .default = b_cancer_nonOva
    ),
    b_cancer_others = case_when(
      b_cancer_nonOva==1 & b_lung==0 & b_panc==0 & b_uter==0 & b_loGI==0 ~ 1,
      .default = 0
    )
  ) %>% 
  # incident 
  mutate(
    cancer_icd = case_when(
      b_cancer_nonOva==1 ~ 0,
      cancer_icd==0 & !is.na(all_date_ncras) & all_date_ncras>=index_date ~ 1, 
      .default = cancer_icd
    ),
    cancer_icd_date = case_when(
      cancer_icd==0 ~ NA,
      is.na(cancer_icd_date) & !is.na(all_date_ncras) & 
        all_date_ncras >= index_date ~ all_date_ncras,
      !is.na(cancer_icd_date) & !is.na(all_date_ncras) & 
        all_date_ncras >= index_date &
        cancer_icd_date > all_date_ncras ~ all_date_ncras,
      .default = cancer_icd_date
    ),
    # 1 year
    cancer_1yr_date = case_when(
      b_cancer_nonOva==1 ~ NA,
      !is.na(cancer_icd_date) & 
        as.numeric(cancer_icd_date - index_date) <=365.25 ~ cancer_icd_date,
      .default = NA
    ),
    cancer_1yr = case_when(
      !is.na(cancer_1yr_date) ~ 1,
      .default = 0
    )
  ) %>% 
  mutate(
  # other incident cancers, not ova, lung, panc, loGI, or uter
    cancer_other_1yr = if_else(
      cancer_1yr==1 & 
        cancer_uter_1yr + cancer_ova_1yr + cancer_lung_1yr + cancer_panc_1yr +
        cancer_loGI_1yr==0, 1, 0),
    cancer_other_1yr_date = if_else(cancer_other_1yr==1, cancer_1yr_date, NA)
  ) %>% 
  mutate(
    cancer_other_icd = if_else(
      cancer_icd==1 & 
        cancer_uter_icd + cancer_ova_icd + cancer_lung_icd + cancer_panc_icd +
        cancer_loGI_icd==0, 1, 0),
    # note, as defined above, cancer_other_icd do not include all cancer_other_1yr, 
    # because some cancer_other_1yr developed non-other cancer e.g. ova, lung, etc later 
    # and can not be included ad cancer_other_icd
    # But we need to include all cancer_other_1yr, so make the patch here
    cancer_other_icd = if_else(cancer_other_1yr==1, 1, cancer_other_icd),
    cancer_other_icd_date = if_else(cancer_other_icd==1, cancer_icd_date, NA)
  )

# change the other cancer date names for easy use later
he_ova <- he_ova %>% 
  rename(
    other_icd_date = cancer_other_icd_date,
    other_1yr_date = cancer_other_1yr_date
  )

# table(he_ova$cancer_uter_icd, is.na(he_ova$uter_icd_date), useNA = "ifany")
# table(he_ova$cancer_uter_1yr, is.na(he_ova$uter_1yr_date), useNA = "ifany")
# table(he_ova$cancer_uter_1yr, he_ova$cancer_uter_icd, useNA = "ifany")
# summary(he_ova$uter_1yr_date)
# 
# 
# table(he_ova$cancer_icd, is.na(he_ova$cancer_icd_date), useNA = "ifany")
# table(he_ova$cancer_1yr, is.na(he_ova$cancer_1yr_date), useNA = "ifany")
# table(he_ova$cancer_1yr, he_ova$cancer_icd, useNA = "ifany")
# summary(he_ova$cancer_1yr_date)


# Recode some vars --------------------------------------------------------

he_ova <- he_ova %>% 
  mutate(
    stage_ova_imp = case_when(
      cancer_ova_1yr==1 & stage_ova==1 ~ "stage 1",
      cancer_ova_1yr==1 & stage_ova==2 ~ "stage 2",
      cancer_ova_1yr==1 & stage_ova==3 ~ "stage 3",
      cancer_ova_1yr==1 & stage_ova==4 ~ "stage 4",
      cancer_ova_1yr==0 ~ "no cancer",
      .default = NA),
    
    stage_loGI_imp = case_when(
      cancer_loGI_1yr==1 & stage_loGI==1 ~ "stage 1",
      cancer_loGI_1yr==1 & stage_loGI==2 ~ "stage 2",
      cancer_loGI_1yr==1 & stage_loGI==3 ~ "stage 3",
      cancer_loGI_1yr==1 & stage_loGI==4 ~ "stage 4",
      cancer_loGI_1yr==0 ~ "no cancer",
      .default = NA),
    
    stage_panc_imp = case_when(
      cancer_panc_1yr==1 & stage_panc==1 ~ "stage 1",
      cancer_panc_1yr==1 & stage_panc==2 ~ "stage 2",
      cancer_panc_1yr==1 & stage_panc==3 ~ "stage 3",
      cancer_panc_1yr==1 & stage_panc==4 ~ "stage 4",
      cancer_panc_1yr==0 ~ "no cancer",
      .default = NA),  
    
    stage_uter_imp = case_when(
      cancer_uter_1yr==1 & stage_uter==1 ~ "stage 1",
      cancer_uter_1yr==1 & stage_uter==2 ~ "stage 2",
      cancer_uter_1yr==1 & stage_uter==3 ~ "stage 3",
      cancer_uter_1yr==1 & stage_uter==4 ~ "stage 4",
      cancer_uter_1yr==0 ~ "no cancer",
      .default = NA),  
    
    stage_lung_imp = case_when(
      cancer_lung_1yr==1 & stage_lung==1 ~ "stage 1",
      cancer_lung_1yr==1 & stage_lung==2 ~ "stage 2",
      cancer_lung_1yr==1 & stage_lung==3 ~ "stage 3",
      cancer_lung_1yr==1 & stage_lung==4 ~ "stage 4",
      cancer_lung_1yr==0 ~ "no cancer",
      .default = NA),
    
    stage_other_imp = case_when(
      cancer_other_1yr==1 & stage_all==1 ~ "stage 1",
      cancer_other_1yr==1 & stage_all==2 ~ "stage 2",
      cancer_other_1yr==1 & stage_all==3 ~ "stage 3",
      cancer_other_1yr==1 & stage_all==4 ~ "stage 4",
      cancer_other_1yr==0 ~ "no cancer",
      .default = NA)
  ) %>% 
  mutate(fu_time = as.numeric(censor_date - index_date)/365.25) %>% 
  # add apporox 1 day to avoid 0 fu time
  mutate(fu_time = if_else(fu_time==0, fu_time + 1/365.25, fu_time),
         age_cent60 = (age-60)/10,
         death_cancer = if_else(is.na(death_cancer), 0, death_cancer),
         death_other = if_else(death==1 & death_cancer==0, 1, 0)) %>% 
  mutate(
    ethn4=case_when(
      ethn_comb_imp=="Mixed or others" ~ "Others",
      .default = ethn_comb_imp
    ),
    ethn4=factor(ethn4, levels=c("White", "Asian", "Black", "Others"))
  ) %>% 
  mutate(
    # other ethnicity has very small number, combine with White
    ethn3 = case_when(ethn4 =="Others" ~"White",
                      .default = ethn4),
    ethn3 = factor(ethn3, levels= c("White", "Asian", "Black"))
  ) %>% 
  mutate(
    town5 = factor(town5, levels=c("Q3", "Q1", "Q2", "Q4", "Q5"))
  )

# for cost analysis, including stage data not limited to diagnosis within 1 year
he_ova <- he_ova %>% 
  mutate(
    stage_ova_allicd = case_when(
      cancer_ova_icd==1 & stage_ova==1 ~ "stage 1",
      cancer_ova_icd==1 & stage_ova==2 ~ "stage 2",
      cancer_ova_icd==1 & stage_ova==3 ~ "stage 3",
      cancer_ova_icd==1 & stage_ova==4 ~ "stage 4",
      cancer_ova_icd==0 ~ "no cancer",
      .default = NA),
    
    stage_loGI_allicd = case_when(
      cancer_loGI_icd==1 & stage_loGI==1 ~ "stage 1",
      cancer_loGI_icd==1 & stage_loGI==2 ~ "stage 2",
      cancer_loGI_icd==1 & stage_loGI==3 ~ "stage 3",
      cancer_loGI_icd==1 & stage_loGI==4 ~ "stage 4",
      cancer_loGI_icd==0 ~ "no cancer",
      .default = NA),
    
    stage_panc_allicd = case_when(
      cancer_panc_icd==1 & stage_panc==1 ~ "stage 1",
      cancer_panc_icd==1 & stage_panc==2 ~ "stage 2",
      cancer_panc_icd==1 & stage_panc==3 ~ "stage 3",
      cancer_panc_icd==1 & stage_panc==4 ~ "stage 4",
      cancer_panc_icd==0 ~ "no cancer",
      .default = NA),  
    
    stage_uter_allicd = case_when(
      cancer_uter_icd==1 & stage_uter==1 ~ "stage 1",
      cancer_uter_icd==1 & stage_uter==2 ~ "stage 2",
      cancer_uter_icd==1 & stage_uter==3 ~ "stage 3",
      cancer_uter_icd==1 & stage_uter==4 ~ "stage 4",
      cancer_uter_icd==0 ~ "no cancer",
      .default = NA),  
    
    stage_lung_allicd = case_when(
      cancer_lung_icd==1 & stage_lung==1 ~ "stage 1",
      cancer_lung_icd==1 & stage_lung==2 ~ "stage 2",
      cancer_lung_icd==1 & stage_lung==3 ~ "stage 3",
      cancer_lung_icd==1 & stage_lung==4 ~ "stage 4",
      cancer_lung_icd==0 ~ "no cancer",
      .default = NA),
    
    stage_other_allicd = case_when(
      cancer_other_icd==1 & stage_all==1 ~ "stage 1",
      cancer_other_icd==1 & stage_all==2 ~ "stage 2",
      cancer_other_icd==1 & stage_all==3 ~ "stage 3",
      cancer_other_icd==1 & stage_all==4 ~ "stage 4",
      cancer_other_icd==0 ~ "no cancer",
      .default = NA)
  )

saveRDS(he_ova, file.path(work_data, "he_ova_step4_withSymptom_prescr_NCRASupdated.RDS"))


