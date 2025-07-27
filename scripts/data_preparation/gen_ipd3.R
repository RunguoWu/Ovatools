rm(list = ls())

source("//Qmcr/dfs/WIPH/Primary Care/OVATOOLS/r_project/scripts/file_path.R")
source(file.path(r_wd, "scripts/fn_data_clean.R"))
library(fastDummies)
library(tidyverse)
library(nnet)
source(file.path(r_wd, "scripts/fn_MI.R"))
he_ova <- readRDS(file.path(work_data, "he_ova_step4_withSymptom_prescr_NCRASupdated.RDS"))
# merge symptoms for MI
he_ova <- merge_symptoms(he_ova)

# Combine benign baseline and 1 yr as history -----------------------------
# baseline  + 1 yr
# for QoL, not mortality
dt <- he_ova %>% 
  mutate(
    ben_gynae_1 = if_else(b_benign_gynae==1 |benign_1yr==1, 1, 0)
  ) %>% 
  select(epatid, age_cent60, ben_gynae_1, age_group_40_80, age_group_50_80, age_group_60_80, 
         ethn3, ethn4, town5, b_uter, b_loGI, b_lung, b_panc, b_cancer_others, ben_gynae_1,
         cancer_ova_1yr, cancer_uter_1yr, cancer_loGI_1yr, cancer_lung_1yr, 
         cancer_panc_1yr, cancer_other_1yr,
         stage_ova_imp, stage_loGI_imp, stage_uter_imp, stage_lung_imp,
         stage_panc_imp, stage_other_imp
  )


# Predict prevalence and stage --------------------------------------------
# do not use summary by categories, because that results in some categories 
# with 0 prevalence

## Ova ----
imp_list <- do_MI(he_ova, cancer_type="ova", N=30)
pred_all <- list()
for (i in 1: length(imp_list)){
  
  x <- imp_list[[i]] %>% select(epatid, stage_late)
  
  d4s <- dt %>% left_join(x) %>% 
    mutate(type_stage = case_when(cancer_ova_1yr==1 & stage_late==0 ~ "early",
                                  cancer_ova_1yr==1 & stage_late==1 ~ "late",
                                  .default = "none"
    ))
  
  pred <- tryCatch({
    fit <- multinom(type_stage ~ age_group_40_80 + ethn4 + town5, d4s)
    predict(fit, dt, type="prob")  
  }, error = function(e){
    return(NULL)})
  
  pred_all[[i]] <- pred
}
pred_all <- pred_all[!sapply(pred_all, is.null)]

dt$ova_early <- rowMeans(sapply(pred_all, function(x){x[,"early"]}))
dt$ova_late <- rowMeans(sapply(pred_all, function(x){x[,"late"]}))

## loGI ----
imp_list <- do_MI(he_ova, cancer_type="loGI", N=30)
pred_all <- list()
for (i in 1: length(imp_list)){
  
  x <- imp_list[[i]] %>% select(epatid, stage_late)
  
  d4s <- dt %>% left_join(x) %>% 
    filter(b_loGI==0) %>% 
    mutate(type_stage = case_when(cancer_loGI_1yr==1 & stage_late==0 ~ "early",
                                  cancer_loGI_1yr==1 & stage_late==1 ~ "late",
                                  .default = "none"
    ))
  
  pred <- tryCatch({
    fit <- multinom(type_stage ~ age_group_50_80 + ethn4 + town5, d4s)
    predict(fit, dt, type="prob")  
  }, error = function(e){
    return(NULL)})
  
  pred_all[[i]] <- pred
}
pred_all <- pred_all[!sapply(pred_all, is.null)]

dt$loGI_early <- rowMeans(sapply(pred_all, function(x){x[,"early"]}))
dt$loGI_late <- rowMeans(sapply(pred_all, function(x){x[,"late"]}))

## uter ----
imp_list <- do_MI(he_ova, cancer_type="uter", N=30)
pred_all <- list()
for (i in 1: length(imp_list)){
  
  x <- imp_list[[i]] %>% select(epatid, stage_late)
  
  d4s <- dt %>% left_join(x) %>% 
    filter(b_uter==0) %>% 
    mutate(type_stage = case_when(cancer_uter_1yr==1 & stage_late==0 ~ "early",
                                  cancer_uter_1yr==1 & stage_late==1 ~ "late",
                                  .default = "none"
    ))
  
  pred <- tryCatch({
    fit <- multinom(type_stage ~ age_group_50_80 + ethn4 + town5, d4s)
    predict(fit, dt, type="prob")  
  }, error = function(e){
    return(NULL)})
  
  pred_all[[i]] <- pred
}
pred_all <- pred_all[!sapply(pred_all, is.null)]

dt$uter_early <- rowMeans(sapply(pred_all, function(x){x[,"early"]}))
dt$uter_late <- rowMeans(sapply(pred_all, function(x){x[,"late"]}))

## lung ---- # use ethn3 due to low obser for others
imp_list <- do_MI(he_ova, cancer_type="lung", N=30)
pred_all <- list()
for (i in 1: length(imp_list)){
  
  x <- imp_list[[i]] %>% select(epatid, stage_late)
  
  d4s <- dt %>% left_join(x) %>% 
    filter(b_lung==0) %>% 
    mutate(type_stage = case_when(cancer_lung_1yr==1 & stage_late==0 ~ "early",
                                  cancer_lung_1yr==1 & stage_late==1 ~ "late",
                                  .default = "none"
    ))
  
  pred <- tryCatch({
    fit <- multinom(type_stage ~ age_group_60_80 + ethn3 + town5, d4s)
    predict(fit, dt, type="prob")  
  }, error = function(e){
    return(NULL)})
  
  pred_all[[i]] <- pred
}
pred_all <- pred_all[!sapply(pred_all, is.null)]

dt$lung_early <- rowMeans(sapply(pred_all, function(x){x[,"early"]}))
dt$lung_late <- rowMeans(sapply(pred_all, function(x){x[,"late"]}))

## panc ---- # use ethn3 due to low obser for others
imp_list <- do_MI(he_ova, cancer_type="panc", N=30)
pred_all <- list()
for (i in 1: length(imp_list)){
  
  x <- imp_list[[i]] %>% select(epatid, stage_late)
  
  d4s <- dt %>% left_join(x) %>% 
    filter(b_panc==0) %>% 
    mutate(type_stage = case_when(cancer_panc_1yr==1 & stage_late==0 ~ "early",
                                  cancer_panc_1yr==1 & stage_late==1 ~ "late",
                                  .default = "none"
    ))
  
  pred <- tryCatch({
    fit <- multinom(type_stage ~ age_group_60_80 + ethn3 + town5, d4s)
    predict(fit, dt, type="prob")  
  }, error = function(e){
    return(NULL)})
  
  pred_all[[i]] <- pred
}
pred_all <- pred_all[!sapply(pred_all, is.null)]

dt$panc_early <- rowMeans(sapply(pred_all, function(x){x[,"early"]}))
dt$panc_late <- rowMeans(sapply(pred_all, function(x){x[,"late"]}))

## other ----
imp_list <- do_MI(he_ova, cancer_type="other", N=30)
pred_all <- list()
for (i in 1: length(imp_list)){
  
  x <- imp_list[[i]] %>% select(epatid, stage_late)
  
  d4s <- dt %>% left_join(x) %>% 
    filter(b_cancer_others==0) %>% 
    mutate(type_stage = case_when(cancer_other_1yr==1 & stage_late==0 ~ "early",
                                  cancer_other_1yr==1 & stage_late==1 ~ "late",
                                  .default = "none"
    ))
  
  pred <- tryCatch({
    fit <- multinom(type_stage ~ age_group_40_80 + ethn4 + town5, d4s)
    predict(fit, dt, type="prob")  
  }, error = function(e){
    return(NULL)})
  
  pred_all[[i]] <- pred
}
pred_all <- pred_all[!sapply(pred_all, is.null)]

dt$other_early <- rowMeans(sapply(pred_all, function(x){x[,"early"]}))
dt$other_late <- rowMeans(sapply(pred_all, function(x){x[,"late"]}))

# Create IPD for export ---------------------------------------------------
ipd <- dt %>% 
  mutate(
    prevalence=ova_early+ova_late,
    pred_uter = if_else(b_uter==1, 0, uter_early + uter_late),
    pred_loGI = if_else(b_loGI==1, 0, loGI_early + loGI_late),
    pred_lung = if_else(b_lung==1, 0, lung_early + lung_late),
    pred_panc = if_else(b_panc==1, 0, panc_early + panc_late),
    pred_other = if_else(b_cancer_others==1, 0, other_early + other_late),
    
    ova_stage34 = ova_late/(prevalence),
    uter_stage34 = if_else(b_uter==1, 0, uter_late/pred_uter),
    loGI_stage34 = if_else(b_loGI==1, 0, loGI_late/pred_loGI),
    lung_stage34 = if_else(b_lung==1, 0, lung_late/pred_lung),
    panc_stage34 = if_else(b_panc==1, 0, panc_late/pred_panc),
    other_stage34 = if_else(b_cancer_others==1, 0, other_late/pred_other),
    
    prevalence_otc1=pred_uter+pred_loGI,
    prevalence_otc2=pred_lung+pred_panc
  ) %>% 
  select(age_cent60, ethn3, ethn4, town5, 
         b_uter, b_loGI, b_lung, b_panc, b_cancer_others, ben_gynae_1,
         prevalence, prevalence_otc1, prevalence_otc2, ova_stage34, 
         pred_uter, pred_lung, pred_panc, pred_loGI, pred_other, 
         uter_stage34, loGI_stage34, lung_stage34, panc_stage34, other_stage34
  ) %>%
  rename(cancers_other_1 = b_cancer_others, # for QoL calculation
         gynae_2 = b_uter,
         upGI_2 = b_panc,
         lung_2 = b_lung,
         loGI_2 = b_loGI
  ) %>% 
  mutate(id=1:nrow(he_ova)) %>% 
  relocate(id, .before = age_cent60)

ipd <- dummy_cols(ipd, select_columns = c("town5","ethn3", "ethn4"),
                  remove_selected_columns = TRUE, remove_first_dummy = TRUE)

saveRDS(ipd, file.path(work_data, "ipd_new20240906.RDS")) # MI used  
  
# ID dictionary -----------------------------------------------------------
id_dic <- he_ova %>% 
  mutate(id=1:nrow(he_ova)) %>% 
  select(id, epatid)

saveRDS(id_dic, file.path(work_data, "id_dic.RDS"))


# gen unique values to export ---------------------------------------------
ipd2 <- ipd %>% 
  mutate(group = group_indices(., across(-id))) %>% 
  distinct(group, .keep_all = T) %>% 
  select(-id)

saveRDS(ipd2, file.path(work_data, "ipd2_20240906.RDS"))

# generate a dictionary for mapping
id_dic2 <- ipd %>% 
  mutate(group = group_indices(., across(-id))) %>% 
  select(id, group)

saveRDS(id_dic2, file.path(work_data, "id_dic2.RDS"))

# generate a dictionary for mapping group id to epatid
id_dic3 <- merge(id_dic, id_dic2)
saveRDS(id_dic3, file.path(work_data, "id_dic3.RDS"))




