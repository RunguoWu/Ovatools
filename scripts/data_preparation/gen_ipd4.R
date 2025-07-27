#################################################################
##                   Generate IPD input file                   ##
##             Do not include prevalence and stage             ##
##                          Suit PSA                           ##
#################################################################

rm(list = ls())

source("//Qmcr/dfs/WIPH/Primary Care/OVATOOLS/r_project/scripts/file_path.R")
source(file.path(r_wd, "scripts/fn_data_clean.R"))
library(fastDummies)
library(tidyverse)
library(nnet)
source(file.path(r_wd, "scripts/fn_MI.R"))
he_ova <- readRDS(file.path(work_data, "he_ova_step4_withSymptom_prescr_NCRASupdated.RDS"))

# Create IPD for export ---------------------------------------------------
ipd <- he_ova %>% 
  mutate(
    ben_gynae_1 = if_else(b_benign_gynae==1 |benign_1yr==1, 1, 0)
  ) %>% 
  select(age_cent60, ethn3, ethn4, town5, 
         b_uter, b_loGI, b_lung, b_panc, b_cancer_others, ben_gynae_1
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

saveRDS(ipd, file.path(work_data, "ipd_new20241126.RDS")) 

# ID dictionary -----------------------------------------------------------
id_dic <- he_ova %>% 
  mutate(id=1:nrow(he_ova)) %>% 
  select(id, epatid)

saveRDS(id_dic, file.path(work_data, "id_dic.RDS"))# same as old one

# gen unique values to export ---------------------------------------------
ipd2 <- ipd %>% 
  mutate(group = group_indices(., across(-id))) %>% 
  distinct(group, .keep_all = T) %>% 
  select(-id)

saveRDS(ipd2, file.path(work_data, "ipd2_20241126.RDS")) 
# same size and order as old one, but group ids changed

# generate a dictionary for mapping
id_dic2 <- ipd %>% 
  mutate(group = group_indices(., across(-id))) %>% 
  select(id, group)

saveRDS(id_dic2, file.path(work_data, "id_dic2_20241126.RDS"))

# generate a dictionary for mapping group id to epatid
id_dic3 <- merge(id_dic, id_dic2)
saveRDS(id_dic3, file.path(work_data, "id_dic3_20241126.RDS"))
