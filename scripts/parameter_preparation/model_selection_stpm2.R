##################################################################
##           Mortality risk equations by cancer types           ##
##################################################################

# v4 use flexible parametric survival model, according to ROyston and Parmar Statist. Med. 2002
# compared to v3, v4 does not split data by stage
# QoL, mortality statistics, decision tree parameters are the same as v2

rm(list = ls())

source("//Qmcr/dfs/WIPH/Primary Care/OVATOOLS/r_project/scripts/file_path.R")
source(file.path(r_wd, "scripts/fn_data_clean.R"))
source(file.path(r_wd, "scripts/fn_MI.R"))

library(flexsurv)
library(survival)
library(fastDummies)
library(lspline)
library(rstpm2)

he_ova <- readRDS(file.path(work_data, "he_ova_step4_withSymptom_prescr_NCRASupdated.RDS"))
mod_cca <- list()
# Ova ---------------------------------------------------------------------
d4s <- he_ova %>% 
  filter(cancer_ova_1yr==1 & !is.na(stage_ova_imp)
  ) %>%
  mutate(stage_late = 
           case_when( 
             stage_ova_imp %in% c("stage 1", "stage 2") ~ 0, 
             stage_ova_imp %in% c("stage 3", "stage 4") ~ 1
           )
  ) %>% 
  mutate(fu_diag2cens = as.numeric(censor_date - ova_1yr_date)/365.25,
         fu_diag2cens = if_else(fu_diag2cens==0, fu_diag2cens + 1/365.25, fu_diag2cens)
  ) %>% 
  select(fu_diag2cens, death_cancer, death_other, age_cent60, town5, ethn3, 
         stage_late, CA125) 

fm <- Surv(fu_diag2cens, death_cancer) ~ age_cent60 + town5 + ethn3+ 
  stage_late

# fit <- stpm2(fm, d4s, df=1)
# AIC(fit)
# fit <- stpm2(fm, d4s, df=2)
# AIC(fit)
# fit <- stpm2(fm, d4s, df=3)
# AIC(fit)
# fit <- stpm2(fm, d4s, df=4)
# AIC(fit)
# fit <- stpm2(fm, d4s, df=5)
# AIC(fit)
# df=5

# test if the model perform in pop with ca125 record as well as the whole
d4s_ca125 <- d4s %>% filter(!is.na(CA125))

# test if the model perform in pop without ca125 record 
d4s_uss <- d4s %>% filter(is.na(CA125))

fit <- stpm2(fm, d4s, df=5)
plot_stpm2(fit, by_stage = 1, max_time = 8, d4s)
plot_stpm2(fit, by_stage = 1, max_time = 8, d4s_ca125) # very good performance
plot_stpm2(fit, by_stage = 1, max_time = 8, d4s_uss)

mod_cca[["ova"]] <- fit

# also check in the whole pop, with/without ca125
x <- he_ova %>% filter(!is.na(CA125))
y <- he_ova %>% filter(is.na(CA125))


# Lung --------------------------------------------------------------------
d4s <- he_ova %>% 
  filter(cancer_lung_1yr==1 & !is.na(stage_lung_imp)) %>%
  mutate(stage_late = 
           case_when( 
             stage_lung_imp %in% c("stage 1", "stage 2") ~ 0, 
             stage_lung_imp %in% c("stage 3", "stage 4") ~ 1
           )
  ) %>% 
  mutate(fu_diag2cens = as.numeric(censor_date - lung_1yr_date)/365.25,
         fu_diag2cens = if_else(fu_diag2cens==0, fu_diag2cens + 1/365.25, fu_diag2cens)
  ) %>% 
  select(fu_diag2cens, death_cancer, death_other, age_cent60, town5, ethn3, stage_late) 

fm <- Surv(fu_diag2cens, death_cancer) ~ age_cent60 + town5 + ethn3+stage_late

fit <- stpm2(fm, d4s, df=1)
AIC(fit)
fit <- stpm2(fm, d4s, df=2)
AIC(fit)
fit <- stpm2(fm, d4s, df=3)
AIC(fit)
fit <- stpm2(fm, d4s, df=4)
AIC(fit)
fit <- stpm2(fm, d4s, df=5)
AIC(fit)

# df=4
fit <- stpm2(fm, d4s, df=4)
plot_stpm2(fit, by_stage = 1, max_time = 8, d4s)
mod_cca[["lung"]] <- fit

# loGI --------------------------------------------------------------------
d4s <- he_ova %>% 
  filter(cancer_loGI_1yr==1 & !is.na(stage_loGI_imp)) %>%
  mutate(stage_late = 
           case_when( 
             stage_loGI_imp %in% c("stage 1", "stage 2") ~ 0, 
             stage_loGI_imp %in% c("stage 3", "stage 4") ~ 1
           )
  ) %>% 
  mutate(fu_diag2cens = as.numeric(censor_date - loGI_1yr_date)/365.25,
         fu_diag2cens = if_else(fu_diag2cens==0, fu_diag2cens + 1/365.25, fu_diag2cens)
  ) %>% 
  select(fu_diag2cens, death_cancer, death_other, age_cent60, town5, ethn3, 
         stage_late) 

fm <- Surv(fu_diag2cens, death_cancer) ~ age_cent60 + town5 + ethn3+stage_late

# fit <- stpm2(fm, d4s, df=1)
# AIC(fit)
# fit <- stpm2(fm, d4s, df=2)
# AIC(fit)
# fit <- stpm2(fm, d4s, df=3)
# AIC(fit)
# fit <- stpm2(fm, d4s, df=4)
# AIC(fit)
# fit <- stpm2(fm, d4s, df=5)
# AIC(fit)

# df=4
fit <- stpm2(fm, d4s, df=4)
plot_stpm2(fit, by_stage = 0, max_time = 8, d4s)
mod_cca[["loGI"]] <- fit

# uter ---------------------------------------------------------------------
d4s <- he_ova %>% 
  filter(cancer_uter_1yr==1 & !is.na(stage_uter_imp)) %>%
  mutate(stage_late = 
           case_when( 
             stage_uter_imp %in% c("stage 1", "stage 2") ~ 0, 
             stage_uter_imp %in% c("stage 3", "stage 4") ~ 1
           )
  ) %>% 
  mutate(fu_diag2cens = as.numeric(censor_date - uter_1yr_date)/365.25,
         fu_diag2cens = if_else(fu_diag2cens==0, fu_diag2cens + 1/365.25, fu_diag2cens)
  ) %>% 
  select(fu_diag2cens, death_cancer, death_other, age_cent60, town5, ethn3, 
         stage_late) 

fm <- Surv(fu_diag2cens, death_cancer) ~ age_cent60 + town5 + ethn3+stage_late

# fit <- stpm2(fm, d4s, df=1)
# AIC(fit)
# fit <- stpm2(fm, d4s, df=2)
# AIC(fit)
# fit <- stpm2(fm, d4s, df=3)
# AIC(fit)
# fit <- stpm2(fm, d4s, df=4)
# AIC(fit)
# fit <- stpm2(fm, d4s, df=5)
# AIC(fit)

#df=2
fit <- stpm2(fm, d4s, df=2)
plot_stpm2(fit, by_stage = 1, max_time = 8, d4s)
mod_cca[["uter"]] <- fit

# panc ---------------------------------------------------------------------
d4s <- he_ova %>% 
  filter(cancer_panc_1yr==1 & !is.na(stage_panc_imp)) %>%
  mutate(stage_late = 
           case_when( 
             stage_panc_imp %in% c("stage 1", "stage 2") ~ 0, 
             stage_panc_imp %in% c("stage 3", "stage 4") ~ 1
           )
  ) %>% 
  mutate(fu_diag2cens = as.numeric(censor_date - panc_1yr_date)/365.25,
         fu_diag2cens = if_else(fu_diag2cens==0, fu_diag2cens + 1/365.25, fu_diag2cens)
  ) %>% 
  select(fu_diag2cens, death_cancer, death_other, age_cent60, town5, ethn3, 
         stage_late) 

fm <- Surv(fu_diag2cens, death_cancer) ~ age_cent60 + town5 + ethn3+stage_late

# fit <- stpm2(fm, d4s, df=1)
# AIC(fit)
# fit <- stpm2(fm, d4s, df=2)
# AIC(fit)
# fit <- stpm2(fm, d4s, df=3)
# AIC(fit)
# fit <- stpm2(fm, d4s, df=4)
# AIC(fit)
# fit <- stpm2(fm, d4s, df=5)
# AIC(fit)

#df=2
fit <- stpm2(fm, d4s, df=2)
plot_stpm2(fit, by_stage = 1, max_time = 8, d4s)
mod_cca[["panc"]] <- fit

# other ---------------------------------------------------------------------
d4s <- he_ova %>% 
  filter(cancer_other_1yr==1 & !is.na(stage_other_imp)) %>%
  mutate(stage_late = 
           case_when( 
             stage_other_imp %in% c("stage 1", "stage 2") ~ 0, 
             stage_other_imp %in% c("stage 3", "stage 4") ~ 1
           )
  ) %>% 
  mutate(fu_diag2cens = as.numeric(censor_date - other_1yr_date)/365.25,
         fu_diag2cens = if_else(fu_diag2cens==0, fu_diag2cens + 1/365.25, fu_diag2cens)
  ) %>% 
  select(fu_diag2cens, death_cancer, death_other, age_cent60, town5, ethn3, 
         stage_late) 

fm <- Surv(fu_diag2cens, death_cancer) ~ age_cent60 + town5 + ethn3+stage_late
# 
# fit <- stpm2(fm, d4s, df=1)
# AIC(fit)
# fit <- stpm2(fm, d4s, df=2)
# AIC(fit)
# fit <- stpm2(fm, d4s, df=3)
# AIC(fit)
# fit <- stpm2(fm, d4s, df=4)
# AIC(fit)
# fit <- stpm2(fm, d4s, df=5)
# AIC(fit)

# df=4
fit <- stpm2(fm, d4s, df=4)
plot_stpm2(fit, by_stage = 0, max_time = 8, d4s)
mod_cca[["other"]] <- fit

# save --------------------------------------------------------------------
saveRDS(mod_cca, file.path(output, "mod_cca_stpm2.rds"))




