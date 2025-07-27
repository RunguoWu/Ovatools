#################################################################
##                  Generation of parameters                   ##
##         Use stpm2 to fit flexible parametric models         ##
#################################################################
library(bannerCommenter)

# Cancer death parameters -------------------------------------------------

rm(list = ls())

library(mice)
library(miceadds)
library(tidyverse)
library(parallel)
library(matrixStats)
library(rstpm2)
library(survival)
library(fastDummies)
library(lspline)
library(dtplyr)
library(rlang)
library(haven)
library(foreign)
library(foreach)
library(doSNOW)
library(snowfall)
library(lspline)
library(parallel)

source("//Qmcr/dfs/WIPH/Primary Care/OVATOOLS/r_project/scripts/file_path.R")
source(file.path(r_wd, "scripts/fn_MI.R"))


he_ova <- readRDS(file.path(work_data, "he_ova_step4_withSymptom_prescr_NCRASupdated.RDS"))
# merge symptoms
he_ova <- merge_symptoms(he_ova)

mod_list <- list()
mod_list2 <- list()

## Ova ----
# get imp_list
imp_list <- do_MI(he_ova, cancer_type="ova", N=50)

fm <- Surv(fu_diag2cens, death_cancer) ~ age_cent60 + 
  town5_Q1 + town5_Q2 + town5_Q4 + town5_Q5 +
  ethn3_Asian + ethn3_Black + 
  stage_late

mod_list[["ova"]] <- fit_MI_mod_stpm2(imp_list, n=30, fm, df=5)

pool(mod_list[["ova"]])

# model on complete data with average stage
dt <- gen_complete(imp_list)
fit <- stpm2(fm, data = dt, df=5)
mod_list2[["ova"]] <- fit

## Lung----
imp_list <- do_MI(he_ova, cancer_type="lung", N=50)

fm <- Surv(fu_diag2cens, death_cancer) ~ age_cent60 + 
  town5_Q1 + town5_Q2 + town5_Q4 + town5_Q5 +
  ethn3_Asian + ethn3_Black + 
  stage_late

mod_list[["lung"]] <- fit_MI_mod_stpm2(imp_list, n=30, fm, df=4)

pool(mod_list[["lung"]])

# model on complete data with average stage
dt <- gen_complete(imp_list)
fit <- stpm2(fm, data = dt, df=4)
mod_list2[["lung"]] <- fit


## loGI----
imp_list <- do_MI(he_ova, cancer_type="loGI", N=50)

fm <- Surv(fu_diag2cens, death_cancer) ~ age_cent60 + 
  town5_Q1 + town5_Q2 + town5_Q4 + town5_Q5 +
  ethn3_Asian + ethn3_Black + 
  stage_late

mod_list[["loGI"]] <- fit_MI_mod_stpm2(imp_list, n=30, fm, df=4)

pool(mod_list[["loGI"]])

# model on complete data with average stage
dt <- gen_complete(imp_list)
fit <- stpm2(fm, data = dt, df=4)
mod_list2[["loGI"]] <- fit

## uter ----
imp_list <- do_MI(he_ova, cancer_type="uter", N =50)

fm <- Surv(fu_diag2cens, death_cancer) ~ age_cent60 + 
  town5_Q1 + town5_Q2 + town5_Q4 + town5_Q5 +
  ethn3_Asian + ethn3_Black + stage_late

mod_list[["uter"]] <- fit_MI_mod_stpm2(imp_list, n=30, fm, df=2)

pool(mod_list[["uter"]])

# model on complete data with average stage
dt <- gen_complete(imp_list)
fit <- stpm2(fm, data = dt, df=2)
mod_list2[["uter"]] <- fit

## panc----
imp_list <- do_MI(he_ova, cancer_type="panc", N=50)

fm <- Surv(fu_diag2cens, death_cancer) ~ age_cent60 + 
  town5_Q1 + town5_Q2 + town5_Q4 + town5_Q5 +
  ethn3_Asian + ethn3_Black +
  stage_late

mod_list[["panc"]] <- fit_MI_mod_stpm2(imp_list, n=30, fm, df=2)

pool(mod_list[["panc"]])

# model on complete data with average stage
dt <- gen_complete(imp_list)
fit <- stpm2(fm, data = dt, df=2)
mod_list2[["panc"]] <- fit

## other ----
imp_list <- do_MI(he_ova, cancer_type="other",N=50)

fm <- Surv(fu_diag2cens, death_cancer) ~ age_cent60 + 
  town5_Q1 + town5_Q2 + town5_Q4 + town5_Q5 +
  ethn3_Asian + ethn3_Black +
  stage_late

mod_list[["other"]] <- fit_MI_mod_stpm2(imp_list, n=30, fm, df=4)

pool(mod_list[["other"]])

# model on complete data with average stage
dt <- gen_complete(imp_list)
fit <- stpm2(fm, data = dt, df=4)
mod_list2[["other"]] <- fit

saveRDS(mod_list, file.path(work_data, "mod_list_MI_stpm2_20240910.rds"))
saveRDS(mod_list2, file.path(work_data, "mod_list_MI_takeMode_stpm2_20240910.rds"))

##> predict using the models----
ipd <- readRDS(file.path(work_data, "ipd2_20240906.RDS"))
mod_list <- readRDS(file.path(work_data, "mod_list_MI_stpm2_20240910.rds"))

# in case extra years of prediction is needed
max_time <- 15
n_cores=4
all <- c("ova", "lung", "loGI", "panc", "uter", "other")
surv_prob_list <- list()

for(c in all){
  surv_prob_list[[paste0(c, "_0")]] <- 
    pred_surv(mod_list, ipd=ipd, c_type=c, stage_late=0, max_time, n_cores)
  
  surv_prob_list[[paste0(c, "_1")]] <- 
    pred_surv(mod_list, ipd=ipd, c_type=c, stage_late=1, max_time, n_cores)
}

saveRDS(surv_prob_list, file.path(work_data, "surv_prob_list.rds"))


# Cost --------------------------------------------------------------------

rm(list = ls())

library(bannerCommenter)
library(mice)
library(miceadds)
library(tidyverse)
library(parallel)
library(matrixStats)
library(fastDummies)
library(lspline)
library(dtplyr)
library(rlang)
library(haven)
library(foreign)
library(foreach)
library(doSNOW)
library(snowfall)

source("//Qmcr/dfs/WIPH/Primary Care/OVATOOLS/r_project/scripts/file_path.R")

source(file.path(r_wd, "scripts/fn_MI.R"))
source(file.path(r_wd, "scripts/cost/fn_model_tests.R"))

he_ova <- readRDS(file.path(work_data, "he_ova_step4_withSymptom_prescr_NCRASupdated.RDS"))

# merge symptoms
he_ova <- merge_symptoms(he_ova)

cost_coef <- list()

## Ova ----

cancer_type="ova"

imp_list <- do_MI_4cost(he_ova, cancer_type=cancer_type, N=30)

var_x <- c("curr_age_cent", 
           "ethn3_Asian", "ethn3_Black", 
           "town5_Q1", "town5_Q2", "town5_Q4", "town5_Q5",
           # "stage_late", "admi_year_fct",
           "stage_late * admi_year_fct",
           "death_c", "death_nc", "fu_prop")

fit <- MI_2p_model(imp_list, var_x, name_mod="poi_log")

cost_coef[[cancer_type]][["p1"]] <- pool(fit$p1)$pooled[, "estimate"]
names(cost_coef[[cancer_type]][["p1"]]) <- pool(fit$p1)$pooled[, "term"]

cost_coef[[cancer_type]][["p2"]] <- pool(fit$p2)$pooled[, "estimate"]
names(cost_coef[[cancer_type]][["p2"]]) <- pool(fit$p2)$pooled[, "term"]

## loGI ----

cancer_type="loGI"

imp_list <- do_MI_4cost(he_ova, cancer_type=cancer_type, N=30)

var_x <- c("curr_age_cent", 
           "ethn3_Asian", "ethn3_Black", 
           "town5_Q1", "town5_Q2", "town5_Q4", "town5_Q5",
           # "stage_late", "admi_year_fct",
           "stage_late * admi_year_fct",
           "death_c", "death_nc", "fu_prop")

fit <- MI_2p_model(imp_list, var_x, name_mod="poi_log")

cost_coef[[cancer_type]][["p1"]] <- pool(fit$p1)$pooled[, "estimate"]
names(cost_coef[[cancer_type]][["p1"]]) <- pool(fit$p1)$pooled[, "term"]

cost_coef[[cancer_type]][["p2"]] <- pool(fit$p2)$pooled[, "estimate"]
names(cost_coef[[cancer_type]][["p2"]]) <- pool(fit$p2)$pooled[, "term"]

## lung ----

cancer_type="lung"

imp_list <- do_MI_4cost(he_ova, cancer_type=cancer_type, N=30)

var_x <- c("curr_age_cent", 
           "ethn3_Asian", "ethn3_Black", 
           "town5_Q1", "town5_Q2", "town5_Q4", "town5_Q5",
           # "stage_late", "admi_year_fct",
           "stage_late * admi_year_fct",
           "death_c", "death_nc", "fu_prop")

fit <- MI_2p_model(imp_list, var_x, name_mod="poi_log")

cost_coef[[cancer_type]][["p1"]] <- pool(fit$p1)$pooled[, "estimate"]
names(cost_coef[[cancer_type]][["p1"]]) <- pool(fit$p1)$pooled[, "term"]

cost_coef[[cancer_type]][["p2"]] <- pool(fit$p2)$pooled[, "estimate"]
names(cost_coef[[cancer_type]][["p2"]]) <- pool(fit$p2)$pooled[, "term"]

## panc ----

cancer_type="panc"

imp_list <- do_MI_4cost(he_ova, cancer_type=cancer_type, N=30)

var_x <- c("curr_age_cent", 
           "ethn3_Asian", "ethn3_Black", 
           "town5_Q1", "town5_Q2", "town5_Q4", "town5_Q5",
           # "stage_late", "admi_year_fct",
           "stage_late * admi_year_fct",
           "death_c", "death_nc", "fu_prop")

fit <- MI_2p_model(imp_list, var_x, name_mod="poi_log")

cost_coef[[cancer_type]][["p1"]] <- pool(fit$p1)$pooled[, "estimate"]
names(cost_coef[[cancer_type]][["p1"]]) <- pool(fit$p1)$pooled[, "term"]

cost_coef[[cancer_type]][["p2"]] <- pool(fit$p2)$pooled[, "estimate"]
names(cost_coef[[cancer_type]][["p2"]]) <- pool(fit$p2)$pooled[, "term"]

## uter ----

cancer_type="uter"

imp_list <- do_MI_4cost(he_ova, cancer_type=cancer_type, N=30)

var_x <- c("curr_age_cent", 
           "ethn3_Asian", "ethn3_Black", 
           "town5_Q1", "town5_Q2", "town5_Q4", "town5_Q5",
           # "stage_late", "admi_year_fct",
           "stage_late * admi_year_fct",
           "death_c", "death_nc", "fu_prop")

fit <- MI_2p_model(imp_list, var_x, name_mod="poi_log")

cost_coef[[cancer_type]][["p1"]] <- pool(fit$p1)$pooled[, "estimate"]
names(cost_coef[[cancer_type]][["p1"]]) <- pool(fit$p1)$pooled[, "term"]

cost_coef[[cancer_type]][["p2"]] <- pool(fit$p2)$pooled[, "estimate"]
names(cost_coef[[cancer_type]][["p2"]]) <- pool(fit$p2)$pooled[, "term"]

## other ----

cancer_type="other"

imp_list <- do_MI_4cost(he_ova, cancer_type=cancer_type, N=30)

var_x <- c("curr_age_cent", 
           "ethn3_Asian", "ethn3_Black", 
           "town5_Q1", "town5_Q2", "town5_Q4", "town5_Q5",
           # "stage_late", "admi_year_fct",
           "stage_late * admi_year_fct",
           "death_c", "death_nc", "fu_prop")

fit <- MI_2p_model(imp_list, var_x, name_mod="poi_log")

cost_coef[[cancer_type]][["p1"]] <- pool(fit$p1)$pooled[, "estimate"]
names(cost_coef[[cancer_type]][["p1"]]) <- pool(fit$p1)$pooled[, "term"]

cost_coef[[cancer_type]][["p2"]] <- pool(fit$p2)$pooled[, "estimate"]
names(cost_coef[[cancer_type]][["p2"]]) <- pool(fit$p2)$pooled[, "term"]

## benign ----
# because no stage data, no missing values anymore
cancer_type <- "benign"

d4s <- model_prep(cancer_type, he_ova)

d4s <- d4s %>% 
  dummy_cols(select_columns = c("town5","ethn3"))


var_x <- c("curr_age_cent", 
           "ethn3_Asian", "ethn3_Black", 
           "town5_Q1", "town5_Q2", "town5_Q4", "town5_Q5",
           "admi_year_fct",
           "death_c", "death_nc", "fu_prop")

fit <- list()
# Define the formula: outcome ~ covariate
var_y <- "cost01"

form <- as.formula(str_c(var_y, "~", str_c(var_x, collapse = " + ")))

# Fit the logistic regression model to the data 
fit$p1 <- glm(data = d4s, 
              formula = form, 
              family = binomial(link = "logit")) 

var_y <- "cost"

form <- as.formula(str_c(var_y, "~", str_c(var_x, collapse = " + ")))

fit$p2 <- glm(data = d4s %>% filter(cost > 0), 
              formula = form, 
              family = Gamma("log"))

cost_coef[[cancer_type]][["p1"]] <- fit$p1$coefficients
cost_coef[[cancer_type]][["p2"]] <- fit$p2$coefficients

##> modify for simulation ----
for (i in 1:length(cost_coef)){
  names(cost_coef[[i]]$p1)[1] <- "Intercept"
  names(cost_coef[[i]]$p2)[1] <- "Intercept"
  
  names(cost_coef[[i]]$p1)[2] <- "age_cent60" # cycle updated in simulation
  names(cost_coef[[i]]$p2)[2] <- "age_cent60" 
}

# saveRDS(cost_coef, file.path(work_data, "cost_coef.rds"))
# add stage year interaction
saveRDS(cost_coef, file.path(work_data, "cost_coef20241112.rds"))


# QoL ---------------------------------------------------------------------
dtqol <- readRDS(file.path(work_data, "qol_data.rds"))

dt <- dtqol %>% 
  drop_na(qol_score) %>% 
  rename(town5=townsend_score_no_na) %>% 
  rename(male=sex) %>% 
  mutate(
    gynae=case_when(cancer_gynae_3level=="No such cancer" ~ "0",
                    cancer_gynae_3level=="<=1 yr" ~ "1",
                    cancer_gynae_3level==">1 yr" ~ "2"
    ),
    gynae=factor(gynae, levels=c("0", "1", "2")),
    
    upGI=case_when(upGI_3level=="No such cancer" ~ "0",
                   upGI_3level=="<=1 yr" ~ "1",
                   upGI_3level==">1 yr" ~ "2"
    ),
    upGI=factor(upGI, levels=c("0", "1", "2")),
    
    lung=case_when(lung_3level=="No such cancer" ~ "0",
                   lung_3level=="<=1 yr" ~ "1",
                   lung_3level==">1 yr" ~ "2"
    ),
    lung=factor(lung, levels=c("0", "1", "2")),
    
    loGI=case_when(loGI_3level=="No such cancer" ~ "0",
                   loGI_3level=="<=1 yr" ~ "1",
                   loGI_3level==">1 yr" ~ "2"
    ),
    loGI=factor(loGI, levels=c("0", "1", "2")),
    
    cancers_other=case_when(cancers_other_2level=="No other cancer" ~ "0",
                            cancers_other_2level=="With other cancer" ~ "1"
    ),
    cancers_other=factor(cancers_other, levels=c("0", "1")),
    
    ben_gynae=case_when(ben_gynae_2level=="No such disease" ~ "0",
                        ben_gynae_2level=="With such disease" ~ "1"
    ),
    ben_gynae=factor(ben_gynae, levels=c("0", "1")),
  ) %>% 
  select(male, qol_score, age_cent60, ethn4, town5,  # use ethn4 including mixed
         gynae, upGI, lung, loGI, cancers_other, ben_gynae) 


dt <- dummy_cols(dt, select_columns = c("town5","ethn4", "gynae", "upGI", "lung", 
                                        "loGI", "cancers_other", "ben_gynae"),
                 remove_selected_columns = TRUE, remove_first_dummy = TRUE)


covar <- colnames(dt)[which(!colnames(dt) %in% c("qol_score", "age_cent60"))]
covar <- paste(covar, collapse = "+") 
# covar <- paste0(covar, "+male*cancers_other_1")


fm <- as.formula(paste0("qol_score ~ lspline(age_cent60, 1, marginal = F) +", covar))

mod_qol <- lm(fm, data = dt)

# summary(mod_qol)
# 
# lm.csv(mod_qol, output)

# add into cf

coef <- mod_qol$coefficients

# cf_b <- coef[which(!names(coef) %in% c("male", "male:cancers_other_1", "(Intercept)",
#                                       "lspline(age_cent60, 1, marginal = F)1",
#                                       "lspline(age_cent60, 1, marginal = F)2"))]

cf_b <- coef[which(names(coef) %in% c("town5_Q1","town5_Q2","town5_Q4","town5_Q5",
                                      "ethn4_Asian","ethn4_Black","ethn4_Others",
                                      "ben_gynae_1", "cancers_other_1"))]

cf_t <- coef[which(names(coef) %in% c("gynae_1","gynae_2","upGI_1", "upGI_2", 
                                      "lung_1", "lung_2", "loGI_1", "loGI_2"))]

cf_b <- c("Intercept"=as.numeric(coef["(Intercept)"]), cf_b)

# add late stage impact from literature
# Chase et al. Gynecologic Oncology 166 2022, 494-502
cf_b <- c(cf_b, stage_late=-0.046)

cf_t <- c(
  "age_cent60_1"=as.numeric(coef["lspline(age_cent60, 1, marginal = F)1"]),
  "age_cent60_2"=as.numeric(coef["lspline(age_cent60, 1, marginal = F)2"]),
  cf_t
)

qol <- list(
  cf_b = cf_b,
  cf_t = cf_t
)

# Create a cf file containing survival, cost and QoL ----------------------
surv_prob_list <- readRDS(file.path(work_data, "surv_prob_list.rds"))

# cost <- readRDS(file.path(work_data, "cost_coef.rds"))
cost <- readRDS(file.path(work_data, "cost_coef20241112.rds"))

cf_t_names <- c("age_cent60", "admi_year_fct1", "admi_year_fct2", 
                "admi_year_fct3", "admi_year_fct4", "admi_year_fct5", 
                "death_c", "death_nc", 
                "stage_late:admi_year_fct1", "stage_late:admi_year_fct2",
                "stage_late:admi_year_fct3", "stage_late:admi_year_fct4",
                "stage_late:admi_year_fct5"
                )

c_type <- names(cost)
for(c in c_type){
  x <- cost[[c]][["p1"]]
  y <- cost[[c]][["p2"]]
  
  cost[[c]][["p1"]] <- list(
    cf_b = x[setdiff(names(x), c(cf_t_names, "fu_prop"))],
    cf_t = x[cf_t_names]
  )
  cost[[c]][["p2"]] <- list(
    cf_b = y[setdiff(names(y), c(cf_t_names, "fu_prop"))],
    cf_t = y[cf_t_names]
  )
}

# remove NA as no interactions in benign
cost$benign$p1$cf_t <- cost$benign$p1$cf_t[!is.na(cost$benign$p1$cf_t)]
cost$benign$p2$cf_t <- cost$benign$p1$cf_t[!is.na(cost$benign$p1$cf_t)]
# use 0 to replace
repl <- c("stage_late:admi_year_fct1" = 0,
          "stage_late:admi_year_fct2" = 0,
          "stage_late:admi_year_fct3" = 0,
          "stage_late:admi_year_fct4" = 0,
          "stage_late:admi_year_fct5" = 0)

cost$benign$p1$cf_t <- c(cost$benign$p1$cf_t, repl)
cost$benign$p2$cf_t <- c(cost$benign$p2$cf_t, repl)

cf <- list(
  surv_prob_list = surv_prob_list,
  qol = qol,
  cost = cost
)

# saveRDS(cf, file.path(work_data, "cf_20240922.rds"))
# add stage, year interaction for cost models
saveRDS(cf, file.path(work_data, "cf_20241112.rds")) 

# Prepare mortality data --------------------------------------------------
library(readxl)
# read in the mortality data by sex and single age in UK
# in 2019, avoid influence from covid 19
cd_prob <- read_excel(file.path(work_data, "mortality_rate_UK_2019.xlsx"), 
                      sheet="cancer_death_prob") %>% 
  select(Age, Female) %>% 
  rename(cd = Female)

ncd_prob <- read_excel(file.path(work_data, "mortality_rate_UK_2019.xlsx"), 
                       sheet="noncancer_death_prob") %>% 
  select(Age, Female) %>% 
  rename(ncd = Female)


mort_prob <- data.frame(
  age=18:110L
)

mort_prob <- mort_prob %>% 
  mutate(
    age_cent60 = (age -60)/10,
    Age = case_when(
      age %in% 15:19 ~ "15-19",
      age %in% 20:24 ~ "20-24",
      age %in% 25:29 ~ "25-29",
      age %in% 30:34 ~ "30-34",
      age %in% 35:39 ~ "35-39",
      age %in% 40:44 ~ "40-44",
      age %in% 45:49 ~ "45-49",
      age %in% 50:54 ~ "50-54",
      age %in% 55:59 ~ "55-59",
      age %in% 60:64 ~ "60-64",
      age %in% 65:69 ~ "65-69",
      age %in% 70:74 ~ "70-74",
      age %in% 75:79 ~ "75-79",
      age %in% 80:84 ~ "80-84",
      age %in% 85:89 ~ "85-89",
      age >=90 ~ "90+"
    )
  ) %>%
  left_join(cd_prob) %>% 
  left_join(ncd_prob) %>% 
  select(-age, -Age) %>% 
  as.matrix()

saveRDS(mort_prob, file.path(work_data, "mort_prob.RDS"))

# Decision tree parameters ------------------------------------------------
# use data from Garth, for invasive ovarian cancer of overall model thresholds.

# overall 
# accuracy <- c(
# # for ovarian
# sens_ca125=0.849, 
# spec_ca125=0.937,
# sens_ovt1=0.902, # Funston et al. 2021
# spec_ovt1=0.880,
# sens_ovt3=0.816,
# spec_ovt3=0.966,
# sens_us= 0.849, # Menon et al. 2009 single US
# spec_us= 0.982,
# sens_us_fl_ca125= 0.894, # CA125 + US
# spec_us_fl_ca125= 0.998,
# sens_us_fl_ovt= 0.894, # copy from above CA125 + US
# spec_us_fl_ovt= 0.998,
# sens_ovt3_fl_us= 0.756, # copy from above
# spec_ovt3_fl_us= 0.969,
# # for other cancer
# # other gynaecological or lower abdominal cancers in false positive
# otc1_us_fp= 0.0134, # US, 5/373 Nyante 2011
# otc1_ovt3_fp= 0.1235, # ovatools 3%; this is CA125 35 (0.212-0.101)/(1-0.101) Funston et al. 2020
# otc1_ca125_us_fp= 0.0217, # ca125+/us+ Nyante 2011
# # lung or pancreas cancers in false positive
# otc2_us_fp= 0,
# otc2_ovt3_fp= 0.1235, # ovatools 3%; this is CA125 35 (0.212-0.101)/(1-0.101) Funston et al. 2020
# otc2_ca125_us_fp= 0
# )

# by stage accuracy for CA125 and Ovatools
# from Kirsten's latest data on 17/06/2024

accuracy <- list(
  # early stage
  accuracy_early = c(
    # for ovarian
    sens_ca125= 0.668, # 0.673,  
    spec_ca125= 0.935, # 0.933,
    sens_ovt1= 0.707, # 0.792, 
    spec_ovt1= 0.923, # 0.876,
    sens_ovt3= 0.548, # 0.584,
    spec_ovt3= 0.976, # 0.962,
    sens_us= 0.85, #0.615, #mixed stage  # 30/33, # Kalsi 2021 table 3 # 0.750, # Menon et al. 2009 single US
    spec_us= 0.83, #0.997, # 0.982,
    sens_us_fl_ca125= 0.85, #0.895, # 16/(16+3), # Menon et al. 2009 table 4, MMS; 0.895, # CA125 + US
    spec_us_fl_ca125= 0.83, #0.998,
    sens_us_fl_ovt= 0.85, #0.895, # 16/(16+3), # 0.895, # copy from above CA125 + US
    spec_us_fl_ovt= 0.83, #0.998,
    sens_ovt3_fl_us= 0.548, # 0.584, # copy from above
    spec_ovt3_fl_us= 0.976, # 0.962,
    
    # for other cancer
    # other gynaecological or lower abdominal cancers in false positive
    # otc1_us_fp= 0.0134, # US, 5/373 Nyante 2011
    # otc1_ovt3_fp= 0.1235, # ovatools 3%; this is CA125 35 (0.212-0.101)/(1-0.101) Funston et al. 2020
    # otc1_ca125_us_fp= 0.0217, # ca125+/us+ Nyante 2011
    
    # pathway involving US+
    fp_uter_us= (4+1)/(373+46), # Nyante 2011
    fp_loGI_us=1/(373+46), # Nyante 2011
    fp_uter_ovt3=(48/384)*0.1235, # ovatools 3%; this is CA125 35 (0.212-0.101)/(1-0.101) Funston et al. 2020
    fp_loGI_ovt3=(58/384)*0.1235, #Funston et al. 2020
    # fp_uter_ca125_us = 1/46, # Nyante 2011
    # fp_loGI_ca125_us = 0, # Nyante 2011
    # lung or pancreas cancers in false positive
    # otc2_us_fp= 0,
    # otc2_ovt3_fp= 0.1235, # ovatools 3%; this is CA125 35 (0.212-0.101)/(1-0.101) Funston et al. 2020
    # otc2_ca125_us_fp= 0,
    fp_lung_us=0,
    fp_panc_us=0,
    fp_lung_ovt3=(49/384)*0.1235, #Funston et al. 2020
    fp_panc_ovt3=(46/384)*0.1235, #Funston et al. 2020
    # fp_lung_ca125_us = 0, # Nyante 2011
    # fp_panc_ca125_us = 0, # Nyante 2011
    
    cost_us=204,
    cost_ca125=10,
    cost_gp1= 37, # face-to-face cons
    cost_gp2=16, # telephone cons after test 
    cost_nurse=9
  ),
  # late stage
  accuracy_late = c(
    # for ovarian
    sens_ca125= 0.933, # 0.947,  
    spec_ca125= 0.935, # 0.935,
    sens_ovt1= 0.946, # 0.971, 
    spec_ovt1= 0.923, # 0.878,
    sens_ovt3= 0.896, # 0.937,
    spec_ovt3=0.976, # 0.965,
    sens_us= 0.85, #0.615, #mixed stage # 50/97, # Kalsi 2021 table 3 # 0.750, # Menon et al. 2009 single US
    spec_us= 0.83, #0.997, # 0.982,
    sens_us_fl_ca125= 0.85, #0.895, # 18/(18+1), # Menon et al. 2009 table 4, MMS; # 0.895, # CA125 + US
    spec_us_fl_ca125= 0.83, #0.998,
    sens_us_fl_ovt= 0.85, #0.895, # 18/(18+1), # 0.895, # copy from above CA125 + US
    spec_us_fl_ovt= 0.83, #0.998,
    sens_ovt3_fl_us= 0.896, # 0.937, # copy from above
    spec_ovt3_fl_us= 0.976, # 0.965,
    
    # for other cancer
    # other gynaecological or lower abdominal cancers in false positive
    # otc1_us_fp= 0.0134, # US, 5/373 Nyante 2011
    # otc1_ovt3_fp= 0.1235, # ovatools 3%; this is CA125 35 (0.212-0.101)/(1-0.101) Funston et al. 2020
    # otc1_ca125_us_fp= 0.0217, # ca125+/us+ Nyante 2011
    fp_uter_us= (4+1)/(373+46), # Nyante 2011
    fp_loGI_us=1/(373+46), # Nyante 2011
    fp_uter_ovt3=(48/384)*0.1235, # ovatools 3%; this is CA125 35 (0.212-0.101)/(1-0.101) Funston et al. 2020
    fp_loGI_ovt3=(58/384)*0.1235, #Funston et al. 2020
    # fp_uter_ca125_us = 1/46, # Nyante 2011
    # fp_loGI_ca125_us = 0, # Nyante 2011
    # lung or pancreas cancers in false positive
    # otc2_us_fp= 0,
    # otc2_ovt3_fp= 0.1235, # ovatools 3%; this is CA125 35 (0.212-0.101)/(1-0.101) Funston et al. 2020
    # otc2_ca125_us_fp= 0,
    fp_lung_us=0,
    fp_panc_us=0,
    fp_lung_ovt3=(49/384)*0.1235, #Funston et al. 2020
    fp_panc_ovt3=(46/384)*0.1235, #Funston et al. 2020
    # fp_lung_ca125_us = 0, # Nyante 2011
    # fp_panc_ca125_us = 0, # Nyante 2011
    
    cost_us=204,
    cost_ca125=10,
    cost_gp1= 37, # face-to-face cons
    cost_gp2=16, # telephone cons after test 
    cost_nurse=9
  )
)

saveRDS(accuracy, file.path(work_data, "accuracy.RDS"))

# Stage shift parameters --------------------------------------------------
# proportional reduction of late stage at diagnosis in the faster detection group

stage_shift <- list(
  # ova, UKCTOCS, symptomatic, screen detected vs. all others
  ova_rr_late = 0.8364,
  ova_rr_early = 1.6194,
  # lung, PLCO
  lung_rr_late = 0.6288,
  lung_rr_early = 1.7109,
  # loGI/colorectal, a cohort study
  loGI_rr_late = 0.8532,
  loGI_rr_early = 1.0665,
  # panc, CAPS5
  panc_rr_late = 0.4761, # 0.3070,
  panc_rr_early = 1.9572, # 5.1579,
  
  # uter, no data temporarily, ajusted using dwelling time from ovarian rr
  uter_rr_late = 0.8831,
  uter_rr_early = 1.0449
)

saveRDS(stage_shift, file.path(work_data, "stage_shift.RDS"))


# false_positive_data -----------------------------------------------------

FP_data <- list(
  
  surg_prob = 0.7452, # data from ROCKet rapid-access sample, without any cancer, post-menopausl
  surg_cost = 3944,
  surg_disu = 0.04,
  surg_ben = 0.008,
  no_surg_cost = 181+204+19,
  ct_cost = 117# CT of one area, without contrast, 19 years and over NHS cost schedule 21/22
)

saveRDS(FP_data, file.path(work_data, "FP_data.RDS"))



