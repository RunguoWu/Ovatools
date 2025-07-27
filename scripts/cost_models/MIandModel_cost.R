##################################################################
##              Multiple Imputation for cost model              ##
##                       & cost modelling                       ##
##################################################################

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

# Ova ---------------------------------------------------------------------

cancer_type="ova"

imp_list <- do_MI_4cost(he_ova, cancer_type=cancer_type)

var_x <- c("curr_age_cent", 
           "ethn3_Asian", "ethn3_Black", 
           "town5_Q1", "town5_Q2", "town5_Q4", "town5_Q5",
           "stage_late", "admi_year_fct",
           "death_c", "death_nc", "fu_prop")

fit <- MI_2p_model(imp_list, var_x, name_mod="poi_log")

cost_coef[[cancer_type]][["p1"]] <- pool(fit$p1)$pooled[, "estimate"]
names(cost_coef[[cancer_type]][["p1"]]) <- pool(fit$p1)$pooled[, "term"]

cost_coef[[cancer_type]][["p2"]] <- pool(fit$p2)$pooled[, "estimate"]
names(cost_coef[[cancer_type]][["p2"]]) <- pool(fit$p2)$pooled[, "term"]

# loGI --------------------------------------------------------------------

cancer_type="loGI"

imp_list <- do_MI_4cost(he_ova, cancer_type=cancer_type)

var_x <- c("curr_age_cent", 
           "ethn3_Asian", "ethn3_Black", 
           "town5_Q1", "town5_Q2", "town5_Q4", "town5_Q5",
           "stage_late", "admi_year_fct",
           "death_c", "death_nc", "fu_prop")

fit <- MI_2p_model(imp_list, var_x, name_mod="poi_log")

cost_coef[[cancer_type]][["p1"]] <- pool(fit$p1)$pooled[, "estimate"]
names(cost_coef[[cancer_type]][["p1"]]) <- pool(fit$p1)$pooled[, "term"]

cost_coef[[cancer_type]][["p2"]] <- pool(fit$p2)$pooled[, "estimate"]
names(cost_coef[[cancer_type]][["p2"]]) <- pool(fit$p2)$pooled[, "term"]

# lung --------------------------------------------------------------------

cancer_type="lung"

imp_list <- do_MI_4cost(he_ova, cancer_type=cancer_type)

var_x <- c("curr_age_cent", 
           "ethn3_Asian", "ethn3_Black", 
           "town5_Q1", "town5_Q2", "town5_Q4", "town5_Q5",
           "stage_late", "admi_year_fct",
           "death_c", "death_nc", "fu_prop")

fit <- MI_2p_model(imp_list, var_x, name_mod="poi_log")

cost_coef[[cancer_type]][["p1"]] <- pool(fit$p1)$pooled[, "estimate"]
names(cost_coef[[cancer_type]][["p1"]]) <- pool(fit$p1)$pooled[, "term"]

cost_coef[[cancer_type]][["p2"]] <- pool(fit$p2)$pooled[, "estimate"]
names(cost_coef[[cancer_type]][["p2"]]) <- pool(fit$p2)$pooled[, "term"]

# panc --------------------------------------------------------------------

cancer_type="panc"

imp_list <- do_MI_4cost(he_ova, cancer_type=cancer_type)

var_x <- c("curr_age_cent", 
           "ethn3_Asian", "ethn3_Black", 
           "town5_Q1", "town5_Q2", "town5_Q4", "town5_Q5",
           "stage_late", "admi_year_fct",
           "death_c", "death_nc", "fu_prop")

fit <- MI_2p_model(imp_list, var_x, name_mod="poi_log")

cost_coef[[cancer_type]][["p1"]] <- pool(fit$p1)$pooled[, "estimate"]
names(cost_coef[[cancer_type]][["p1"]]) <- pool(fit$p1)$pooled[, "term"]

cost_coef[[cancer_type]][["p2"]] <- pool(fit$p2)$pooled[, "estimate"]
names(cost_coef[[cancer_type]][["p2"]]) <- pool(fit$p2)$pooled[, "term"]

# uter --------------------------------------------------------------------

cancer_type="uter"

imp_list <- do_MI_4cost(he_ova, cancer_type=cancer_type)

var_x <- c("curr_age_cent", 
           "ethn3_Asian", "ethn3_Black", 
           "town5_Q1", "town5_Q2", "town5_Q4", "town5_Q5",
           "stage_late", "admi_year_fct",
           "death_c", "death_nc", "fu_prop")

fit <- MI_2p_model(imp_list, var_x, name_mod="poi_log")

cost_coef[[cancer_type]][["p1"]] <- pool(fit$p1)$pooled[, "estimate"]
names(cost_coef[[cancer_type]][["p1"]]) <- pool(fit$p1)$pooled[, "term"]

cost_coef[[cancer_type]][["p2"]] <- pool(fit$p2)$pooled[, "estimate"]
names(cost_coef[[cancer_type]][["p2"]]) <- pool(fit$p2)$pooled[, "term"]

# other -------------------------------------------------------------------

cancer_type="other"

imp_list <- do_MI_4cost(he_ova, cancer_type=cancer_type)

var_x <- c("curr_age_cent", 
           "ethn3_Asian", "ethn3_Black", 
           "town5_Q1", "town5_Q2", "town5_Q4", "town5_Q5",
           "stage_late", "admi_year_fct",
           "death_c", "death_nc", "fu_prop")

fit <- MI_2p_model(imp_list, var_x, name_mod="poi_log")

cost_coef[[cancer_type]][["p1"]] <- pool(fit$p1)$pooled[, "estimate"]
names(cost_coef[[cancer_type]][["p1"]]) <- pool(fit$p1)$pooled[, "term"]

cost_coef[[cancer_type]][["p2"]] <- pool(fit$p2)$pooled[, "estimate"]
names(cost_coef[[cancer_type]][["p2"]]) <- pool(fit$p2)$pooled[, "term"]

# benign ------------------------------------------------------------------
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


for (i in 1:length(cost_coef)){
  names(cost_coef[[i]]$p1)[1] <- "Intercept"
  names(cost_coef[[i]]$p2)[1] <- "Intercept"
  
  names(cost_coef[[i]]$p1)[2] <- "age_cent60" # cycle updated in simulation
  names(cost_coef[[i]]$p2)[2] <- "age_cent60" 
}

saveRDS(cost_coef, file.path(work_data, "cost_coef.rds"))


