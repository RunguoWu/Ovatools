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

list_test <- list(gau_id = gaussian("identity"),
                  gau_log = gaussian("log"),
                  poi_id = poisson("identity"),
                  poi_log = poisson("log"),
                  gam_id = Gamma("identity"),
                  gam_log = Gamma("log"),
                  gam_sqrt = Gamma("sqrt")
)


cancer_type="ova"

# prepare the analytical data for test
# model_1p_2p will prepare the data automatically
d4s <- model_prep(cancer_type, he_ova)

# p1
var_x <- c("ethn3", "town5", "curr_age_cent", "stage_late*admi_year_fct",
           "death_c", "death_nc", "fu_prop")

var_y <- "cost01"
form <- as.formula(str_c(var_y, "~", str_c(var_x, collapse = " + ")))

p1 <- glm(data = d4s, 
                 formula = form, 
                 family = binomial(link = "logit")) 

#p2
ana <- d4s %>% filter(cost > 0)

var_y <- "cost"

form <- as.formula(str_c(var_y, "~", str_c(var_x, collapse = " + ")))

# > Gaussian - Identity ----
p2 <- glm(data = ana, 
                     formula = form, 
                     family = list_test[["poi_log"]])


summary(p1); summary(p2)