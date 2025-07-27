##################################################################
##      Master of multi-state model using stpm2 prediction      ##
##                            For PSA                           ##
##################################################################

rm(list = ls())

library(tidyverse)
library(bannerCommenter)
library(data.table)
library(readxl)
library(parallel)
library(doSNOW)
library(rstpm2)

source("./github/file_path.R")
source(file.path(scripts, "fn_msm_stpm2.R"))

run_list <- list(
  list(type = "ova", value = 1),
  list(type = "ova", value = 0),
  list(type = "lung", value = 1),
  list(type = "lung", value = 0),
  list(type = "uter", value = 1),
  list(type = "uter", value = 0),
  list(type = "panc", value = 1),
  list(type = "panc", value = 0),
  list(type = "loGI", value = 1),
  list(type = "loGI", value = 0),
  list(type = "other", value = 1),
  list(type = "other", value = 0),
  list(type = "benign", value = 0)
)

for (s in 1:13) {
  
  c_type <- run_list[[s]][["type"]]
  stage_late <- run_list[[s]][["value"]]
  
  use_mort <- TRUE # TRUE means use UK mortality (CD/NCD) after 8 years
  
  n_cores <- 2
  
  # total PSA number
  n_sim <- 1000
  
  # read in the mortality data by sex and single age in UK
  mort_prob <- readRDS(file.path(work_data, "mort_prob.RDS"))
  
  age_stop <- expression(as.numeric(floor(110 - (vP_0[["age_cent60"]]*10 + 60) + 1)))
  
  # QoL coefficient PSA list
  qol_list <- readRDS(file.path(work_data, "qol_list.rds"))
  # cancer death probability PSA list 
  surv_prob <- if (c_type == "benign") NULL else 
    surv_prob <- readRDS(file.path(work_data, paste0("surv_prob_", c_type, "_all.rds")))
  
  # cost model coefficient PSA file
  cost_coef <- readRDS(file.path(work_data, paste0("cost_coef_int_psa_all.rds")))
  
  set.seed((s+100)*5) 
  # ensure each cancer type and stage has a separate random seed to generate the random number matrix
  mat_rand <- sample(1: n_sim^2, n_sim, replace = T)
  
  ptm <- proc.time()
  for (num in 1:n_sim) {
    
    # num <- 1
    
    # form the new cf file for each sim
    surv_prob_c <- list()
    surv_prob_c[[paste0(c_type, "_", stage_late)]] <- surv_prob[[paste0(c_type, "_", stage_late)]][[num]]
    
    cost_c <- list()
    cost_c[[c_type]] <- cost_coef[[c_type]][[num]]
    
    cf <- list(
      surv_prob_list = surv_prob_c,
      qol = qol_list[[num]],
      cost = cost_c
    )
    
    # output
    output_folder <- file.path(r_wd, "output", paste0(c_type, "_", stage_late))
    output_name <- paste0("psa_", c_type, "_stage", stage_late,"_useMort_", as.character(use_mort), "_", num)
    save_cycle <- F
    
    # patients
    IPD <- readRDS(file.path(work_data, "ipd2_pseudo.RDS")) %>% 
      # generate pseudo id
      mutate(id=1:n()) %>% 
      relocate(id, .before = age_cent60) %>% as.matrix()
    
    set.seed(mat_rand[num])
    # run model ---------------------------------------------------------------
    
    alpha <- master(
      IPD = IPD,
      mort_prob = mort_prob,
      age_stop = age_stop,
      n_cores = n_cores,
      output_folder = output_folder,
      output_name = output_name,
      cf = cf,
      c_type = c_type,
      stage_late = stage_late,
      use_mort = use_mort,
      save_cycle = save_cycle
    )
    
  }
  
  print(proc.time() - ptm)
  print(Sys.time())
  
}
