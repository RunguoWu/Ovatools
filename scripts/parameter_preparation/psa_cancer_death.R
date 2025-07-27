##################################################################
##                  Prepare PSA parameters sets                 ##
##                for cancer death probabilities                ##
##################################################################

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
ipd <- readRDS(file.path(work_data, "ipd2_20240906.RDS"))

mod_list <- list()
# surv_prob_list <- list()
# saveRDS(surv_prob_list, file.path(work_data, "surv_prob_list_psa.rds"))
max_time <- 10
n_cores=10

# random sample with replacement
N <- 1000

df_vec <- c(ova = 5, lung = 4, loGI = 4, uter= 2, panc = 2, other=4)

fm <- Surv(fu_diag2cens, death_cancer) ~ age_cent60 + 
  town5_Q1 + town5_Q2 + town5_Q4 + town5_Q5 +
  ethn3_Asian + ethn3_Black + 
  stage_late

for(c_type in c("other")){
  
  # c_type <- "ova"
  # get imp_list
  
  df <- df_vec[c_type]
  
  imp_list <- do_MI(he_ova, cancer_type=c_type, N=50)
  
  set.seed(1234)
  boots_ids <- replicate(N, sample(1:nrow(imp_list[[1]]), nrow(imp_list[[1]]), replace = TRUE))
  
  set.seed(1014)
  boots_ids2 <- replicate(N, sample(1:nrow(imp_list[[1]]), nrow(imp_list[[1]]), replace = TRUE))
  
  set.seed(1017)
  boots_ids3 <- replicate(N, sample(1:nrow(imp_list[[1]]), nrow(imp_list[[1]]), replace = TRUE))
  
  for(i in 101:248){
    
    ptm <- proc.time()
    
    surv_prob_list <- list()
    
    imp_list2 <- imp_list
    
    for (imp in 1:length(imp_list)) {
      imp_list2[[imp]] <- imp_list2[[imp]][boots_ids[, i], ]
    }
    
    mod_list[[c_type]] <- fit_MI_mod_stpm2(imp_list2, n=30, fm, df=df)
    
    # try second boots ids
    if (is.null(mod_list[[c_type]][[1]])) {
      
      # try the backup bootstrap ids
      imp_list2 <- imp_list
      
      for (imp in 1:length(imp_list)) {
        imp_list2[[imp]] <- imp_list2[[imp]][boots_ids2[, i], ]
      }
      
      mod_list[[c_type]] <- fit_MI_mod_stpm2(imp_list2, n=30, fm, df=df)
    }
    
    # try third boots ids
    if (is.null(mod_list[[c_type]][[1]])) {
      
      # try the backup bootstrap ids
      imp_list2 <- imp_list
      
      for (imp in 1:length(imp_list)) {
        imp_list2[[imp]] <- imp_list2[[imp]][boots_ids3[, i], ]
      }
      
      mod_list[[c_type]] <- fit_MI_mod_stpm2(imp_list2, n=30, fm, df=df)
    }
    
    surv_prob_list[[paste0(c_type, "_0")]][[i]] <- 
      pred_surv(mod_list, ipd=ipd, c_type=c_type, stage_late=0, max_time, n_cores)
    
    surv_prob_list[[paste0(c_type, "_1")]][[i]] <- 
      pred_surv(mod_list, ipd=ipd, c_type=c_type, stage_late=1, max_time, n_cores)
    
    saveRDS(surv_prob_list, file.path(work_data, "psa", paste0("surv_prob_list_", c_type, "_", i, ".rds")))
    
    print(proc.time() - ptm)
    print(Sys.time())
  }
  
}





# Combine predicted survival probabilities --------------------------------
# to combine predicted survival probabilities for psa

source("//Qmcr/dfs/WIPH/Primary Care/OVATOOLS/r_project/scripts/file_path.R")
library(tidyverse)

# at beginning, some files were saved in the single file
surv_prob_list_psa <- readRDS(file.path(work_data, paste0("surv_prob_list_psa.rds")))

cf <- readRDS(file.path(work_data, "cf_20240922.RDS"))

combine_psa <- function(surv_prob_list_psa, c_type, cf){
  
  surv_prob <- list()
  
  surv_prob[[paste0(c_type, "_0")]] <- surv_prob_list_psa[[paste0(c_type, "_0")]]
  surv_prob[[paste0(c_type, "_1")]] <- surv_prob_list_psa[[paste0(c_type, "_1")]]
  
  for(i in (length(surv_prob_list_psa[[paste0(c_type, "_0")]])+1):1000){
    
    x <- readRDS(file.path(work_data, "psa", paste0("surv_prob_list_", c_type, "_",i, ".rds")))
    
    surv_prob[[paste0(c_type, "_0")]][[i]] <- x[[paste0(c_type, "_0")]][[i]]
    surv_prob[[paste0(c_type, "_1")]][[i]] <- x[[paste0(c_type, "_1")]][[i]]
  }
  
  # add the deterministic one as 1001
  surv_prob[[paste0(c_type, "_0")]][[1001]] <- cf$surv_prob_list[[paste0(c_type, "_0")]]
  surv_prob[[paste0(c_type, "_1")]][[1001]] <- cf$surv_prob_list[[paste0(c_type, "_1")]]
  
  saveRDS(surv_prob, file.path(work_data, "psa", paste0("surv_prob_", c_type, "_all.rds")))
}

combine_psa(surv_prob_list_psa, c_type="other", cf)


# test if CI is around the mean
test_psa <- function(surv_prob_psa, c_type, if_late=0){
  
  det <- surv_prob_psa[[paste0(c_type, "_", if_late)]][[1001]]
  det_avg <- colMeans(det)
  psa_mat <- matrix(nrow = 1000, ncol = ncol(surv_prob_psa[[paste0(c_type, "_", if_late)]][[1]]))
  
  for (i in 1:1000){
    psa_avg <- colMeans(surv_prob_psa[[paste0(c_type, "_", if_late)]][[i]])
    psa_mat[i, ] <- psa_avg
  }
  
  library(matrixStats)
  CI <- colQuantiles(psa_mat, probs = c(0.025, 0.975))
  
  to_plot <- data.frame(
    year = 1:10,
    mean = det_avg[1:10],
    ci_l = CI[,1],
    ci_h = CI[,2]
  )
  
  p <- ggplot(to_plot, aes(x = year, y = mean))+
    geom_line()+
    geom_errorbar(aes(ymin = ci_l, ymax = ci_h), width=0.2)+
    theme_minimal()
  
  return(p)
}

c_type="other"

surv_prob_psa <- readRDS(file.path(work_data, "psa", paste0( "surv_prob_", c_type, "_all.rds")))

test_psa(surv_prob_psa, c_type=c_type, if_late=0)





