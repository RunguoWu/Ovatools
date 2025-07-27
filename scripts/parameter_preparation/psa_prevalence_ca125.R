##################################################################
##                Prepare PSA set for prevalence                ##
##################################################################

# copy from psa_prevalence.R
# Just develop the models on patients with CA125 records
# but use the models to predict prevalence and stage for all
# just pick up those with CA125 in post-simulation analysis step

rm(list = ls())

source("//Qmcr/dfs/WIPH/Primary Care/OVATOOLS/r_project/scripts/file_path.R")
source(file.path(r_wd, "scripts/fn_data_clean.R"))
library(bannerCommenter)
library(fastDummies)
library(tidyverse)
library(nnet)

source(file.path(r_wd, "scripts/fn_MI.R"))
he_ova <- readRDS(file.path(work_data, "he_ova_step4_withSymptom_prescr_NCRASupdated.RDS"))
# merge symptoms for MI
he_ova <- merge_symptoms(he_ova)

dt0 <- he_ova %>% 
  mutate(
    ben_gynae_1 = if_else(b_benign_gynae==1 |benign_1yr==1, 1, 0)
  ) %>% 
  select(epatid, age_cent60, ben_gynae_1, age_group_40_80, age_group_50_80, age_group_60_80, 
         ethn3, ethn4, town5, b_uter, b_loGI, b_lung, b_panc, b_cancer_others, ben_gynae_1,
         cancer_ova_1yr, cancer_uter_1yr, cancer_loGI_1yr, cancer_lung_1yr, 
         cancer_panc_1yr, cancer_other_1yr,
         stage_ova_imp, stage_loGI_imp, stage_uter_imp, stage_lung_imp,
         stage_panc_imp, stage_other_imp, CA125
  ) 

dt1 <- dt0 %>% filter(!is.na(CA125))


# it is too slow to use multinomial model and MI to estimate the PSA set
# estimate the prevalence first without stage information 
# could avoid MI
# based on the MI result of cancer stage data, and use beta distribution to sample 1000 
N <- 1000
ncores <- 20
#> prevalence: all----

c_list <- c("lung", "loGI", "uter", "panc")

fm_vec <- c(ova = as.formula("cancer_ova_1yr ~ age_group_40_80 + ethn4 + town5"), 
            lung = as.formula("cancer_lung_1yr ~ age_group_60_80 + ethn3 + town5"), 
            loGI = as.formula("cancer_loGI_1yr ~ age_group_50_80 + ethn3 + town5"), 
            uter= as.formula("cancer_uter_1yr ~ age_group_50_80 + ethn3 + town5"), # ethn4 have trouble with others
            panc = as.formula("cancer_panc_1yr ~ age_group_60_80 + ethn3 + town5"), 
            other= as.formula("cancer_other_1yr ~ age_group_40_80 + ethn4 + town5")
)

# cancer_type <- "lung"

for (cancer_type in c_list) {
  
  fm <- fm_vec[[cancer_type]]
  
  set.seed(1234)
  boots_ids <- replicate(N, sample(1:nrow(dt1), nrow(dt1), replace = TRUE))
  
  sfInit(cpus = ncores, parallel = T)
  cl <- sfGetCluster()
  # clusterExport(cl=cl, list = list())
  clusterEvalQ(cl, library(tidyverse))
  clusterEvalQ(cl, library(nnet))
  clusterEvalQ(cl, library(mice))
  registerDoSNOW(cl)
  
  pred <- foreach (i = 1:N) %dopar%{
    
    dt2 <- dt1[boots_ids[, i], ]
    
    if(cancer_type %in% c("loGI", "panc", "uter", "lung")) {
      dt2 <- dt2[dt2[[paste0("b_", cancer_type)]]==0, ]
    }
    if(cancer_type =="other" ) {
      dt2 <- dt2[dt2[["b_cancer_others"]]==0, ]
    }
    
    fit <- glm(fm, dt2, family="binomial") 
    pred <- predict(fit, dt0, type = "response")
    
    return(pred)
  }
  sfStop()
  
  # deterministic prediction
  
  if(cancer_type %in% c("loGI", "panc", "uter", "lung")) {
    dt2 <- dt1[dt1[[paste0("b_", cancer_type)]]==0, ]
  } else if(cancer_type =="other" ) {
    dt2 <- dt1[dt1[["b_cancer_others"]]==0, ]
  } else {
    dt2 <- dt1}
  
  fit <- glm(fm, dt2, family="binomial")
  pred_det <-  predict(fit, dt0, type = "response")
  pred[[length(pred) + 1]] <- pred_det
  saveRDS(pred, file.path(work_data, "psa/prev", paste0(cancer_type, "_prev_ca125.rds")))
} 

# print the deterministic results
for (cancer_type in c("ova", "loGI", "panc", "uter", "lung", "other")) {
  fm <- fm_vec[[cancer_type]]
  
  if(cancer_type %in% c("loGI", "panc", "uter", "lung")) {
    dt2 <- dt1[dt1[[paste0("b_", cancer_type)]]==0, ]
  } else if(cancer_type =="other" ) {
    dt2 <- dt1[dt1[["b_cancer_others"]]==0, ]
  } else {
    dt2 <- dt1}
  
  fit <- glm(fm, dt2, family="binomial")
  write.csv(print_glm_model(fit), file.path(output, paste0(cancer_type, "_prevMod_ca125.csv")))
}


## check distribution
library(matrixStats)
psa <- pred[-1001]

psa_pred <- do.call(cbind, psa)
ci <- rowQuantiles(psa_pred, probs = c(0.025, 0.975))

det <- pred[[1001]]

to_plot <- data.frame(
  mean = det,
  ci_l = ci[,1],
  ci_h = ci[,2]
) %>%
  arrange(mean) %>%
  distinct() %>%
  mutate(var = 1:n())

ggplot(to_plot, aes(x = var, y = mean))+
  geom_line()+
  geom_errorbar(aes(ymin = ci_l, ymax = ci_h), width=0.2)+
  theme_minimal()


#> stage: on cancer cases only----
N <- 1000
ncores <- 20
c_list <- c("ova", "lung", "loGI", "uter", "panc")

fm_vec <- c(ova = as.formula("stage_late ~ age_group_40_80 + ethn4 + town5"), 
            lung = as.formula("stage_late ~ age_group_60_80 + ethn3 + town5"), # no late stage in other ethnicity 
            loGI = as.formula("stage_late ~ age_group_50_80 + ethn3 + town5"), 
            uter= as.formula("stage_late ~ age_group_50_80 + ethn3 + town5"), 
            panc = as.formula("stage_late ~ age_group_60_80 + ethn3 + town5"), 
            other= as.formula("stage_late ~ age_group_40_80 + ethn4 + town5")
)

for (cancer_type in c_list) {
  
  pred <- list()
  
  fm <- fm_vec[[cancer_type]]
  
  imp_list <- do_MI(he_ova, cancer_type=cancer_type, N=30)
  
  # for those with CA125
  for (i in 1:length(imp_list)) {
    imp_list[[i]] <- imp_list[[i]] %>% filter(!is.na(CA125))
  }
  
  # data frames in imp_list are split by years
  # so need to sample in the patient ids with replace
  # and map the id into the original data frame to extract all split rows for the id
  # rows for one id could be extract more than once, 
  # depending on the frequency of the id in the sample
  # then combine them together to form the new boots data frame
  set.seed(1234)
  boots_ids <- replicate(N, sample(1:nrow(imp_list[[1]]), nrow(imp_list[[1]]), replace = TRUE))
  
  set.seed(1014)
  boots_ids2 <- replicate(N, sample(1:nrow(imp_list[[1]]), nrow(imp_list[[1]]), replace = TRUE))
  
  set.seed(1017)
  boots_ids3 <- replicate(N, sample(1:nrow(imp_list[[1]]), nrow(imp_list[[1]]), replace = TRUE))
  
  set.seed(2019)
  boots_ids4 <- replicate(N, sample(1:nrow(imp_list[[1]]), nrow(imp_list[[1]]), replace = TRUE))
  
  set.seed(2022)
  boots_ids5 <- replicate(N, sample(1:nrow(imp_list[[1]]), nrow(imp_list[[1]]), replace = TRUE))
  
  set.seed(3022)
  boots_ids6 <- replicate(N, sample(1:nrow(imp_list[[1]]), nrow(imp_list[[1]]), replace = TRUE))
  
  set.seed(876)
  boots_ids7 <- replicate(N, sample(1:nrow(imp_list[[1]]), nrow(imp_list[[1]]), replace = TRUE))
  
  set.seed(992)
  boots_ids8 <- replicate(N, sample(1:nrow(imp_list[[1]]), nrow(imp_list[[1]]), replace = TRUE))
  
  set.seed(12)
  boots_ids9 <- replicate(N, sample(1:nrow(imp_list[[1]]), nrow(imp_list[[1]]), replace = TRUE))
  
  set.seed(67895)
  boots_ids10 <- replicate(N, sample(1:nrow(imp_list[[1]]), nrow(imp_list[[1]]), replace = TRUE))
  
  sfInit(cpus = ncores, parallel = T)
  cl <- sfGetCluster()
  clusterExport(cl=cl, list = list("boots_ids2", "boots_ids3", "boots_ids4", 
                                   "boots_ids5", "boots_ids6", "boots_ids7" , 
                                   "boots_ids8", "boots_ids9", "boots_ids10" ))
  clusterEvalQ(cl, library(tidyverse))
  clusterEvalQ(cl, library(mice))
  
  registerDoSNOW(cl)
  
  pred <- foreach (i = 1:N) %dopar%{ # too large
    
    pred_imp <- c()
    
    for (imp in 1:length(imp_list)) {
      
      pred_imp_i <- tryCatch({
        ids <- boots_ids[, i]
        dtt <- imp_list[[imp]][ids, ]
        fit <- glm(fm, dtt, family="binomial")
        predict(fit, type = "response", dt0)
      }, error = function(e){
        return(NULL)
      })
      
      for (x in 2:10){
        
        if(is.null(pred_imp_i)){
          pred_imp_i <- tryCatch({
            ids <- get(paste0("boots_ids", x))[, i]
            dtt <- imp_list[[imp]][ids, ]
            fit <- glm(fm, dtt, family="binomial")
            predict(fit, type = "response", dt0)
          }, error = function(e){
            return(NULL)
          })
        } else break
      }
      
      pred_imp <- cbind(pred_imp, pred_imp_i)
      
    }
    pred_avg <- rowMeans(pred_imp)
    
    return(pred_avg)
    
  }
  sfStop()
  
  # deterministic prediction
  pred_imp <- c()
  for (imp in 1:length(imp_list)) {
    fit <- glm(fm, imp_list[[imp]], family="binomial")
    pred_imp <- cbind(pred_imp, predict(fit, type = "response", dt0))
  }
  pred_avg <- rowMeans(pred_imp)
  pred[[length(pred) + 1]] <- pred_avg
  
  saveRDS(pred, file.path(work_data, "psa/stage_rate", paste0(cancer_type, "_stage_ca125.rds")))
  
}

# print the deterministic results
mod_list <- list()
for (cancer_type in c("ova", "loGI", "panc", "uter", "lung", "other")) {
  fm <- fm_vec[[cancer_type]]
  imp_list <- do_MI(he_ova, cancer_type=cancer_type, N=30)
  
  # for those with CA125
  for (i in 1:length(imp_list)) {
    imp_list[[i]] <- imp_list[[i]] %>% filter(!is.na(CA125))
  }
  
  mod_list[[cancer_type]] <- list()
  for (imp in 1:length(imp_list)) {
    fit <- glm(fm, imp_list[[imp]], family="binomial")
    mod_list[[cancer_type]][[imp]] <- fit
  }  
  
  x <- pool_alter(mod_list[[cancer_type]])
  est <- round(x$Estimate, 2)
  se <- round(x$se, 2)
  nam <- rownames(x)
  
  coef <- paste0(est, " (", se, ")")
  names(coef) <- nam
  
  write.csv(coef, file.path(output, paste0(cancer_type, "_stageMod_ca125.csv")))
}


## check distribution
psa <- ova_stage[-1001]

psa_pred <- do.call(cbind, psa)
ci <- rowQuantiles(psa_pred, probs = c(0.025, 0.975))

det <- ova_stage[[1001]]

to_plot <- data.frame(
  mean = det,
  ci_l = ci[,1],
  ci_h = ci[,2]
) %>%
  arrange(mean) %>%
  distinct() %>%
  mutate(var = 1:n())

ggplot(to_plot, aes(x = var, y = mean))+
  geom_line()+
  geom_errorbar(aes(ymin = ci_l, ymax = ci_h), width=0.2)+
  theme_minimal()

#> save in csv----
library(foreach)
library(snowfall)

source("//Qmcr/dfs/WIPH/Primary Care/OVATOOLS/r_project/scripts/file_path.R")

# ova
wd_prev <- file.path(work_data, "psa/prev")
wd_stage <- file.path(work_data, "psa/stage_rate")

prev <- readRDS(file.path(wd_prev, "ova_prev_ca125.rds"))
stage <- readRDS(file.path(wd_stage, "ova_stage_ca125.rds"))

sfInit(cpus = 10, parallel = T)
cl <- sfGetCluster()
clusterExport(cl=cl, list = list())
registerDoSNOW(cl)

foreach(i = 1:length(prev))%dopar%{
  
  x <- data.frame(prevalence = prev[[i]])
  saveRDS(x, file.path(wd_prev, "ova", paste0("ova_prev_ca125", i, ".rds")))
  
  y <- data.frame(ova_stage34 = stage[[i]])
  saveRDS(y, file.path(wd_stage, "ova", paste0("ova_stage_ca125", i, ".rds")))
  
}
sfStop()

# for others

for (c_type in c("lung", "panc", "uter", "loGI", "other")) {
  
  if(c_type=="other") base <- "b_cancer_others" else
    base <- paste0("b_", c_type)
  
  dtt <- dt0 %>% select(!!sym(base))
  stage <- readRDS(file.path(wd_stage, paste0(c_type, "_stage_ca125.rds")))
  prev <- readRDS(file.path(wd_prev, paste0(c_type, "_prev_ca125.rds")))
  
  sfInit(cpus = 10, parallel = T)
  cl <- sfGetCluster()
  clusterExport(cl=cl, list = list())
  clusterEvalQ(cl, library(tidyverse))
  registerDoSNOW(cl)
  
  foreach(i = 1:length(stage))%dopar%{
    
    x <- prev[[i]]
    
    if(c_type=="other") base <- "b_cancer_others" else
      base <- paste0("b_", c_type)
    
    rt <- cbind(dtt, x) %>% 
      mutate(!!sym(paste0("pred_", c_type)) := if_else(!!sym(base)==1, 0, x)) %>% 
      select(!!sym(paste0("pred_", c_type)))
    
    y <- stage[[i]]
    
    rt2 <- cbind(dtt, y) %>% 
      mutate(!!sym(paste0(c_type, "_stage34")) := if_else(!!sym(base)==1, 0, y)) %>% 
      select(!!sym(paste0(c_type, "_stage34")))
    
    saveRDS(rt, file.path(wd_prev, c_type, paste0(c_type, "_prev_ca125", i, ".rds")))
    saveRDS(rt2, file.path(wd_stage, c_type, paste0(c_type, "_stage_ca125", i, ".rds")))
    
  }
  sfStop()
}






