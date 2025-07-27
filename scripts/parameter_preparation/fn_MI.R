#################################################################
##                Multiple Imputation functions                ##
#################################################################

library(mice)
library(miceadds)
library(tidyverse)
library(parallel)
library(matrixStats)
# library(flexsurv)
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


# Merge in symptoms&prescribe
merge_symptoms <- function(he_ova){
  
  # read in wide-form symptom data
  # prepared in fn_data_clean.R
  # because the original symptom data is too large, I split it and generate two datasets
  symptoms_wide1 <- readRDS(file.path(r_wd, "working_data", "symptoms_prescr_wide1.RDS"))
  symptoms_wide2 <- readRDS(file.path(r_wd, "working_data", "symptoms_prescr_wide2.RDS"))
  
  # only keep records within 3 months before and after the index date
  
  dt1 <- he_ova %>% 
    inner_join(symptoms_wide1)%>% 
    lazy_dt() %>% 
    # calculate the gap between US date and symptom date
    # only keep ids with at least one gap within 3 months before US date
    filter(if_any(starts_with("symptom_obs_"), ~ as.numeric(index_date - .x)>-93 & 
                    as.numeric(index_date - .x)<93)) %>% 
    mutate(symptom = 1) %>% 
    select(-starts_with("symptom_obs_")) %>%
    select(epatid, symptom) %>% 
    as.data.frame()
  
  dt2 <- he_ova %>% 
    inner_join(symptoms_wide2)%>% 
    lazy_dt() %>% 
    # calculate the gap between US date and symptom date
    # only keep ids with at least one gap within 3 months before US date
    filter(if_any(starts_with("symptom_obs_"), ~ as.numeric(index_date - .x)>-93 & 
                    as.numeric(index_date - .x)<93)) %>% 
    mutate(symptom = 1) %>% 
    select(-starts_with("symptom_obs_")) %>%
    select(epatid, symptom) %>% 
    as.data.frame()
  
  dt <- bind_rows(dt1, dt2)
  
  he_ova <- he_ova %>%
    left_join(dt, by = "epatid") %>% 
    mutate(symptom = if_else(is.na(symptom), 0, symptom))
  
  return(he_ova)
}

# cancer_type <- "ova"
# Multiple imputation
# cancer_type: "ova", "lung", "panc", "loGI", "uter", "other"
do_MI <- function(he_ova, cancer_type, N, rand_seed=12345, for_cost=F){

  if (for_cost){
    filter_name <- paste0("cancer_", cancer_type, "_icd")
    stage_name <- paste0("stage_", cancer_type, "_allicd")
    date_name <- paste0(cancer_type, "_icd_date")
  } else {
    filter_name <- paste0("cancer_", cancer_type, "_1yr")
    stage_name <- paste0("stage_", cancer_type, "_imp")
    date_name <- paste0(cancer_type, "_1yr_date")
  }
  
  behaviour_name <- paste0("behaviour_", cancer_type)
  
  if(cancer_type=="other") behaviour_name <- "behaviour_all"
  
  # select vars for imputation
  dt_imp <- he_ova %>% 
    filter(!!sym(filter_name)==1) %>% 
    mutate(fu_diag2cens = as.numeric(censor_date - !!sym(date_name))/365.25,
           fu_diag2cens = if_else(fu_diag2cens<=0, 1/365.25, fu_diag2cens)
           ) %>% # for other cancer, one case = -1/365.25, so use fu_diag2cens<=0 here
    select(epatid, fu_diag2cens, death_cancer, death_other, age_cent60, 
           town5, ethn3, ethn4, CA125,symptom, all_of(c(behaviour_name, stage_name))
    ) %>% 
    mutate(!!sym(behaviour_name) := case_when(!!sym(behaviour_name)==1 ~ "malignant",
                                 !!sym(behaviour_name)==0 ~ "borderline",
                                 !!sym(behaviour_name)==2~ "secondary",
                                 !!sym(behaviour_name)==3 ~ "unknown",
                                 is.na(!!sym(behaviour_name))~ "unknown"))
  
  # simplify auxiliary variables
  dt_imp[[stage_name]] <- factor(dt_imp[[stage_name]], ordered = T)
  dt_imp$death_cancer <- as.factor(dt_imp$death_cancer)
  dt_imp$death_other <- as.factor(dt_imp$death_other)
  dt_imp$town5 <- factor(dt_imp$town5, ordered = T)
  dt_imp$ethn3 <- factor(dt_imp$ethn3, ordered = F)
  dt_imp$symptom <- factor(dt_imp$symptom)
  dt_imp[[behaviour_name]] <- factor(dt_imp[[behaviour_name]], ordered = F)
  
  # pre-define the matrix ----
  imp <- mice(dt_imp, maxit=0)
  predM = imp$predictorMatrix
  meth = imp$method

  # ID and ethn4 won't predict other var.
  predM[, c("epatid", "ethn4", "CA125")] <- 0
  predM["CA125", ] <- 0 # do impute CA125
  meth["CA125"] <- ""
  
  
  # start impute ----
  imp <- futuremice(data = dt_imp, m=N, n.core = 4, parallelseed = rand_seed,
                    predictorMatrix = predM, method = meth)
  
  imp_list <- miceadds::mids2datlist(imp)
  
  # add age group
  add_age <- he_ova %>% 
    filter(!!sym(filter_name)==1) %>% 
    select(epatid, age, age_group10, age_group_40_80, age_group_50_80, age_group_60_80)
  
  for (i in 1:N) {
    imp_list[[i]] <- imp_list[[i]] %>% 
      mutate(stage_late = 
               case_when( 
                 !!sym(stage_name) %in% c("stage 1", "stage 2") ~ 0, 
                 !!sym(stage_name) %in% c("stage 3", "stage 4") ~ 1
               )
      ) %>% 
      dummy_cols(select_columns = c("town5","ethn3"),
                 remove_selected_columns = FALSE, remove_first_dummy = TRUE) %>% 
      mutate(death_cancer = if_else(death_cancer=="1", 1, 0)) %>% 
      # add age group
      left_join(add_age)
  }
  
  return(imp_list)
}

# output is the original mice output
do_MI2 <- function(he_ova, cancer_type, N, rand_seed=12345, for_cost=F){
  
  if (for_cost){
    filter_name <- paste0("cancer_", cancer_type, "_icd")
    stage_name <- paste0("stage_", cancer_type, "_allicd")
    date_name <- paste0(cancer_type, "_icd_date")
  } else {
    filter_name <- paste0("cancer_", cancer_type, "_1yr")
    stage_name <- paste0("stage_", cancer_type, "_imp")
    date_name <- paste0(cancer_type, "_1yr_date")
  }
  
  behaviour_name <- paste0("behaviour_", cancer_type)
  
  if(cancer_type=="other") behaviour_name <- "behaviour_all"
  
  # select vars for imputation
  dt_imp <- he_ova %>% 
    filter(!!sym(filter_name)==1) %>% 
    mutate(fu_diag2cens = as.numeric(censor_date - !!sym(date_name))/365.25,
           fu_diag2cens = if_else(fu_diag2cens<=0, 1/365.25, fu_diag2cens)
    ) %>% # for other cancer, one case = -1/365.25, so use fu_diag2cens<=0 here
    select(epatid, fu_diag2cens, death_cancer, death_other, age_cent60, 
           town5, ethn3, ethn4, symptom, all_of(c(behaviour_name, stage_name))
    ) %>% 
    mutate(!!sym(behaviour_name) := case_when(!!sym(behaviour_name)==1 ~ "malignant",
                                              !!sym(behaviour_name)==0 ~ "borderline",
                                              !!sym(behaviour_name)==2~ "secondary",
                                              !!sym(behaviour_name)==3 ~ "unknown",
                                              is.na(!!sym(behaviour_name))~ "unknown"))
  
  # simplify auxiliary variables
  dt_imp[[stage_name]] <- factor(dt_imp[[stage_name]], ordered = T)
  dt_imp$death_cancer <- as.factor(dt_imp$death_cancer)
  dt_imp$death_other <- as.factor(dt_imp$death_other)
  dt_imp$town5 <- factor(dt_imp$town5, ordered = T)
  dt_imp$ethn3 <- factor(dt_imp$ethn3, ordered = F)
  dt_imp$symptom <- factor(dt_imp$symptom)
  dt_imp[[behaviour_name]] <- factor(dt_imp[[behaviour_name]], ordered = F)
  
  # pre-define the matrix ----
  imp <- mice(dt_imp, maxit=0)
  predM = imp$predictorMatrix
  meth = imp$method
  
  # ID and ethn4 won't predict other var.
  predM[, c("epatid", "ethn4")]=0
  
  # start impute ----
  imp <- futuremice(data = dt_imp, m=N, n.core = 4, parallelseed = rand_seed,
                    predictorMatrix = predM, method = meth)
  
  imp2 <- complete(imp, "long", include = T)
  
  # add age group
  add_age <- he_ova %>% 
    filter(!!sym(filter_name)==1) %>% 
    select(epatid, age, age_group10, age_group_40_80, age_group_50_80, age_group_60_80)
  
  imp2 <- imp2 %>% mutate(stage_late = 
                            case_when( 
                              !!sym(stage_name) %in% c("stage 1", "stage 2") ~ 0, 
                              !!sym(stage_name) %in% c("stage 3", "stage 4") ~ 1
                            )) %>%
    dummy_cols(select_columns = c("town5","ethn3"),
               remove_selected_columns = FALSE, remove_first_dummy = TRUE) %>% 
    mutate(death_cancer = if_else(death_cancer=="1", 1, 0)) %>% 
    # add age group
    left_join(add_age)
  
  imp2 <- as.mids(imp2)
  
  return(imp2)
}





# This function is to address failure in some imputations when fitting flexible models
# it will fit more models than needed, and go on if error happen, 
# finally select the front n successful models
# imp_list: output from do_MI
# n: the number of successful models needed
# k: k in flexible parametric model
# scale: scale in flexible parametric model

fit_MI_mod <- function(imp_list, n, fm, k, scale){
  
  sfInit(cpus = 4, parallel = T)
  cl <- sfGetCluster()
  clusterExport(cl=cl, list = list())
  clusterEvalQ(cl, library(tidyverse))
  clusterEvalQ(cl, library(survival))
  clusterEvalQ(cl, library(flexsurv))
  clusterEvalQ(cl, library(mice))
  registerDoSNOW(cl)
  
  mods <- foreach(i = c(1:length(imp_list))) %dopar%{ 
    
    # fit model
    fit <- tryCatch({
      flexsurvspline(fm, data = imp_list[[i]], k = k, scale = scale)
    }, error = function(e){
      return(NULL)
    })
    
    return(fit)
  }
  sfStop()
  
  mods2 <- mods[!sapply(mods, is.null)]
  
  rt <- mods2[1:min(n, length(mods2))]
  
  return(rt)
}

# use stpm2
fit_MI_mod_stpm2 <- function(imp_list, n, fm, df){
  
  sfInit(cpus = 32, parallel = T)
  cl <- sfGetCluster()
  clusterExport(cl=cl, list = list())
  clusterEvalQ(cl, library(tidyverse))
  clusterEvalQ(cl, library(survival))
  clusterEvalQ(cl, library(rstpm2))
  clusterEvalQ(cl, library(mice))
  registerDoSNOW(cl)
  
  mods <- foreach(i = c(1:length(imp_list))) %dopar%{ 
    
    # fit model
    fit <- tryCatch({
      stpm2(fm, data = imp_list[[i]], df=df)
    }, error = function(e){
      return(NULL)
    })
    
    return(fit)
  }
  sfStop()
  
  mods2 <- mods[!sapply(mods, is.null)]
  
  rt <- mods2[1:min(n, length(mods2))]
  
  return(rt)
}



# generate a complete data set with stage_late imputed as the mode of all imputations
gen_complete <- function(imp_list){
  
  x <- sapply(imp_list, "[[", "stage_late")
  
  complete_case <- imp_list[[1]] %>% 
    mutate(mean_stage = rowMeans(x)) %>% 
    mutate(stage_late = if_else(mean_stage>=0.5, 1, 0))
  
  return(complete_case)
}


# prepare for cost model
do_MI_4cost <- function(he_ova, cancer_type, N, rand_seed=12345){
  
  imp_list <- do_MI(he_ova, cancer_type=cancer_type,N=N, rand_seed=rand_seed, for_cost=T)
  
  rt <- list()
  
  for (i in 1:length(imp_list)) {
    
    dt <- imp_list[[i]]
    
    rt[[i]] <- readRDS(file.path(work_data, paste0("cost_", cancer_type, ".rds"))) %>% 
      left_join(dt, by= "epatid") %>% 
      # Convert cost outcome to 1 or 0 for the part 1 model
      mutate(cost01 = if_else(cost > 0, 1, 0), 
             diag_time = as.numeric(diag_date - index_date),
             # convert into factor, where 5 means 5+
             admi_year_fct = as.factor(if_else(admi_year>=5, 5, admi_year)),
             curr_age_cent=(curr_age-60)/10
      ) %>% 
      # specify death reason
      mutate(death_c = if_else(death == 1 & death_cancer==1, 1, 0),
             death_nc = if_else(death == 1 & death_cancer==0, 1, 0)
      )
  }
  
  return(rt)
}


# alternative way to get combined coefficients and variance across multiple imputations
# e.g. mod_list_i <- mod_list[["loGI"]]
pool_alter <- function(mod_list_i){
  
  estimates <- lapply(mod_list_i, function(model) coef(model))
  variances <- lapply(mod_list_i, function(model) vcov(model))
  
  m <- length(mod_list_i)
  Q_bar <- colMeans(do.call(rbind, estimates), na.rm = TRUE)
  U_bar <- Reduce("+", variances)/m
  B <- cov(do.call(rbind, estimates))
  T <- U_bar + (1 + 1/m) * B
  se <- sqrt(diag(T))
  ci_lower <- Q_bar - qt(0.975, df = Inf)*se
  ci_upper <- Q_bar + qt(0.975, df = Inf)*se
  results <- data.frame(Estimate = Q_bar, se, ci_lower, ci_upper)
  return(results)
}

# Fit two-part model based on MI data
# imp_list: list of multiply imputed data, output of do_MI_4cost
# var_x: independent variables
# name_mod: the model type of the second part

MI_2p_model <- function(imp_list, var_x, name_mod){
  
  list_dist <- list(gau_id = gaussian("identity"),
                    gau_log = gaussian("log"),
                    poi_id = poisson("identity"),
                    poi_log = poisson("log"),
                    gam_id = Gamma("identity"),
                    gam_log = Gamma("log")#,
                    # gam_sqrt = Gamma("sqrt")
  )
  
  fit <- list(p1 = list(), p2 = list())
  for (i in 1:length(imp_list)){
    # Part one
    var_y <- "cost01"
    form <- as.formula(str_c(var_y, "~", str_c(var_x, collapse = " + ")))
    fit[["p1"]][[i]] <- glm(data =imp_list[[i]], form, family = binomial(link = "logit"))
    
    # Part two
    var_y <- "cost"
    form <- as.formula(str_c(var_y, "~", str_c(var_x, collapse = " + ")))
    fit[["p2"]][[i]] <- glm(data =imp_list[[i]] %>% 
                              filter(cost > 0) %>% 
                              mutate(cost = round(cost, 0)), 
                            form, family = list_dist[[name_mod]])
  }
  
  return(fit)
}

# print MI risk models with est and se
print_coef <- function(mod_list, cancer_type){
  
  rt <- tryCatch({
    x <- pool(mod_list[[cancer_type]])
    est <- round(x$pooled[["estimate"]], 2)
    se <- round(sqrt(x$pooled[["t"]]), 2)
    nam <- as.character(x$pooled[["term"]])
    list(est=est, se=se, nam=nam)
  }, error = function(e){
    x <- pool_alter(mod_list[[cancer_type]])
    est <- round(x$Estimate, 2)
    se <- round(x$se, 2)
    nam <- rownames(x)
    list(est=est, se=se, nam=nam)
  })
  
  coef <- paste0(rt$est, " (", rt$se, ")")
  names(coef) <- rt$nam
  return(coef)
}

# print ordinary model with est and se
print_coef2 <- function(mod_list, cancer_type){
  
  x <- mod_list[[cancer_type]]
  est <- round(x$coefficients, 2)
  se <- round(sqrt(diag(x$cov)), 2)
  coef <- paste0(est, " (", se, ")")
  names(coef) <- names(est)
  return(coef)
}

# print stpm2 model with est and se
print_coef_stpm2 <- function(mod_list, cancer_type){
  
  x <- mod_list[[cancer_type]]
  x <- summary(x)@coef
  est <- round(x[, "Estimate"], 2)
  se <- round(x[, "Std. Error"], 2)
  coef <- paste0(est, " (", se, ")")
  names(coef) <- names(est)
  return(coef)
}

# print cost model: two part model
print_coef_cost <- function(fit, benign=FALSE){ # fit from MI_2p_model
  
  if (benign) {
    x <- summary(fit$p1)
    est <- round(x$coefficients[, "Estimate"], 2)
    se <- round(x$coefficients[, "Std. Error"], 2)
    nam <- as.character(rownames(x$coefficients))
    rt1 <- list(est=est, se=se, nam=nam)
    
    x <- summary(fit$p2)
    est <- round(x$coefficients[, "Estimate"], 2)
    se <- round(x$coefficients[, "Std. Error"], 2)
    nam <- as.character(rownames(x$coefficients))
    rt2 <- list(est=est, se=se, nam=nam)
  } else {
    x <- pool(fit$p1)
    est <- round(x$pooled[["estimate"]], 2)
    se <- round(sqrt(x$pooled[["t"]]), 2)
    nam <- as.character(x$pooled[["term"]])
    rt1 <- list(est=est, se=se, nam=nam)
    
    x <- pool(fit$p2)
    est <- round(x$pooled[["estimate"]], 2)
    se <- round(sqrt(x$pooled[["t"]]), 2)
    nam <- as.character(x$pooled[["term"]])
    rt2 <- list(est=est, se=se, nam=nam)
  }

  coef1 <- paste0(rt1$est, " (", rt1$se, ")")
  names(coef1) <- rt1$nam
  
  coef2 <- paste0(rt2$est, " (", rt2$se, ")")
  names(coef2) <- rt2$nam
  
  coef <- cbind(coef1, coef2)
  
  return(coef)
}

print_glm_model <- function(fit){
  
  x <- summary(fit)
  est <- round(x$coefficients[, "Estimate"], 2)
  se <- round(x$coefficients[, "Std. Error"], 2)
  nam <- as.character(rownames(x$coefficients))
  rt <- paste0(est, " (", se, ")")
  names(rt) <- nam
  return(rt)
}

plot_stpm2 <- function(fit, by_stage=0, max_time=8, d4s){
  
  time_seq <- 0:max_time
  
  if (by_stage==0){
    pred_mat <- matrix(NA, nrow = length(time_seq), ncol = nrow(d4s))
    
    for(i in 1:length(time_seq)){
      newdata <- d4s %>% mutate(fu_diag2cens=time_seq[i])
      
      pred_mat[i, ] <- predict(fit, newdata=newdata)
    }
    avg_pred <- rowMeans(pred_mat)
    
    km <- survfit(Surv(fu_diag2cens, death_cancer) ~ 1, d4s)
    plot(km, col="black")
    lines(time_seq, avg_pred, col="red")  
  } else {
    d4s0 <- subset(d4s, stage_late==0)
    d4s1 <- subset(d4s, stage_late==1)
    
    pred_mat0 <- matrix(NA, nrow = length(time_seq), ncol = nrow(d4s0))
    for(i in 1:length(time_seq)){
      newdata <- d4s0 %>% mutate(fu_diag2cens=time_seq[i])
      pred_mat0[i, ] <- predict(fit, newdata=newdata)
    }
    avg_pred0 <- rowMeans(pred_mat0)
    
    pred_mat1 <- matrix(NA, nrow = length(time_seq), ncol = nrow(d4s1))
    for(i in 1:length(time_seq)){
      newdata <- d4s1 %>% mutate(fu_diag2cens=time_seq[i])
      pred_mat1[i, ] <- predict(fit, newdata=newdata)
    }
    avg_pred1 <- rowMeans(pred_mat1)
    
    km0 <- survfit(Surv(fu_diag2cens, death_cancer) ~ 1, d4s0)
    km1 <- survfit(Surv(fu_diag2cens, death_cancer) ~ 1, d4s1)
    plot(km0, col="grey")
    lines(km1, col="black")
    lines(time_seq, avg_pred0, col="blue")  
    lines(time_seq, avg_pred1, col="red") 
  }
}

# c_type: ova, lung, loGI, panc, uter, other
# stage_late: 0, 1
# predict survival probability until 110 years of age
pred_surv <- function(mod_list, ipd, c_type, stage_late, max_time, n_cores){
  
  mod <- mod_list[[c_type]]
  
  ipd2 <- ipd %>% mutate(stage_late=stage_late)
  
  cl <- makeSOCKcluster(n_cores)
  clusterExport(cl, c("c_type", "max_time"))
  clusterEvalQ(cl, library(tidyverse))
  clusterEvalQ(cl, library(rstpm2))
  registerDoSNOW(cl)
  
  ts <- 1:max_time
  
  surv_prob <- parLapply(cl = cl, ts, function(j){
    
    ipd2 <- ipd2 %>% mutate(fu_diag2cens = j)
    
    pred_mat <- matrix(NA, nrow = nrow(ipd2), ncol = length(mod))
    for (i in 1:length(mod)) {
      pred <- predict(mod[[i]], ipd2)
      pred_mat[, i] <- pred
    }
    avg_pred <- rowMeans(pred_mat)
    
    return(avg_pred)
  }) 
  stopCluster(cl)
  
  ret <- do.call("cbind", surv_prob)
  colnames(ret) <- paste0("year_", 1:max_time)
  
  return(ret)
}


# reduce the size using the actual simulation years
# however, it does not reduce the size in effect...
reduce_size <- function(surv_prob, ipd){
  
  sim_years <- as.numeric(floor(110 - (ipd$age_cent60*10 + 60) + 1))
  
  red_surv_prob <- matrix(NA, nrow = nrow(surv_prob), ncol = ncol(surv_prob))
  
  for (i in 1:nrow(surv_prob)) {
    red_surv_prob[i, ][1: sim_years[i]] <- surv_prob[i, ][1: sim_years[i]]
  }
  colnames(red_surv_prob) <- colnames(surv_prob)
  
  return(red_surv_prob)
}







