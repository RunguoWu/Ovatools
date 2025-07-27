##################################################################
##           Create PSA parameter sets for parameters           ##
##                 other than cancer death risk                 ##
##################################################################

library(bannerCommenter)
library(tidyverse)
library(fastDummies)
library(lspline)
library(matrixStats)
source("//Qmcr/dfs/WIPH/Primary Care/OVATOOLS/r_project/scripts/file_path.R")

source(file.path(scripts, "fn_data_clean.R"))

cf <- readRDS(file.path(work_data, "cf_20241112.rds"))

N <- 1000

# QoL ---------------------------------------------------------------------
dtqol <- readRDS(file.path(work_data, "qol_data.rds")) %>% 
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


dtqol <- dummy_cols(dtqol, select_columns = c("town5","ethn4", "gynae", "upGI", "lung", 
                                        "loGI", "cancers_other", "ben_gynae"),
                 remove_selected_columns = TRUE, remove_first_dummy = TRUE)


covar <- colnames(dtqol)[which(!colnames(dtqol) %in% c("qol_score", "age_cent60"))]
covar <- paste(covar, collapse = "+") 
# covar <- paste0(covar, "+male*cancers_other_1")

fm <- as.formula(paste0("qol_score ~ lspline(age_cent60, 1, marginal = F) +", covar))

set.seed(1234)
boots_ids <- replicate(N, sample(1:nrow(dtqol), nrow(dtqol), replace = TRUE))

qol_list <- list()
for (i in 1:N) {
  
  dt <- dtqol[boots_ids[, i], ]
  
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
  
  qol_list[[i]] <- qol
}

saveRDS(qol_list, file.path(work_data, "psa", "qol_list.rds"))

test_psa2(qol_list, cf)

## add the deterministic coefficient as 1001
qol_list[[1001]] <- cf$qol

saveRDS(qol_list, file.path(work_data, "psa", "qol_list.rds"))


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

N <- 1000
ncores <- 20

#> cancers----
var_x <- c("curr_age_cent", 
           "ethn3_Asian", "ethn3_Black", 
           "town5_Q1", "town5_Q2", "town5_Q4", "town5_Q5",
           # "stage_late", "admi_year_fct",
           "stage_late * admi_year_fct",
           "death_c", "death_nc", "fu_prop")

# cancer_type <- "ova"
cost_coef <- list()

for (cancer_type in c("panc")) {
  
  imp_list <- do_MI_4cost(he_ova, cancer_type=cancer_type, N=30)
  
  # data frames in imp_list are split by years
  # so need to sample in the patient ids with replace
  # and map the id into the original data frame to extract all split rows for the id
  # rows for one id could be extract more than once, 
  # depending on the frequency of the id in the sample
  # then combine them together to form the new boots data frame
  set.seed(1234)
  boots_ids <- replicate(N, sample(unique(imp_list[[1]]$epatid), length(unique(imp_list[[1]]$epatid)), replace = TRUE))

  sfInit(cpus = ncores, parallel = T)
  cl <- sfGetCluster()
  clusterExport(cl=cl, list = list())
  clusterEvalQ(cl, library(tidyverse))
  clusterEvalQ(cl, library(mice))

  registerDoSNOW(cl)
  
  cost_coef[[cancer_type]] <- foreach (i = 1:N) %dopar%{
    
    imp_list2 <- list()
    
    ids <- boots_ids[, i]
    
    # for efficiency, use the method mentioned above for the imp 1
    # extract id and admission year and map to other imps
    # because the two are unchanged across imps
    z <- do.call(rbind, lapply(ids, function(x) imp_list[[1]][imp_list[[1]]$epatid==x, ]))
    rownames(z) <- NULL
    
    imp_list2[[1]] <- z
    
    keys <- z %>% select(epatid, admi_year)
    
    for (imp in 2:length(imp_list)) {
      zz <- keys %>% left_join(imp_list[[imp]])
      imp_list2[[imp]] <- zz
    }

    # for (imp in 1:length(imp_list)) {
    #   imp_list2[[imp]] <- imp_list2[[imp]][boots_ids[, i], ]
    # }
    
    fit <- MI_2p_model(imp_list2, var_x, name_mod="poi_log")
    
    x <- list(
      p1 = pool(fit$p1)$pooled[, "estimate"],
      p2 = pool(fit$p2)$pooled[, "estimate"]
    )
    
    if(any(is.na(x$p1))) x$p1 <- pool_alter(fit$p1)[, "Estimate"]
    if(any(is.na(x$p2))) x$p2 <- pool_alter(fit$p2)[, "Estimate"]
    
    # i=636 stage_late:admi_year_fct5 =NA in all imputations
    # recode it as 0
    
    names(x$p1) <- pool(fit$p1)$pooled[, "term"]
    names(x$p2) <- pool(fit$p2)$pooled[, "term"]
    
    names(x$p1)[1] <- names(x$p2)[1] <- "Intercept"
    names(x$p1)[2] <- names(x$p2)[2] <- "age_cent60"
    
    # to avoid NA coefficient in any imputation
    # only happened to panc, 4-5 years after diagnosis

    
    return(x)
    # cost_coef[[cancer_type]][[i]] <- x
  }
  sfStop()

  # split into cf_b and cf_t
  cf_t_names <- c("age_cent60", "admi_year_fct1", "admi_year_fct2", 
                  "admi_year_fct3", "admi_year_fct4", "admi_year_fct5", 
                  "death_c", "death_nc",
                  "stage_late:admi_year_fct1", "stage_late:admi_year_fct2",
                  "stage_late:admi_year_fct3", "stage_late:admi_year_fct4",
                  "stage_late:admi_year_fct5"
                  )
  
  for (i in 1:N) {
    x <- cost_coef[[cancer_type]][[i]][["p1"]]
    y <- cost_coef[[cancer_type]][[i]][["p2"]]
    
    cost_coef[[cancer_type]][[i]][["p1"]] <- list(
      cf_b = x[setdiff(names(x), c(cf_t_names, "fu_prop"))],
      cf_t = x[cf_t_names]
    )
    cost_coef[[cancer_type]][[i]][["p2"]] <- list(
      cf_b = y[setdiff(names(y), c(cf_t_names, "fu_prop"))],
      cf_t = y[cf_t_names]
    )
  }

  # saveRDS(cost_coef, file.path(work_data, "psa", paste0("cost_coef_int_psa_", cancer_type, ".rds")))
  saveRDS(cost_coef, file.path(work_data, "psa", paste0("cost_coef_int_psa_", cancer_type, "_new.rds")))
  
}

#> benign----
cancer_type <- "benign"

d4s <- model_prep(cancer_type, he_ova)

d4s <- d4s %>% 
  dummy_cols(select_columns = c("town5","ethn3"))

var_x <- c("curr_age_cent", 
           "ethn3_Asian", "ethn3_Black", 
           "town5_Q1", "town5_Q2", "town5_Q4", "town5_Q5",
           "admi_year_fct",
           "death_c", "death_nc", "fu_prop")

# Define the formula: outcome ~ covariate
var_y <- "cost01"

form <- as.formula(str_c(var_y, "~", str_c(var_x, collapse = " + ")))

# Fit the logistic regression model to the data 
p1 <- glm(data = d4s, 
              formula = form, 
              family = binomial(link = "logit")) 

# Use Cholesky decomposition to generate parameter uncertainty
coef_est <- coef(p1)
vcov_mat <- vcov(p1)
chol_decomp <- chol(vcov_mat)

n_sim <- 1000

set.seed(5678)
random_err <- matrix(rnorm(n_sim * length(coef_est)), nrow = n_sim)
p1_psa <- t(apply(random_err, 1, function(x) coef_est + chol_decomp %*% x))
p1_psa <- as.data.frame(p1_psa)
colnames(p1_psa) <- names(coef_est)

# Fit the cost part model
var_y <- "cost"

form <- as.formula(str_c(var_y, "~", str_c(var_x, collapse = " + ")))

p2 <- glm(data = d4s %>% filter(cost > 0), 
              formula = form, 
              family = Gamma("log"))

coef_est <- coef(p2)
vcov_mat <- vcov(p2)
chol_decomp <- chol(vcov_mat)

random_err <- matrix(rnorm(n_sim * length(coef_est)), nrow = n_sim)

p2_psa <- t(apply(random_err, 1, function(x) coef_est + chol_decomp %*% x))
p2_psa <- as.data.frame(p2_psa)
colnames(p2_psa) <- names(coef_est)

colnames(p1_psa)[1] <- colnames(p2_psa)[1] <- "Intercept"
colnames(p1_psa)[2] <- colnames(p2_psa)[2] <- "age_cent60"

fit <- list()
cf_t_names <- c("age_cent60", "admi_year_fct1", "admi_year_fct2", 
                "admi_year_fct3", "admi_year_fct4", "admi_year_fct5", 
                "death_c", "death_nc")
for (i in 1: n_sim) {
  
  x <- as.numeric(p1_psa[i, ])
  y <- as.numeric(p2_psa[i, ])
  names(x) <- colnames(p1_psa)
  names(y) <- colnames(p2_psa)
  
  p <- list(
    p1 = list(
      cf_b = x[setdiff(names(x), c(cf_t_names, "fu_prop"))],
      cf_t = x[cf_t_names]
    ),
    p2 = list(
      cf_b = y[setdiff(names(y), c(cf_t_names, "fu_prop"))],
      cf_t = y[cf_t_names]
    )
  )
  
  repl <- c("stage_late:admi_year_fct1" = 0,
            "stage_late:admi_year_fct2" = 0,
            "stage_late:admi_year_fct3" = 0,
            "stage_late:admi_year_fct4" = 0,
            "stage_late:admi_year_fct5" = 0)
  p$p1$cf_t <- c(p$p1$cf_t, repl)
  p$p2$cf_t <- c(p$p2$cf_t, repl)
  
  fit[[i]] <- p
}

cost_coef2 <- readRDS(file.path(work_data, "psa", "cost_coef_int_psa_uter.rds"))

cost_coef2[["other"]] <- cost_coef[["other"]]

cost_coef2[[cancer_type]] <- fit

saveRDS(cost_coef2, file.path(work_data, "psa", "cost_coef_int_psa_all.rds"))


## add the deterministic coefficient as 1001
for (i in names(cost_coef_int_psa_all)) {
  cost_coef_int_psa_all[[i]][[1001]] <- cf$cost[[i]]
}

saveRDS(cost_coef_int_psa_all, file.path(work_data, "psa", "cost_coef_int_psa_all.rds"))

cancer_type <- "ova"
test_psa2(cost_coef_int_psa_all,  cancer_type)

# cost_coef$panc[[1001]] <- cf$cost$panc
# cost_coef_int_psa_panc_new$panc[[636]]$p2$cf_t["stage_late:admi_year_fct5"] <- 0
# saveRDS(cost_coef_int_psa_panc_new, file.path(work_data, "psa", paste0("cost_coef_int_psa_", cancer_type, "_new.rds")))
cost_coef_int_psa_all$panc <- cost_coef_int_psa_panc_new$panc
saveRDS(cost_coef_int_psa_all, file.path(work_data, "psa", "cost_coef_int_psa_all.rds"))

# Stage shift parameters --------------------------------------------------
# ova, UKCTOCS, symptomatic, screen detected vs. all others
set.seed(1234)

ova_rr_late <- exp(rnorm(1000, -0.17864, 0.0647449))
ova_rr_early <- exp(rnorm(1000, 0.482043, 0.1402807))

lung_rr_late <- exp(rnorm(1000, -0.463914335, 0.0695154))
lung_rr_early <- exp(rnorm(1000, 0.537014062, 0.0542976))

loGI_rr_late <- exp(rnorm(1000, -0.15878, 0.070747))
loGI_rr_early <- exp(rnorm(1000, 0.06442, 0.029856))
 
panc_rr_late <- exp(rnorm(1000, -0.742127, 0.294792))
panc_rr_early <- exp(rnorm(1000, 0.6715152, 0.13531))

uter_rr_late <- exp(rnorm(1000, -0.124260248, 0.0647449))
uter_rr_early <- exp(rnorm(1000, 0.0439478, 0.1402807))

N <- 1000
stage_shift <- list()
for (i in 1:N) {
  
  stage_shift[[i]] <- list(
    
    ova_rr_late = ova_rr_late[i], 
    ova_rr_early = ova_rr_early[i],

    lung_rr_late = lung_rr_late[i],
    lung_rr_early = lung_rr_early[i],

    loGI_rr_late = loGI_rr_late[i],
    loGI_rr_early = loGI_rr_early[i],

    panc_rr_late = panc_rr_late[i], 
    panc_rr_early = panc_rr_early[i], 
    
    uter_rr_late = uter_rr_late[i],
    uter_rr_early = uter_rr_early[i]
  )
}

# deterministic
stage_shift[[1001]] <- list(
  # ova, UKCTOCS, symptomatic, screen detected vs. all others
  ova_rr_late = 0.8364,
  ova_rr_early = 1.6194,
  # sensitivity analysis using upper and lower bounds
  ova_rr_late_lo = 0.7367,
  ova_rr_late_up = 0.9496,
  ova_rr_early_lo = 1.2301,
  ova_rr_early_up = 2.1318,
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

saveRDS(stage_shift, file.path(work_data, "psa", "stage_shift_psa.rds"))

test_psa4(stage_shift)


# False positive data -----------------------------------------------------
set.seed(1234)

surg_prob <- sample_beta2((78+155), (97+208))# data from ROCKet rapid-access sample, without any cancer, pre- and post-menopausl

surg_disu <- sample_norm(-0.04, -0.06, -0.01)

surg_ben <- c(rnorm(1000, 0.008, sqrt(0.084^2 + 0.087^2)/sqrt(311)), 0.008) 

N <- 1001
FP_data <- list()
for (i in 1:N) {
  
  FP_data[[i]] <- list(
    
    surg_prob = surg_prob[i], 
    surg_cost = 3944,
    surg_disu = surg_disu[i],
    surg_ben = surg_ben[i], 
    no_surg_cost = 181+204+19,
    ct_cost = 117# CT of one area, without contrast, 19 years and over NHS cost schedule 21/22
  )
}

saveRDS(FP_data, file.path(work_data, "psa", "FP_data_psa.rds"))

test_psa5(FP_data)

