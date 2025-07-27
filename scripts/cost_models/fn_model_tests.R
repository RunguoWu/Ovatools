#################################################################
##                       model test codes                      ##
##                  adapted from Tutorial code 4               ##
#################################################################

# Source packages
library(tidyverse)
library(sandwich) # f_clx
library(lmtest) # f_clx
library(aod)  # f_clx
select <- dplyr::select

# Prepare function ----

# > Cluster robust standard error ----

f_clx <- function(mod,cluster){
  # library(sandwich)
  # library(lmtest)
  
  ## https://rdrr.io/cran/ivpack/src/R/clx.R
  ## https://www.ne.su.se/polopoly_fs/1.216115.1426234213!/menu/standard/file/clustering1.pdf
  dfcw=1
  M <- length(unique(cluster))
  N <- length(cluster)
  dfc <- (M/(M-1))*((N-1)/(N-mod$rank))
  u <- apply(estfun(mod),2,
             function(x) tapply(x, cluster, sum))
  vcovCL <- dfc*sandwich(mod, meat.=crossprod(u)/N)*dfcw
  test <- coeftest(mod, vcovCL)
  list(vcovCL = vcovCL, coeftest = test)
}

# > Specification test ----

f_test_spe_hl_clx <- function(epatid, y_hat, res){
  
  # Required package
  # library(sandwich) # f_clx
  # library(lmtest) # f_clx
  # library(aod)  # f_clx
  
  # Required function: f_clx
  
  # - checks fit on raw scale for systematic bias 
  ana <- tibble(epatid = epatid, 
                res = res)
  ana$y_hat_decile <- factor(dplyr::ntile(y_hat, n = 10), levels = 1:10)
  hl_mod <- lm(res ~ y_hat_decile - 1, data = ana)
  hl_test <- aod::wald.test(b = coef(hl_mod), Sigma = f_clx(hl_mod, ana$epatid)$vcovCL, Terms = 1:10)
  output <- tibble(hl_p = hl_test$result[[1]][3])
  return(output)
}

f_test_spe_pl_clx <- function(family, epatid, y, y_link){
  
  # - determines linearity of response 
  ana <- tibble(y = y,
                y_link = y_link,
                epatid = epatid)
  pl_mod <- glm(formula = y ~ + y_link + I(y_link^2), data = ana, 
                family = family, start = c(0,1,0))
  pl_test <- f_clx(pl_mod, ana$epatid)$coeftest
  output <- tibble(pl_p = pl_test[3,4], pl_est = pl_test[3,1])
  return(output)  
}

f_test_spe_mp_clx <- function(epatid, y_hat, res){
  
  # Modified Park test
  # - determines family 
  ana <- tibble(y_hat = y_hat, 
                res = res, 
                epatid = epatid)
  n_le0 <- sum(ana$y_hat <= 0)
  
  output <- tibble(mp_slope = NA)
  if(n_le0 == 0){ # it will not work if there is any negative predicted value
    mp_mod <- glm(I(res^2) ~ I(log(y_hat)), data = ana %>% filter(res != 0), 
                  family = Gamma(link="log"), start = rep(2, 2))
    mp_test <- f_clx(mp_mod, ana %>% filter(res !=0) %>% pull(epatid))
    
    output$mp_slope <- mp_test$coeftest[2]
    # Slope
    # 0: Gaussian
    # 1: Poisson
    # 2: Gamma
    # 3: Inverse Gaussian or Wald
  }
  return(output)  
}

# Performance tests ----

## Prepare function ----

f_test_gof <- function(res){
  
  # ME: mean error
  # MAE: mean absolute error
  # RMSE: root mean squared error
  output <- tibble(me = round(mean(res)),
                   mae = round(mean(abs(res))), 
                   rmse = round(sqrt(mean(res^2))))
  return(output)
  
}

## prepare data for analysis
model_prep <- function(cancer_type, he_ova){
  
  if (cancer_type != "benign") {
    # use newly created he_ova_step4_withSymptom_prescr_NCRASupdated.RDS
    stage_name <- paste0("stage_", cancer_type, "_allicd")
    
    dt <- he_ova %>% 
      mutate(stage_late = 
               case_when( 
                 !!sym(stage_name) %in% c("stage 1", "stage 2") ~ 0, 
                 !!sym(stage_name) %in% c("stage 3", "stage 4") ~ 1,
                 is.na(!!sym(stage_name)) ~ NA
               )
      ) %>% 
      select(epatid, e_pracid, death_cancer, death_other, age_cent60, town5, ethn3, 
             stage_late)
  } else {
    
    dt <- he_ova %>% 
      select(epatid, e_pracid, death_cancer, death_other, age_cent60, town5, ethn3) %>% 
      mutate(stage_late = 0)
  }
  
  d4s <- readRDS(file.path(work_data, paste0("cost_", cancer_type, ".rds"))) %>% 
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
    ) %>% 
    filter(!is.na(stage_late))
  
  return(d4s)
}    
 

# Fit one part and two part models ----------------------------------------
# fit a one-part model and a two-part model, using the same model type

# cancer_type = "ova", "lung", "panc", "loGI", "uter", "benign"
# mod_type is one type in the list above
# he_ova is the study population
# name_mod <- c("gau_id", "gau_log", "poi_id", "poi_log", "gam_id", "gam_log", "gam_sqrt") 
# above: the model type to be tested for the second part of 2PM

model_1p_2p <- function(cancer_type, he_ova, var_x, name_mod){
  
  list_test <- list(gau_id = gaussian("identity"),
                    gau_log = gaussian("log"),
                    poi_id = poisson("identity"),
                    poi_log = poisson("log"),
                    gam_id = Gamma("identity"),
                    gam_log = Gamma("log"),
                    gam_sqrt = Gamma("sqrt")
  )
  # do not fit poi_id, gam_id and gam_sqrt, because they will incur starting value issue
  # the possible reason is the very negative coefficient of fu_prop. 
  
  d4s <- model_prep(cancer_type, he_ova)
    
  # One-part model  ----
  
  # Model the costs incurring in the period directly 
  
  var_y <- "cost"
  
  form <- as.formula(str_c(var_y, "~", str_c(var_x, collapse = " + ")))
  
  # for one-part, only use gau_id, because still lots of 0, converge issues happens
  mod_lm <- glm(data = d4s,
                formula = form, 
                family = gaussian("identity"))
  
  lm_coef <- coef(mod_lm) # Keep the coefficient from linear regression for later use
  # # > Poisson - Identity (LM start) ----
  # # this is because the coef. of fu_prop is extremely negative,
  # # as the start value, which can make estimated cost to below 0
  # # that will make some glm generate error message
  # # so change it to 0
  # lm_coef["fu_prop"] <- 10000
  # 
  # mod_gl <- glm(data = d4s,
  #               formula = form, 
  #               family = poisson("identity"),
  #               start = lm_coef
  # )
  # 
  mod_list <- list(
    one_part = list(mod_1p_gau_id = mod_lm)
  )
  
  # Two-part model ----
  
  # Model two parts 
  # > Part 1: probability of any costs incurring in the period
  # > Part 2: costs conditional on any costs incurring in the period
  
  ## Part 1 ----
  
  # Probability of any costs incurring in the period
  
  # > Logistic regression ----
  
  # Define the formula: outcome ~ covariate
  var_y <- "cost01"
  
  form <- as.formula(str_c(var_y, "~", str_c(var_x, collapse = " + ")))

  # Fit the logistic regression model to the data 
  mod_2p_p1 <- glm(data = d4s, 
                formula = form, 
                family = binomial(link = "logit")) 
  
  mod_list[["two_part"]][["p1"]] <- mod_2p_p1
  
  ## Part 2 ----
  
  # Costs conditional on any costs incurring in the period
  
  # Select the records with positive cost outcome for the part 2 model
  ana <- d4s %>% filter(cost > 0)
  
  # Define the formula: outcome ~ covariate
  var_y <- "cost"
  
  form <- as.formula(str_c(var_y, "~", str_c(var_x, collapse = " + ")))
  
  # > Gaussian - Identity ----
  name_test <- "gau_id"
  mod_p2_gau_id <- glm(data = ana, 
                   formula = form, 
                   family = list_test[[name_test]])
  mod_list[["two_part"]][["p2"]][["mod_p2_gau_id"]] <- mod_p2_gau_id
  
  lm_coef <- coef(mod_p2_gau_id) # Keep the coefficient from linear regression for later use
  # # this is because the coef. of fu_prop is extremely negative,
  # # as the start value, which can make estimated cost to below 0
  # # that will make some glm generate error message
  # # so change it to 0
  # lm_coef["fu_prop"] <- 0
  
  # > Gaussian - LOG ----
  if ("gau_log" %in% name_mod){
    name_test <- "gau_log"
    mod_p2_gau_log <- glm(data = ana, 
                          formula = form, 
                          family = list_test[[name_test]])
    mod_list[["two_part"]][["p2"]][["mod_p2_gau_log"]] <- mod_p2_gau_log
  }
  
  # > Poisson - Identity ----
  if ("poi_id" %in% name_mod){
    name_test <- "poi_id"
    mod_p2_poi_id <- glm(data = ana %>% mutate(cost = round(cost, 0)),
                         # Poisson regression usually require outcome to be rounded to 1
                         formula = form,
                         family = list_test[[name_test]],
                         start = lm_coef
                         # Sometimes poisson - identity GLM will need a good initial coefficients
                         # as identity link is not the default link for poisson regression
    )
    mod_list[["two_part"]][["p2"]][["mod_p2_poi_id"]] <- mod_p2_poi_id
  }
  
  # > Poisson - Log ----
  if ("poi_log" %in% name_mod){
    name_test <- "poi_log"
    mod_p2_poi_log <- glm(data = ana %>% mutate(cost = round(cost, 0)), 
                          # Poisson regression usually require outcome to be rounded to 1
                          formula = form, 
                          family = list_test[[name_test]])
    mod_list[["two_part"]][["p2"]][["mod_p2_poi_log"]] <- mod_p2_poi_log
  }
  
  # > Gamma - Identity (LM start) ----
  if ("gam_id" %in% name_mod){
    name_test <- "gam_id"
    mod_p2_gam_id <- glm(data = ana,
                         formula = form,
                         family = list_test[[name_test]],
                         start = lm_coef
                         # Sometimes Gamma - identity GLM will need a good initial coefficients
                         # as identity link is not the default link for Gamma regression
    )
    mod_list[["two_part"]][["p2"]][["mod_p2_gam_id"]] <- mod_p2_gam_id
  }
  
  # > Gamma - LOG ----
  if ("gam_log" %in% name_mod){
    name_test <- "gam_log"
    mod_p2_gam_log <- glm(data = ana, 
                          formula = form, 
                          family = list_test[[name_test]])
    mod_list[["two_part"]][["p2"]][["mod_p2_gam_log"]] <- mod_p2_gam_log
  }
  
  # > Gamma - squared root ----
  if ("gam_sqrt" %in% name_mod){
    name_test <- "gam_sqrt"
    mod_p2_gam_sqrt <- glm(data = ana,
                           formula = form,
                           family = list_test[[name_test]],
                           start = lm_coef)
    mod_list[["two_part"]][["p2"]][["mod_p2_gam_sqrt"]] <- mod_p2_gam_sqrt
  }
  
  return(mod_list)
}


# Do tests ----------------------------------------------------------------

do_tests <- function(mod_list){
  
  ## within Two-part Test----
  # name_mod <- c("gau_id", "gau_log", "poi_log",  "gam_log")
  name_mod <- c("gau_id", "gau_log", "poi_id", "poi_log", "gam_id", "gam_log", "gam_sqrt")
  n <- length(name_mod)
  tmp <- rep(list(NA), n)
  
  # Select the records with positive cost outcome for the part 2 model
  ana <- d4s %>% filter(cost > 0)
  
  for(i in 1:n){
    mod <- mod_list[["two_part"]][["p2"]][[paste0("mod_p2_",name_mod[[i]])]] 
    ana2 <- with(mod, tibble(epatid = ana$epatid,
                             y = ana$cost,
                             y_hat = fitted.values,
                             y_link = linear.predictors)) %>% 
      mutate(res = y - y_hat)
    family <- mod$family
    test_hl <- with(ana2, f_test_spe_hl_clx(epatid, y_hat, res))
    test_pl <- with(ana2, f_test_spe_pl_clx(family = family, epatid, y, y_link))
    test_mp <- with(ana2, f_test_spe_mp_clx(epatid, y_hat, res))
    test_gof <- f_test_gof(ana2$res)
    tmp[[i]] <- bind_cols(test_hl, test_pl, test_mp, test_gof)
  }
  
  names(tmp) <- name_mod
  tp1 <- bind_rows(tmp, .id = "mod")
  rt1 <- tp1 %>% 
    mutate_at(c("hl_p", "pl_p", "mp_slope"), 
              ~ifelse(. < 0.01, "<0.01", as.character(round(., 2)))) %>% 
    select(mod,mp_slope, hl_p, pl_p, me, mae, rmse)
  
  
  ## Between one-part and two-part models----
  
  # One-part model
  name_mod <- c("gau_id")
  tmp1 <- list()
  # for one-part, only use gau_id, because still lots of 0, converge issues happens
  mod_1p <- mod_list[["one_part"]][[paste0("mod_1p_",name_mod[[1]])]] 
  
  ana <- with(mod_1p, tibble(y = data$cost,
                             y_hat = fitted.values)) %>%
    mutate(res = y_hat - y)
  tmp1[[1]] <- f_test_gof(ana$res)
  names(tmp1) <- str_c("1p_", name_mod)
  
  
  # Two-part model
  # this time: poi_log 
  mod_2p_p1 <- mod_list[["two_part"]][["p1"]] 
  
  ana2_p1 <- d4s %>% 
    select(epatid, admi_year, cost) %>% 
    bind_cols(cost_p1 = mod_2p_p1$fitted.values)
  
  # name_mod <- c("gau_id", "gau_log", "poi_log", "gam_log")
  name_mod <- c("gau_id", "gau_log", "poi_id", "poi_log", "gam_id", "gam_log", "gam_sqrt")
  n <- length(name_mod)
  tmp2 <- rep(list(NA), n)
  for(i in 1:n){
    mod_2p_p2 <- mod_list[["two_part"]][["p2"]][[paste0("mod_p2_",name_mod[[i]])]]  
    ana2_p2 <- predict(mod_2p_p2, newdata= d4s, type = "response")
    ana2 <- bind_cols(ana2_p1, tibble(cost_p2 = ana2_p2)) %>% 
      mutate(y_hat = cost_p1 * cost_p2,
             res = y_hat - cost)
    tmp2[[i]] <- f_test_gof(ana2$res)
  }
  
  names(tmp2) <- str_c("2p_p2_", name_mod)
  rt2 <- bind_rows(c(tmp1, tmp2), .id = "mod")

  rt_list <- list(rt1, rt2)
  
  return(rt_list)
}


# Wrap model and test -----------------------------------------------------
# name_mod <- c("gau_id", "gau_log", "poi_id", "poi_log", "gam_id", "gam_log", "gam_sqrt")
model_test <- function(cancer_type, he_ova, var_x, name_mod){
  
  mod_list <- model_1p_2p(cancer_type, he_ova, var_x)
  
  rt <- do_tests(mod_list)
  
  return(rt)
  
}






















