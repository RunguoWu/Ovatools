##################################################################
##                 R code for the decision tree                 ##
##################################################################

# we allow the decision tree parameters vary across individuals
# parameters depend on individual characteristics
# e.g. age, ethnicity and Townsend score
# we can pre-estimate all parameters for each individual before entering decision tree

library(tidyverse)
# library(bannerCommenter)
library(dtplyr)
library(data.table)

# Derive decision tree parameters -----------------------------------------
# prevalence: probabilities of ovarian cancer (prevalence) using individual characteristics
# positive test results/PTR: probability of test positive, e.g. CA125>=35 or US+
# PPV: positive predictive value
# NPV: negative predictive value

# PTR = (sensitivity * prevalence) + (1-specificity)*(1-prevalence)
# PPV = (sensitivity * prevalence)/PTR
# NPV = (specificity * (1-prevalence))/(1-PTR)

# formula to calculate true positive, false positive, false negative and True negative
# TP = PTR * PPV
# FP = PTR * (1 - PPV)
# FN = (1 - PTR) * (1 - NPV)
# TN = (1 - PTR) * NPV

# alternative formula
# TP = prevalence * sensitivity
# FP = (1 - prevalence) * (1 - specificity)
# FN = prevalence * (1 - sensitivity)
# TN = (1 - prevalence) * specificity

# Decision tree 1 ---------------------------------------------------------
# patients receive CA125 first, followed by Ultrasound if CA125 positive
# each test followed by a telephone consultation, i.e. cost_gp2
# before referral there is one full consultation, i.e. cost_gp1
# ca125 test along with nurse cost for drawing blood

## Decision tree 1.0----
# 35 U/ml threshold 

decision_tree100 <- function(dt=dt, accuracy=accuracy, rho=0, 
                             cost_us_cdc=FALSE, cost_us_op=FALSE, cost_us_input=NA,
                             us_sens_adj=0, us_spec_adj=0
                             ){

  # specificity is always the same for early and late
  sens_ca125_early <- accuracy$accuracy_early["sens_ca125"]
  sens_ca125_late <- accuracy$accuracy_late["sens_ca125"]
  spec_ca125 <- accuracy$accuracy_early["spec_ca125"] # the same as late

  sens_us_early <- accuracy$accuracy_early["sens_us"] + us_sens_adj
  sens_us_late <- accuracy$accuracy_late["sens_us"] + us_sens_adj
  spec_us <- accuracy$accuracy_early["spec_us"] + us_spec_adj
  
  sens_us_ca125_early <- accuracy$accuracy_early["sens_us_fl_ca125"] + us_sens_adj
  sens_us_ca125_late <- accuracy$accuracy_late["sens_us_fl_ca125"] + us_sens_adj
  spec_us_ca125 <- accuracy$accuracy_early["spec_us_fl_ca125"] + us_spec_adj
  
  # either from accuracy_early or late, the same
  if(cost_us_cdc) cost_us=accuracy$accuracy_early["cost_us_cdc"] else
    if(cost_us_op) cost_us=accuracy$accuracy_early["cost_us_op"] else
      cost_us=accuracy$accuracy_early["cost_us"]
  
  if(!is.na(cost_us_input)) cost_us = cost_us_input
  
  cost_ca125=accuracy$accuracy_early["cost_ca125"]
  cost_gp1=accuracy$accuracy_early["cost_gp1"]
  cost_gp2=accuracy$accuracy_early["cost_gp2"]
  cost_nurse=accuracy$accuracy_early["cost_nurse"]
  
  fp_uter_us = accuracy$accuracy_early["fp_uter_us"]
  fp_loGI_us = accuracy$accuracy_early["fp_loGI_us"]
  fp_lung_us = accuracy$accuracy_early["fp_lung_us"]
  fp_panc_us = accuracy$accuracy_early["fp_panc_us"]

  # expected probabilities for each pathway
  ep <- dt %>%
    lazy_dt() %>%
    mutate(prev12 = prevalence * (1-ova_stage34),
           prev34 = prevalence * ova_stage34) %>% 
    mutate(
      tp_ca125_12 = prev12 * sens_ca125_early,
      tp_ca125_34 = prev34 * sens_ca125_late, 
      fp_ca125 = (1 - prev12 - prev34) * (1 - spec_ca125),
      ptr_ca125 = tp_ca125_12 + tp_ca125_34 + fp_ca125,

      e_TP12 = tp_ca125_12 * sens_us_ca125_early,
      e_TP34 = tp_ca125_34 * sens_us_ca125_late,
      
      e_FP = fp_ca125 * (1 - spec_us_ca125), 
      
      ptr_us_fl_ca125 = e_TP12 + e_TP34 + e_FP,
      
      e_TN = (1 - prev12 - prev34) * spec_ca125 + fp_ca125 * spec_us_ca125, 
      
      e_FN12 = prev12 * (1 - sens_ca125_early) + tp_ca125_12 * (1 - sens_us_ca125_early),
      e_FN34 = prev34 * (1 - sens_ca125_late) + tp_ca125_34 * (1 - sens_us_ca125_late),
      e_FN = e_FN12 + e_FN34,
      # check = e_TP12+e_TP34+e_FP+e_TN+e_FN12+e_FN34,
  
      # other cancer 1: uterus, lower GI
      # e_otc1_FP = pmin(prevalence_otc1, e_FP * accuracy["otc1_us_fp"]),
      # e_otc1_TN = prevalence_otc1 - e_otc1_FP,
      e_uter_FP = pmin(pred_uter, e_FP * fp_uter_us),
      e_uter_TN = pred_uter - e_uter_FP,
      e_loGI_FP = pmin(pred_loGI, e_FP * fp_loGI_us),
      e_loGI_TN = pred_loGI - e_loGI_FP,
      # other cancer 2: lung, pancreas
      # e_otc2_FP = pmin(prevalence_otc2, e_FP * accuracy["otc2_us_fp"]),
      # e_otc2_TN = prevalence_otc2 - e_otc2_FP,
      e_lung_FP = pmin(pred_lung, e_FP * fp_lung_us),
      e_lung_TN = pred_lung - e_lung_FP,
      e_panc_FP = pmin(pred_panc, e_FP * fp_panc_us),
      e_panc_TN = pred_panc - e_panc_FP,
      
      # record those referred without US+ or CA125+, in this pathway all 0
      e_TP_woUS = 0,
      e_FP_woUS = 0,
      prev_us_woUS_FP = 0,
      e_TP_woCA125= 0,
      e_FP_woCA125= 0,
      prev_ca125_woCA125_FP= 0,
      
      # cost
      cost = cost_gp1 + cost_nurse + cost_ca125 +  
        ptr_ca125*(cost_gp2 + cost_us) + ptr_us_fl_ca125 * cost_gp1,
      # assume those with cancer (ovarian and others) not detect initially were finally diganosed but
      # repeated the whole process so the cost was doubled
      cost = 2 * cost * (e_FN + e_uter_TN + e_loGI_TN + e_lung_TN + e_panc_TN) + 
        cost * (1 - e_FN - e_uter_TN - e_loGI_TN - e_lung_TN - e_panc_TN)
      
    ) %>% 
    mutate(sum = rowSums(across(c(e_TN, e_FN12, e_FN34, e_TP12, e_TP34, e_FP)))) %>% 
    select(id, e_TN, e_FN12, e_FN34, e_TP12, e_TP34, e_FP, sum, 
           e_TP_woUS, e_FP_woUS, prev_us_woUS_FP,
           e_TP_woCA125, e_FP_woCA125, prev_ca125_woCA125_FP,
           e_uter_FP, e_uter_TN, e_loGI_FP, e_loGI_TN, 
           e_lung_FP, e_lung_TN, e_panc_FP, e_panc_TN, cost) %>% 
    as.data.frame()
  
  return(ep)
}

## Decision tree 1.1----
# Ovatool with threshold
# default 1% and 3%

decision_tree101 <- function(dt=dt, accuracy=accuracy, rho=0,  
                             cost_us_op=FALSE, cost_us_cdc=FALSE,cost_us_input=NA,
                             us_sens_adj=0, us_spec_adj=0){

  # specificity is always the same for early and late
  sens_ovt1_early <- accuracy$accuracy_early["sens_ovt1"]
  sens_ovt1_late <- accuracy$accuracy_late["sens_ovt1"]
  spec_ovt1 <- accuracy$accuracy_early["spec_ovt1"] 
  
  sens_ovt3_early <- accuracy$accuracy_early["sens_ovt3"]
  sens_ovt3_late <- accuracy$accuracy_late["sens_ovt3"]
  spec_ovt3 <- accuracy$accuracy_early["spec_ovt3"] 
  
  sens_us_early <- accuracy$accuracy_early["sens_us"] + us_sens_adj
  sens_us_late <- accuracy$accuracy_late["sens_us"] + us_sens_adj
  spec_us <- accuracy$accuracy_early["spec_us"] + us_spec_adj
  
  sens_us_ca125_early <- accuracy$accuracy_early["sens_us_fl_ca125"] + us_sens_adj
  sens_us_ca125_late <- accuracy$accuracy_late["sens_us_fl_ca125"] + us_sens_adj
  spec_us_ca125 <- accuracy$accuracy_early["spec_us_fl_ca125"] + us_spec_adj
  
  # either from accuracy_early or late, the same
  if(cost_us_cdc) cost_us=accuracy$accuracy_early["cost_us_cdc"] else
    if(cost_us_op) cost_us=accuracy$accuracy_early["cost_us_op"] else
      cost_us=accuracy$accuracy_early["cost_us"]
  
  if(!is.na(cost_us_input)) cost_us = cost_us_input
  
  cost_ca125=accuracy$accuracy_early["cost_ca125"]
  cost_gp1=accuracy$accuracy_early["cost_gp1"]
  cost_gp2=accuracy$accuracy_early["cost_gp2"]
  cost_nurse=accuracy$accuracy_early["cost_nurse"]
  
  fp_uter_us = accuracy$accuracy_early["fp_uter_us"]
  fp_loGI_us = accuracy$accuracy_early["fp_loGI_us"]
  fp_lung_us = accuracy$accuracy_early["fp_lung_us"]
  fp_panc_us = accuracy$accuracy_early["fp_panc_us"]
  fp_uter_ovt3 = accuracy$accuracy_early["fp_uter_ovt3"]
  fp_loGI_ovt3 = accuracy$accuracy_early["fp_loGI_ovt3"]
  fp_lung_ovt3 = accuracy$accuracy_early["fp_lung_ovt3"]
  fp_panc_ovt3 = accuracy$accuracy_early["fp_panc_ovt3"]
  
  ep <- dt %>%
    lazy_dt() %>%
    mutate(prev12 = prevalence * (1-ova_stage34),
           prev34 = prevalence * ova_stage34) %>% 
    mutate(
      
      tp_ovt1_12 = prev12 * sens_ovt1_early,
      tp_ovt1_34 = prev34 * sens_ovt1_late,
      fp_ovt1 = (1 - prev12 - prev34) * (1 - spec_ovt1),
      ptr_ovt1 = tp_ovt1_12 + tp_ovt1_34 + fp_ovt1,
      
      tp_ovt3_12 = prev12 * sens_ovt3_early,
      tp_ovt3_34 = prev34 * sens_ovt3_late,
      fp_ovt3 = (1 - prev12 - prev34) * (1 - spec_ovt3),
      ptr_ovt3 = tp_ovt3_12 + tp_ovt3_34 + fp_ovt3,
      
      e_TP12 = (tp_ovt1_12 - tp_ovt3_12) * sens_us_ca125_early + tp_ovt3_12,
      e_TP34 = (tp_ovt1_34 - tp_ovt3_34) * sens_us_ca125_late + tp_ovt3_34,
      
      e_FP1 = (fp_ovt1 - fp_ovt3) * (1 - spec_us_ca125),
      e_FP2 = fp_ovt3,
      e_FP = e_FP1 + e_FP2,
      
      e_TN = (1 - prev12 - prev34) * spec_ovt1 + (fp_ovt1 - fp_ovt3) * spec_us_ca125, 
      
      e_FN12 = prev12 * (1 - sens_ovt1_early) + (tp_ovt1_12 - tp_ovt3_12) * (1 - sens_us_ca125_early),
      e_FN34 = prev34 * (1 - sens_ovt1_late) + (tp_ovt1_34 - tp_ovt3_34) * (1 - sens_us_ca125_late),
      e_FN = e_FN12 + e_FN34,
      
      ptr_us_fl_ovt = e_TP12 + e_TP34 + e_FP,
      
      check = e_TP12+e_TP34+e_FP+e_TN+e_FN12+e_FN34,
 
      # record those referred without US+
      e_TP_woUS = (tp_ovt3_12 + tp_ovt3_34),
      e_FP_woUS = fp_ovt3,
      # FP in the hospital US check for FP from OVatools 3% of primary care pathway
      prev_us_woUS_FP = (e_TP_woUS+e_FP_woUS) * (1 - e_TP_woUS/(e_TP_woUS+e_FP_woUS)) * (1 - spec_us_ca125),

      # no one referred without ca125+ (or Ovatools 1/3%+)
      e_TP_woCA125= 0,
      e_FP_woCA125= 0,
      prev_ca125_woCA125_FP= 0,
      
      e_uter_FP1 = e_FP1 * fp_uter_us,
      e_uter_FP2 = e_FP2 * fp_uter_ovt3,
      e_uter_FP = pmin(pred_uter, e_uter_FP1 + e_uter_FP2),
      e_uter_TN = pred_uter - e_uter_FP,
      
      e_loGI_FP1 = e_FP1 * fp_loGI_us,
      e_loGI_FP2 = e_FP2 * fp_loGI_ovt3,
      e_loGI_FP = pmin(pred_loGI, e_loGI_FP1 + e_loGI_FP2),
      e_loGI_TN = pred_loGI - e_loGI_FP,

      e_lung_FP1 = e_FP1 * fp_lung_us,
      e_lung_FP2 = e_FP2 * fp_lung_ovt3,
      e_lung_FP = pmin(pred_lung, e_lung_FP1 + e_lung_FP2),
      e_lung_TN = pred_lung - e_lung_FP,
      
      e_panc_FP1 = e_FP1 * fp_panc_us,
      e_panc_FP2 = e_FP2 * fp_panc_ovt3,
      e_panc_FP = pmin(pred_panc, e_panc_FP1 + e_panc_FP2),
      e_panc_TN = pred_panc - e_panc_FP,
      
      # cost
      cost=cost_gp1 + cost_nurse + cost_ca125 +
        (ptr_ovt1 - ptr_ovt3)*(cost_us + cost_gp2) + ptr_us_fl_ovt * cost_gp1,
      # assume those with cancer (ovarian and others) not detect initially were finally diganosed but
      # repeated the whole process so the cost was doubled
      cost = 2 * cost * (e_FN + e_uter_TN + e_loGI_TN + e_lung_TN + e_panc_TN) + 
        cost * (1 - e_FN - e_uter_TN - e_loGI_TN - e_lung_TN - e_panc_TN)
      
    ) %>% 
    mutate(sum = rowSums(across(c(e_TN, e_FN12, e_FN34, e_TP12, e_TP34, e_FP)))) %>% 
    select(id, e_TN, e_FN12, e_FN34, e_TP12, e_TP34, e_FP, sum, 
           e_TP_woUS, e_FP_woUS, prev_us_woUS_FP,
           e_TP_woCA125, e_FP_woCA125, prev_ca125_woCA125_FP,
           e_uter_FP, e_uter_TN, e_loGI_FP, e_loGI_TN, 
           e_lung_FP, e_lung_TN, e_panc_FP, e_panc_TN, cost) %>% 
    as.data.frame()
  
  return(ep)
}

# Decision tree 2 ---------------------------------------------------------
# patients receive US first and in Ovatools followed by Ovatools if US positive

## Decision tree 2.0----
# US only  

# Decision tree 3 ---------------------------------------------------------
# patients receive US and CA125 at the same time
## Decision tree 3.0----
# US+Ovatools, referral = US+ or CA125+ 
decision_tree110 <- function(dt=dt, accuracy=accuracy, rho = 0,
                             cost_us_cdc=FALSE, cost_us_op=FALSE, cost_us_input=NA,
                             us_sens_adj=0, us_spec_adj=0){
  
  # specificity is always the same for early and late
  sens_ca125_early <- accuracy$accuracy_early["sens_ca125"]
  sens_ca125_late <- accuracy$accuracy_late["sens_ca125"]
  spec_ca125 <- accuracy$accuracy_early["spec_ca125"] # the same as late
  
  sens_us_early <- accuracy$accuracy_early["sens_us"] + us_sens_adj
  sens_us_late <- accuracy$accuracy_late["sens_us"] + us_sens_adj
  spec_us <- accuracy$accuracy_early["spec_us"] + us_spec_adj
  
  sens_us_ca125_early <- accuracy$accuracy_early["sens_us_fl_ca125"] + us_sens_adj
  sens_us_ca125_late <- accuracy$accuracy_late["sens_us_fl_ca125"] + us_sens_adj
  spec_us_ca125 <- accuracy$accuracy_early["spec_us_fl_ca125"] + us_spec_adj
  
  # either from accuracy_early or late, the same
  if(cost_us_cdc) cost_us=accuracy$accuracy_early["cost_us_cdc"] else
    if(cost_us_op) cost_us=accuracy$accuracy_early["cost_us_op"] else
      cost_us=accuracy$accuracy_early["cost_us"]
  
  if(!is.na(cost_us_input)) cost_us = cost_us_input
  
  cost_ca125=accuracy$accuracy_early["cost_ca125"]
  cost_gp1=accuracy$accuracy_early["cost_gp1"]
  cost_gp2=accuracy$accuracy_early["cost_gp2"]
  cost_nurse=accuracy$accuracy_early["cost_nurse"]
  
  fp_uter_us = accuracy$accuracy_early["fp_uter_us"]
  fp_loGI_us = accuracy$accuracy_early["fp_loGI_us"]
  fp_lung_us = accuracy$accuracy_early["fp_lung_us"]
  fp_panc_us = accuracy$accuracy_early["fp_panc_us"]
  fp_uter_ovt3 = accuracy$accuracy_early["fp_uter_ovt3"]
  fp_loGI_ovt3 = accuracy$accuracy_early["fp_loGI_ovt3"]
  fp_lung_ovt3 = accuracy$accuracy_early["fp_lung_ovt3"]
  fp_panc_ovt3 = accuracy$accuracy_early["fp_panc_ovt3"]
  
  # NOTE: the big assumption is that CA125 and US results are independent
  ep <- dt %>%
    lazy_dt() %>%
    mutate(prev12 = prevalence * (1-ova_stage34),
           prev34 = prevalence * ova_stage34) %>% 
    mutate(
      # US- & Ovatools<3%
      # calculate joint accuarcy under a certain rho - correlation between the two tests
      # the basica assumption is rho = 0
      joint_sens_adj_early = rho * sqrt((1 - sens_ca125_early) * sens_ca125_early * (1 - sens_us_ca125_early) * sens_us_ca125_early),
      joint_sens_adj_late = rho * sqrt((1 - sens_ca125_late) * sens_ca125_late * (1 - sens_us_ca125_late) * sens_us_ca125_late),
      joint_spec_adj = rho * sqrt((1 - spec_ca125) * spec_ca125 * (1 - spec_us_ca125) * spec_us_ca125),
      
      e_TN = (1 - prev12 - prev34) * (spec_ca125 * spec_us_ca125 + joint_spec_adj),
      
      e_FN12 = prev12 * ((1 - sens_ca125_early) * (1 - sens_us_ca125_early) + joint_sens_adj_early),
      e_FN34 = prev34 * ((1 - sens_ca125_late) * (1 - sens_us_ca125_late) + joint_sens_adj_late),
      e_FN = e_FN12 + e_FN34,
      
      # US+ or Ovatools>=3%
      e_TP12 = prev12 - e_FN12,
      e_TP34 = prev34 - e_FN34,
      e_TP = e_TP12 + e_TP34,
      
      e_FP = (1 - prev12 - prev34) - e_TN,
      e_FP2 = (1 - prev12 - prev34) * (1 - spec_ca125),
      e_FP1 = pmax(e_FP - e_FP2, 0),
      
      # for those referred without US+
      # consider US check in hospital before applying surgery in FPs
      e_TP_woUS = e_TP - prev12 * sens_us_ca125_early - prev34 * sens_us_ca125_late,
      e_FP_woUS = e_FP - (1 - prev12 - prev34) * (1 - spec_us_ca125),
      prev_us_woUS_FP = (e_TP_woUS+e_FP_woUS) * (1 - e_TP_woUS/(e_TP_woUS+e_FP_woUS)) * (1 - spec_us_ca125),
      
      # for those referred without CA125
      # consider CA125 check in hospital before apply surgery in FPs
      e_TP_woCA125 = e_TP - prev12*sens_ca125_early - prev34*sens_ca125_late,
      e_FP_woCA125 = e_FP1,
      prev_ca125_woCA125_FP = (e_TP_woCA125+e_FP_woCA125) * (1 - e_TP_woCA125/(e_TP_woCA125+e_FP_woCA125)) * (1 - spec_ca125),
      
      # other cancer 1: uterus, lower GI
      e_uter_FP1 = e_FP1 * fp_uter_us,
      e_uter_FP2 = e_FP2 * fp_uter_ovt3,
      e_uter_FP = pmin(pred_uter, e_uter_FP1 + e_uter_FP2),
      e_uter_TN = pred_uter - e_uter_FP,
      
      e_loGI_FP1 = e_FP1 * fp_loGI_us,
      e_loGI_FP2 = e_FP2 * fp_loGI_ovt3,
      e_loGI_FP = pmin(pred_loGI, e_loGI_FP1 + e_loGI_FP2),
      e_loGI_TN = pred_loGI - e_loGI_FP,
      # other cancer 2: lung, pancreas
      e_lung_FP1 = e_FP1 * fp_lung_us,
      e_lung_FP2 = e_FP2 * fp_lung_ovt3,
      e_lung_FP = pmin(pred_lung, e_lung_FP1 + e_lung_FP2),
      e_lung_TN = pred_lung - e_lung_FP,
      
      e_panc_FP1 = e_FP1 * fp_panc_us,
      e_panc_FP2 = e_FP2 * fp_panc_ovt3,
      e_panc_FP = pmin(pred_panc, e_panc_FP1 + e_panc_FP2),
      e_panc_TN = pred_panc - e_panc_FP,
      
      # cost
      cost = cost_gp1 + cost_us + cost_nurse + cost_ca125 + (e_TP+e_FP) * cost_gp1,
      # assume those with cancer (ovarian and others) not detect initially were finally diganosed but
      # repeated the whole process so the cost was doubled
      cost = 2 * cost * (e_FN + e_uter_TN + e_loGI_TN + e_lung_TN + e_panc_TN) + 
        cost * (1 - e_FN - e_uter_TN - e_loGI_TN - e_lung_TN - e_panc_TN)
    ) %>% 
    mutate(sum = rowSums(across(c(e_TN, e_FN12, e_FN34, e_TP12, e_TP34, e_FP)))) %>% 
    select(id, e_TN, e_FN12, e_FN34, e_TP12, e_TP34, e_FP, sum, 
           e_TP_woUS, e_FP_woUS, prev_us_woUS_FP,
           e_TP_woCA125, e_FP_woCA125, prev_ca125_woCA125_FP,
           e_uter_FP, e_uter_TN, e_loGI_FP, e_loGI_TN, 
           e_lung_FP, e_lung_TN, e_panc_FP, e_panc_TN, cost) %>% 
    as.data.frame()
  return(ep)
}


## Decision tree 3.2----
# US+Ovatools, referral = US+ or Ovatools>=3% 
decision_tree112 <- function(dt=dt, accuracy=accuracy, rho = 0, 
                             cost_us_op=FALSE, cost_us_cdc=FALSE, cost_us_input=NA,
                             us_sens_adj=0, us_spec_adj=0){
  
  # specificity is always the same for early and late
  spec_ca125 <- accuracy$accuracy_early["spec_ca125"] # the same as late
  
  sens_ovt1_early <- accuracy$accuracy_early["sens_ovt1"]
  sens_ovt1_late <- accuracy$accuracy_late["sens_ovt1"]
  spec_ovt1 <- accuracy$accuracy_early["spec_ovt1"] 
  
  sens_ovt3_early <- accuracy$accuracy_early["sens_ovt3"]
  sens_ovt3_late <- accuracy$accuracy_late["sens_ovt3"]
  spec_ovt3 <- accuracy$accuracy_early["spec_ovt3"] 
  
  sens_us_ca125_early <- accuracy$accuracy_early["sens_us_fl_ca125"] + us_sens_adj
  sens_us_ca125_late <- accuracy$accuracy_late["sens_us_fl_ca125"] + us_sens_adj
  spec_us_ca125 <- accuracy$accuracy_early["spec_us_fl_ca125"] + us_spec_adj
  
  # either from accuracy_early or late, the same
  if(cost_us_cdc) cost_us=accuracy$accuracy_early["cost_us_cdc"] else
    if(cost_us_op) cost_us=accuracy$accuracy_early["cost_us_op"] else
      cost_us=accuracy$accuracy_early["cost_us"]
  
  if(!is.na(cost_us_input)) cost_us = cost_us_input
  
  cost_ca125=accuracy$accuracy_early["cost_ca125"]
  cost_gp1=accuracy$accuracy_early["cost_gp1"]
  cost_gp2=accuracy$accuracy_early["cost_gp2"]
  cost_nurse=accuracy$accuracy_early["cost_nurse"]
  
  fp_uter_us = accuracy$accuracy_early["fp_uter_us"]
  fp_loGI_us = accuracy$accuracy_early["fp_loGI_us"]
  fp_lung_us = accuracy$accuracy_early["fp_lung_us"]
  fp_panc_us = accuracy$accuracy_early["fp_panc_us"]
  fp_uter_ovt3 = accuracy$accuracy_early["fp_uter_ovt3"]
  fp_loGI_ovt3 = accuracy$accuracy_early["fp_loGI_ovt3"]
  fp_lung_ovt3 = accuracy$accuracy_early["fp_lung_ovt3"]
  fp_panc_ovt3 = accuracy$accuracy_early["fp_panc_ovt3"]
  
  # NOTE: the big assumption is that CA125 and US results are independent
  ep <- dt %>%
    lazy_dt() %>%
    mutate(prev12 = prevalence * (1-ova_stage34),
           prev34 = prevalence * ova_stage34) %>% 
    mutate(
      # US- & Ovatools<3%
      # calculate joint accuarcy under a certain rho - correlation between the two tests
      # the basica assumption is rho = 0
      joint_sens_adj_early = rho * sqrt((1 - sens_ovt3_early) * sens_ovt3_early * (1 - sens_us_ca125_early) * sens_us_ca125_early),
      joint_sens_adj_late = rho * sqrt((1 - sens_ovt3_late) * sens_ovt3_late * (1 - sens_us_ca125_late) * sens_us_ca125_late),
      joint_spec_adj = rho * sqrt((1 - spec_ovt3) * spec_ovt3 * (1 - spec_us_ca125) * spec_us_ca125),
      
      e_TN = (1 - prev12 - prev34) * (spec_ovt3 * spec_us_ca125 + joint_spec_adj),
      
      e_FN12 = prev12 * ((1 - sens_ovt3_early) * (1 - sens_us_ca125_early) + joint_sens_adj_early),
      e_FN34 = prev34 * ((1 - sens_ovt3_late) * (1 - sens_us_ca125_late) + joint_sens_adj_late),
      e_FN = e_FN12 + e_FN34,
      
      # US+ or Ovatools>=3%
      e_TP12 = prev12 - e_FN12,
      e_TP34 = prev34 - e_FN34,
      e_TP = e_TP12 + e_TP34,
      
      e_FP = (1 - prev12 - prev34) - e_TN,
      e_FP2 = (1 - prev12 - prev34) * (1 - spec_ovt3),
      e_FP1 = pmax(e_FP - e_FP2, 0),
      
      # for those referred without US+
      # consider US check in hospital before applying surgery in FPs
      e_TP_woUS = e_TP - prev12 * sens_us_ca125_early - prev34 * sens_us_ca125_late,
      e_FP_woUS = e_FP - (1 - prev12 - prev34) * (1 - spec_us_ca125),
      prev_us_woUS_FP = (e_TP_woUS+e_FP_woUS) * (1 - e_TP_woUS/(e_TP_woUS+e_FP_woUS)) * (1 - spec_us_ca125),
      
      # for those referred without CA125
      # consider CA125 check in hospital before apply surgery in FPs
      e_TP_woCA125 = e_TP - prev12*sens_ovt3_early - prev34*sens_ovt3_late,
      e_FP_woCA125 = e_FP1,
      prev_ca125_woCA125_FP = (e_TP_woCA125+e_FP_woCA125) * (1 - e_TP_woCA125/(e_TP_woCA125+e_FP_woCA125)) * (1 - spec_ca125),
      
      # other cancer 1: uterus, lower GI
      e_uter_FP1 = e_FP1 * fp_uter_us,
      e_uter_FP2 = e_FP2 * fp_uter_ovt3,
      e_uter_FP = pmin(pred_uter, e_uter_FP1 + e_uter_FP2),
      e_uter_TN = pred_uter - e_uter_FP,
      
      e_loGI_FP1 = e_FP1 * fp_loGI_us,
      e_loGI_FP2 = e_FP2 * fp_loGI_ovt3,
      e_loGI_FP = pmin(pred_loGI, e_loGI_FP1 + e_loGI_FP2),
      e_loGI_TN = pred_loGI - e_loGI_FP,
      # other cancer 2: lung, pancreas
      e_lung_FP1 = e_FP1 * fp_lung_us,
      e_lung_FP2 = e_FP2 * fp_lung_ovt3,
      e_lung_FP = pmin(pred_lung, e_lung_FP1 + e_lung_FP2),
      e_lung_TN = pred_lung - e_lung_FP,
      
      e_panc_FP1 = e_FP1 * fp_panc_us,
      e_panc_FP2 = e_FP2 * fp_panc_ovt3,
      e_panc_FP = pmin(pred_panc, e_panc_FP1 + e_panc_FP2),
      e_panc_TN = pred_panc - e_panc_FP,
      
      # cost
      cost = cost_gp1 + cost_us + cost_nurse + cost_ca125 + (e_TP+e_FP) * cost_gp1,
      # assume those with cancer (ovarian and others) not detect initially were finally diganosed but
      # repeated the whole process so the cost was doubled
      cost = 2 * cost * (e_FN + e_uter_TN + e_loGI_TN + e_lung_TN + e_panc_TN) + 
        cost * (1 - e_FN - e_uter_TN - e_loGI_TN - e_lung_TN - e_panc_TN)
    ) %>% 
    mutate(sum = rowSums(across(c(e_TN, e_FN12, e_FN34, e_TP12, e_TP34, e_FP)))) %>% 
    select(id, e_TN, e_FN12, e_FN34, e_TP12, e_TP34, e_FP, sum, 
           e_TP_woUS, e_FP_woUS, prev_us_woUS_FP,
           e_TP_woCA125, e_FP_woCA125, prev_ca125_woCA125_FP,
           e_uter_FP, e_uter_TN, e_loGI_FP, e_loGI_TN, 
           e_lung_FP, e_lung_TN, e_panc_FP, e_panc_TN, cost) %>% 
    as.data.frame()
  return(ep)
}

# combined decision tree --------------------------------------------------

decision_tree_family <- list(decision_tree100=decision_tree100, 
                             decision_tree110=decision_tree110, 
                             decision_tree101=decision_tree101,
                             decision_tree112=decision_tree112
) 

# by stage accuracy
# multiply the percentage of stage (early/late) with the result
# this approach in effect splits patients without ovarian cancer into two groups
# , by the same percentage for stage of ovarian cancer patients
# , and apply different accuracy data to the two groups
# However, for those without ovarian cancer, only specificity matters.
# And specificities are basically the same for the two groups.
# therefore this approach is feasible.

decision_tree <- function(ca125_first=1, us_first=0, ovatools=0, dt=dt, 
                          accuracy_tot=accuracy_tot, rho=0, 
                          cost_us_cdc=FALSE, cost_us_op=FALSE, cost_us_input=NA,
                          us_sens_adj=0, us_spec_adj=0) {
  
  # 1 means yes, 0 means no
  tree_fn <- str_c("decision_tree", ca125_first, us_first, ovatools)

  ep <- decision_tree_family[[tree_fn]](dt, accuracy_tot, rho=rho, 
                                        cost_us_cdc=cost_us_cdc, cost_us_op=cost_us_op, cost_us_input=cost_us_input,
                                        us_sens_adj=us_sens_adj, us_spec_adj=us_spec_adj)

  ep <- ep %>% 
    left_join(dt[, c("id", "group", "ova_stage34", "uter_stage34", "loGI_stage34", 
                     "lung_stage34", "panc_stage34", "other_stage34",
                     "pred_uter", "pred_loGI", "pred_lung", "pred_panc", "pred_other"
    )], by=join_by(id)) 

  ep <- ep %>% 
    mutate(
    e_uter_FP12 = e_uter_FP * (1-uter_stage34), 
    e_uter_FP34 = e_uter_FP * uter_stage34,
    
    e_loGI_FP12 = e_loGI_FP * (1-loGI_stage34), 
    e_loGI_FP34 = e_loGI_FP * loGI_stage34,
    
    e_lung_FP12 = e_lung_FP * (1-lung_stage34),
    e_lung_FP34 = e_lung_FP * lung_stage34,
    
    e_panc_FP12 = e_panc_FP * (1-panc_stage34),
    e_panc_FP34 = e_panc_FP * panc_stage34,

    stage_early_uter = pred_uter * (1-uter_stage34),
    stage_late_uter = pred_uter * uter_stage34,
    stage_early_loGI = pred_loGI * (1-loGI_stage34),
    stage_late_loGI = pred_loGI * loGI_stage34,
    stage_early_lung = pred_lung * (1-lung_stage34),
    stage_late_lung = pred_lung * lung_stage34,
    stage_early_panc = pred_panc * (1-panc_stage34),
    stage_late_panc = pred_panc * panc_stage34,
    stage_early_other = pred_other * (1-other_stage34),
    stage_late_other = pred_other * other_stage34,
  ) %>% 
    select(id, group, e_TP12, e_TP34, e_FN12, e_FN34, e_TN, e_FP, e_uter_FP12,
           e_uter_FP34, e_uter_TN, e_loGI_FP12, e_loGI_FP34, e_loGI_TN,
           e_lung_FP12, e_lung_FP34, e_lung_TN, e_panc_FP12, e_panc_FP34,
           e_panc_TN, stage_early_uter, stage_late_uter, stage_early_loGI,
           stage_late_loGI, stage_early_lung, stage_late_lung, stage_early_panc,
           stage_late_panc, stage_early_other, stage_late_other, cost, e_TP_woUS,
           e_FP_woUS, prev_us_woUS_FP, e_TP_woCA125, e_FP_woCA125, prev_ca125_woCA125_FP
           )
    
  
  return(ep)
}

























