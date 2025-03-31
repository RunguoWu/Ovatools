# adapte from the same name file in post_analysis folder

# use real distribution of pathways, combine the 3 current pathways to 1
with_real_dist <- function(decitree100, decitree010, decitree110,# three current pathways  
                                p_100, p_010, p_110 # the corresponding proportions
                                ){
  
  deci_list <- c("100", "010", "110")
  for (deci in deci_list) {
    
    x <- get(paste0("decitree", deci)) %>% 
      rename(
        !!sym(paste0("TP12_", deci)) := e_TP12, 
        !!sym(paste0("FN12_", deci)) := e_FN12,
        !!sym(paste0("TP34_", deci)) := e_TP34,
        !!sym(paste0("FN34_", deci)) := e_FN34, 
        !!sym(paste0("FP_", deci)) := e_FP,
        !!sym(paste0("TP_woUS_", deci)) := e_TP_woUS,
        !!sym(paste0("FP_woUS_", deci)) := e_FP_woUS,
        # !!sym(paste0("ptr_us_woUS_", deci)) := ptr_us_woUS,
        !!sym(paste0("prev_us_woUS_FP_", deci)) := prev_us_woUS_FP,
        
        !!sym(paste0("TP_woCA125_", deci)) := e_TP_woCA125,
        !!sym(paste0("FP_woCA125_", deci)) := e_FP_woCA125,
        !!sym(paste0("prev_ca125_woCA125_FP_", deci)) := prev_ca125_woCA125_FP,
        
        !!sym(paste0("uter_FP12_", deci)) := e_uter_FP12, 
        !!sym(paste0("loGI_FP12_", deci)) := e_loGI_FP12, 
        !!sym(paste0("lung_FP12_", deci)) := e_lung_FP12,
        !!sym(paste0("panc_FP12_", deci)) := e_panc_FP12, 
        !!sym(paste0("uter_FP34_", deci)) := e_uter_FP34,
        !!sym(paste0("loGI_FP34_", deci)) := e_loGI_FP34, 
        !!sym(paste0("lung_FP34_", deci)) := e_lung_FP34,
        !!sym(paste0("panc_FP34_", deci)) := e_panc_FP34,
        !!sym(paste0("cost_pc_", deci)) := cost
      ) %>% 
      mutate(
        !!sym(paste0("FP_true_", deci)) := !!sym(paste0("FP_", deci)) - 
          !!sym(paste0("uter_FP12_", deci)) - 
          !!sym(paste0("loGI_FP12_", deci)) - 
          !!sym(paste0("lung_FP12_", deci)) -
          !!sym(paste0("panc_FP12_", deci)) - 
          !!sym(paste0("uter_FP34_", deci)) -
          !!sym(paste0("loGI_FP34_", deci)) - 
          !!sym(paste0("lung_FP34_", deci)) -
          !!sym(paste0("panc_FP34_", deci))
      ) %>% 
      select(-e_TN, -e_uter_TN, -e_lung_TN, -e_loGI_TN,-e_panc_TN,
             -stage_early_uter, -stage_early_loGI, -stage_early_lung, -stage_early_panc,
             -stage_late_uter, -stage_late_loGI, -stage_late_lung, -stage_late_panc, 
             -stage_early_other, -stage_late_other) %>% 
      mutate(across(!!sym(paste0("TP12_", deci)) : !!sym(paste0("FP_true_", deci)), ~ .x*get(paste0("p_", deci))))
    
    assign(paste0("dt", deci), x)
  }
  
  dt0 <- dt100 %>% left_join(dt010) %>% left_join(dt110) %>% 
    mutate(
      TP12_0 = TP12_100+TP12_010+TP12_110, 
      FN12_0 = FN12_100+FN12_010+FN12_110,
      TP34_0 = TP34_100+TP34_010+TP34_110, 
      FN34_0 = FN34_100+FN34_010+FN34_110,
      FP_0 = FP_100+FP_010+FP_110,
      
      TP_woUS_0 = TP_woUS_100 + TP_woUS_010 + TP_woUS_110,
      FP_woUS_0 = FP_woUS_100 + FP_woUS_010 + FP_woUS_110,
      # ptr_us_woUS_0 = ptr_us_woUS_100 + ptr_us_woUS_010 +ptr_us_woUS_110,
      prev_us_woUS_FP_0 = prev_us_woUS_FP_100 + prev_us_woUS_FP_010 +prev_us_woUS_FP_110, 
      
      TP_woCA125_0 = TP_woCA125_100 + TP_woCA125_010 + TP_woCA125_110,
      FP_woCA125_0 = FP_woCA125_100 + FP_woCA125_010 + FP_woCA125_110,
      prev_ca125_woCA125_FP_0 = prev_ca125_woCA125_FP_100 + prev_ca125_woCA125_FP_010 +prev_ca125_woCA125_FP_110, 
      
      uter_FP12_0 = uter_FP12_100+uter_FP12_010+uter_FP12_110,
      loGI_FP12_0 = loGI_FP12_100+loGI_FP12_010+loGI_FP12_110,
      lung_FP12_0 = lung_FP12_100+lung_FP12_010+lung_FP12_110,
      panc_FP12_0 = panc_FP12_100+panc_FP12_010+panc_FP12_110,
      uter_FP34_0 = uter_FP34_100+uter_FP34_010+uter_FP34_110,
      loGI_FP34_0 = loGI_FP34_100+loGI_FP34_010+loGI_FP34_110,
      lung_FP34_0 = lung_FP34_100+lung_FP34_010+lung_FP34_110,
      panc_FP34_0 = panc_FP34_100+panc_FP34_010+panc_FP34_110,
      cost_pc_0 = cost_pc_100+cost_pc_010+cost_pc_110,
      FP_true_0 = FP_true_100+FP_true_010+FP_true_110
    ) %>% 
    select(id, TP12_0:FP_true_0)

  return(dt0)
}


# calculate the effect on stage shift
# decitree0 and decitree1 are results from decision trees, one for the old, one for the new.
# stage_shift is the file for the effect of early detection on stage
# with_real_dist: whether decitree0 is the combined result from real distribution

stage_shift_effect <- function(decitree0, decitree1, with_real_dist = FALSE, stage_shift,
                               shift_other = TRUE, 
                               shift_late_sens = "", shift_early_sens=""){
  
  shift_effects <- stage_shift
  
  if (shift_late_sens == "upper") shift_effects$ova_rr_late <- stage_shift$ova_rr_late_up
  if (shift_late_sens == "lower") shift_effects$ova_rr_late <- stage_shift$ova_rr_late_lo
  if (shift_early_sens == "upper") shift_effects$ova_rr_early <- stage_shift$ova_rr_early_up
  if (shift_early_sens == "lower") shift_effects$ova_rr_early <- stage_shift$ova_rr_early_lo
  
  if (!shift_other){
    shift_effects[["lung_rr_late"]] <- shift_effects[["lung_rr_early"]] <- 1
    shift_effects[["loGI_rr_late"]] <- shift_effects[["loGI_rr_early"]] <- 1
    shift_effects[["panc_rr_late"]] <- shift_effects[["panc_rr_early"]] <- 1
    shift_effects[["uter_rr_late"]] <- shift_effects[["uter_rr_early"]] <- 1
  }
  
  
  # import decision tree results
  
  if(with_real_dist == FALSE){
    dt0 <- decitree0 %>% 
      # select(id, e_TP12, e_FN12, e_TP34, e_FN34, e_FP, stage_early, stage_late, no_ova, cost) %>% 
      rename(
        TP12_0 = e_TP12, 
        FN12_0 = e_FN12,
        TP34_0 = e_TP34, 
        FN34_0 = e_FN34,
        FP_0 = e_FP,
        TP_woUS_0 = e_TP_woUS,
        FP_woUS_0 = e_FP_woUS,
        # ptr_us_woUS_0 = ptr_us_woUS,
        prev_us_woUS_FP_0 = prev_us_woUS_FP, 
        
        TP_woCA125_0 = e_TP_woCA125,
        FP_woCA125_0 = e_FP_woCA125,
        prev_ca125_woCA125_FP_0 = prev_ca125_woCA125_FP, 
        
        uter_FP12_0 = e_uter_FP12,
        loGI_FP12_0 = e_loGI_FP12,
        lung_FP12_0 = e_lung_FP12,
        panc_FP12_0 = e_panc_FP12,
        uter_FP34_0 = e_uter_FP34,
        loGI_FP34_0 = e_loGI_FP34,
        lung_FP34_0 = e_lung_FP34,
        panc_FP34_0 = e_panc_FP34,
        cost_pc_0 = cost
      ) %>% 
      mutate(FP_true_0 = FP_0 - uter_FP12_0 - loGI_FP12_0 - lung_FP12_0 - panc_FP12_0
             - uter_FP34_0 - loGI_FP34_0 - lung_FP34_0 - panc_FP34_0) %>% 
      select(-e_TN, -e_uter_TN, -e_lung_TN, -e_loGI_TN,-e_panc_TN,
             -stage_early_uter, -stage_early_loGI, -stage_early_lung, -stage_early_panc,
             -stage_late_uter, -stage_late_loGI, -stage_late_lung, -stage_late_panc, 
             -stage_early_other, -stage_late_other, -epatid)
  } else dt0 <- decitree0
  
  dt1 <- decitree1 %>% 
    mutate(
      stage_early = e_TP12 + e_FN12,
      stage_late = e_TP34 + e_FN34,
      no_ova = e_TN + e_FP
    ) %>% 
    rename(
      TP12_1 = e_TP12, 
      FN12_1 = e_FN12,
      TP34_1 = e_TP34, 
      FN34_1 = e_FN34, 
      FP_1 = e_FP,
      TP_woUS_1 = e_TP_woUS,
      FP_woUS_1 = e_FP_woUS,
      # ptr_us_woUS_1 = ptr_us_woUS,
      prev_us_woUS_FP_1 = prev_us_woUS_FP, 
      
      TP_woCA125_1 = e_TP_woCA125,
      FP_woCA125_1 = e_FP_woCA125,
      prev_ca125_woCA125_FP_1 = prev_ca125_woCA125_FP, 
      
      uter_FP12_1 = e_uter_FP12,
      loGI_FP12_1 = e_loGI_FP12,
      lung_FP12_1 = e_lung_FP12,
      panc_FP12_1 = e_panc_FP12,
      uter_FP34_1 = e_uter_FP34,
      loGI_FP34_1 = e_loGI_FP34,
      lung_FP34_1 = e_lung_FP34,
      panc_FP34_1 = e_panc_FP34,
      cost_pc_1 = cost
    ) %>% 
    mutate(FP_true_1 = FP_1 - uter_FP12_1 - loGI_FP12_1 - lung_FP12_1 - panc_FP12_1
           - uter_FP34_1 - loGI_FP34_1 - lung_FP34_1 - panc_FP34_1) %>% 
    select(-e_TN, -e_uter_TN, -e_lung_TN, -e_loGI_TN,-e_panc_TN, -epatid)
  
  # stage shift
  comp <- merge(dt0, dt1) %>% 
    mutate(
      # late diagnosis in FN potentially moved to TP
      TP34_dif = TP34_1 - TP34_0,
      TP12_dif = TP12_1 - TP12_0,
      
      uter_FP12_dif = uter_FP12_1 - uter_FP12_0,
      loGI_FP12_dif = loGI_FP12_1 - loGI_FP12_0,
      lung_FP12_dif = lung_FP12_1 - lung_FP12_0,
      panc_FP12_dif = panc_FP12_1 - panc_FP12_0,
      
      uter_FP34_dif = uter_FP34_1 - uter_FP34_0,
      loGI_FP34_dif = loGI_FP34_1 - loGI_FP34_0,
      lung_FP34_dif = lung_FP34_1 - lung_FP34_0,
      panc_FP34_dif = panc_FP34_1 - panc_FP34_0,
      
      FP_dif = FP_1 - FP_0,
      FP_true_dif = FP_true_1 - FP_true_0,
      cost_pc_dif = cost_pc_1 - cost_pc_0
    ) %>% 
    mutate(
      # calculate stage shift
      # stage_change >0: from late to early
      # stage_change <0: from early to late
      
      # ova
      stage_change = case_when(
        
        # move from FN to TP, from late to early
        # only have effect on original late stage (34)
        TP34_dif>=0 & TP12_dif>=0 ~ TP34_dif - TP34_dif * shift_effects$ova_rr_late,
        # 
        # move from TP to FN, from early to late
        # only have effect on original early stage (12)
        TP34_dif<0 & TP12_dif<0 ~ -(abs(TP12_dif) - (abs(TP12_dif) / shift_effects$ova_rr_early)),
        # 
        # some move from early to late, some from late to early
        TP34_dif>=0 & TP12_dif<0 ~
          (TP34_dif - TP34_dif * shift_effects$ova_rr_late) -
          (abs(TP12_dif) - (abs(TP12_dif) / shift_effects$ova_rr_early)),
        # 
        # # no effect at all 
        TP34_dif<0 & TP12_dif>=0 ~ 0
      ),
      stage_early_upd = stage_early + stage_change,
      stage_late_upd = stage_late - stage_change
    ) %>% 
    mutate(
      ref_0 = TP12_0 + TP34_0 + FP_0,
      ref_1 = TP12_1 + TP34_1 + FP_1
    )
  
  # for other cancers
  oc_list <- c("uter", "loGI", "lung", "panc") 
  for (oc in oc_list){
    
    comp <- comp %>% 
      mutate(
        !!sym(paste0("stage_change_", oc)) := case_when(
          
          # move from FN to TP, from late to early
          # only have effect on original late stage (34)
          !!sym(paste0(oc, "_FP34_dif"))>=0 & !!sym(paste0(oc, "_FP12_dif"))>=0 ~
            !!sym(paste0(oc, "_FP34_dif")) - !!sym(paste0(oc, "_FP34_dif")) * shift_effects[[paste0(oc, "_rr_late")]],
          
          # move from TP to FN, from early to late
          # only have effect on original early stage (12)
          !!sym(paste0(oc, "_FP34_dif"))<0 & !!sym(paste0(oc, "_FP12_dif"))<0 ~
            -(abs(!!sym(paste0(oc, "_FP12_dif"))) - (abs(!!sym(paste0(oc, "_FP12_dif"))) / shift_effects[[paste0(oc, "_rr_early")]])),
          
          # some move from early to late, some from late to early  
          !!sym(paste0(oc, "_FP34_dif"))>=0 & !!sym(paste0(oc, "_FP12_dif"))<0 ~  
            (!!sym(paste0(oc, "_FP34_dif")) - !!sym(paste0(oc, "_FP34_dif")) * shift_effects[[paste0(oc, "_rr_late")]])-
            (abs(!!sym(paste0(oc, "_FP12_dif"))) - (abs(!!sym(paste0(oc, "_FP12_dif"))) / shift_effects[[paste0(oc, "_rr_early")]])), 
          
          # no effect at all 
          !!sym(paste0(oc, "_FP34_dif"))<0 & !!sym(paste0(oc, "_FP12_dif"))>=0 ~ 0
        )
      ) %>% 
      mutate(
        !!sym(paste0("stage_early_", oc, "_upd")) := !!sym(paste0("stage_early_", oc)) + !!sym(paste0("stage_change_", oc)),
        !!sym(paste0("stage_late_", oc, "_upd")) := !!sym(paste0("stage_late_", oc)) - !!sym(paste0("stage_change_", oc))
      )
  }
  
  # comp <- comp %>% 
  #   select(id, stage_early, stage_late, no_ova, stage_early_upd, stage_late_upd)
  
  return(comp)
}

# link the result from above to Markov results
# benign: Markov result for benign gynae disease, life expectancy equal to general population
# ova_early: Markov result for ovarian cancer patients diagnosed early
# ova_late: Markov result for ovarian cancer patients diagnosed late
# he_ova: linking to baseline characteristics
link_markov <- function(he_ova, comp, 
                        benign, 
                        ova_early, ova_late,
                        lung_early, lung_late,
                        panc_early, panc_late,
                        uter_early, uter_late,
                        loGI_early, loGI_late,
                        other_early, other_late,
                        FP_data, disuti, uti_ben, extra_cost,
                        surg_prob_adj=0
                        ){

  surg_prob <- FP_data$surg_prob + surg_prob_adj

  # link old id to new id
  id_dic <- readRDS(file.path(work_data, "id_dic.RDS"))
  
  benign <- benign %>% 
    rename(alive_ben = alive,
           cd_ben = cd,
           ncd_ben = ncd,
           qaly_ben = qaly, 
           cost_ben = cost, 
           alive_disc_ben = alive_disc,
           qaly_disc_ben = qaly_disc, 
           cost_disc_ben = cost_disc
           ) %>% 
    select(-age)
  
  c_list <- list(benign = benign)
  
  cancers <- c("ova", "lung", "panc", "uter", "loGI", "other")
  for (can in cancers) {
    for (stage in c("early", "late")){
      x <- get(paste0(can, "_", stage)) %>% 
        rename_with( ~ paste0(., "_", stage, "_", can), 
                     c(alive, cd, ncd, qaly, cost, alive_disc, qaly_disc, cost_disc)) %>% 
        select(-age)
      
      c_list[[paste0(can, "_", stage)]] <- x
    }
  }

  bsl <- he_ova %>% 
    left_join(id_dic) %>% 
    select(id,age, town5, ethn3, ethn4, age_group10, age_group_40_80, CA125)
  
  # to suit simulation with just a sample  
  dt <- comp %>% 
    right_join(bsl) %>% 
    right_join(c_list[["benign"]]) %>% 
    right_join(c_list[["ova_early"]]) %>% 
    right_join(c_list[["ova_late"]]) %>% 
    right_join(c_list[["lung_early"]]) %>% 
    right_join(c_list[["lung_late"]]) %>% 
    right_join(c_list[["panc_early"]]) %>% 
    right_join(c_list[["panc_late"]]) %>% 
    right_join(c_list[["uter_early"]]) %>% 
    right_join(c_list[["uter_late"]]) %>% 
    right_join(c_list[["loGI_early"]]) %>% 
    right_join(c_list[["loGI_late"]]) %>% 
    right_join(c_list[["other_early"]]) %>% 
    right_join(c_list[["other_late"]])
  
  
  # derive the fixed risk of other cancers from he_ova
  # assume missing stage as advanced
  # exclude those with ova
  #> use individual prevalence data instead 
  # n <- nrow(he_ova[he_ova[["cancer_ova_1yr"]]==0,])
  # c_rates <- c()
  # for (can in setdiff(cancers, "ova")) {
  # 
  #   x <- nrow(he_ova[he_ova[[paste0("cancer_", can, "_1yr")]]==1 & 
  #                      he_ova[[paste0("stage_", can, "_imp")]] %in% c("stage 1", "stage 2") &
  #                      he_ova[["cancer_ova_1yr"]]==0,])
  #   
  #   y <- nrow(he_ova[he_ova[[paste0("cancer_", can, "_1yr")]]==1 & 
  #                      (he_ova[[paste0("stage_", can, "_imp")]] %in% c("stage 3", "stage 4") | 
  #                         is.na(he_ova[[paste0("stage_", can, "_imp")]])) &
  #                      he_ova[["cancer_ova_1yr"]]==0,])
  #   
  #   c_rates[paste0(can, "_rate_early")] <- x / n
  #   c_rates[paste0(can, "_rate_late")] <- y / n
  # }
  # 
  # c_rates["all_cancer"] <- sum(c_rates)
  
  oc_list <- c("ly", "qaly", "cost", "ly_disc", "qaly_disc", "cost_disc", "cd", "ncd")

  # generate the non-cancer rate
  dt <- dt %>% 
    mutate(
      non_cancer_rate = no_ova - stage_early_lung - stage_late_lung - 
        stage_early_panc - stage_late_panc - 
        stage_early_uter - stage_late_uter -
        stage_early_loGI - stage_late_loGI -
        stage_early_other - stage_late_other,
      non_cancer_rate_FP_0 = (FP_0 - 
        uter_FP12_0 - loGI_FP12_0 - lung_FP12_0 - panc_FP12_0 - 
        uter_FP34_0 - loGI_FP34_0 - lung_FP34_0 - panc_FP34_0 -
        # not estimate other cancer in FP/TN, so assume proportional
        (FP_0/no_ova)*(stage_early_other + stage_late_other))/FP_0,
      non_cancer_rate_FP_1 = (FP_1 - 
        uter_FP12_1 - loGI_FP12_1 - lung_FP12_1 - panc_FP12_1 - 
        uter_FP34_1 - loGI_FP34_1 - lung_FP34_1 - panc_FP34_1 -
        # not estimate other cancer in FP/TN, so assume proportional
        (FP_1/no_ova)*(stage_early_other + stage_late_other))/FP_1
    )
  

  for (oc in oc_list) {
    
    if (oc=="ly") oc2 <- "alive" else 
      if (oc=="ly_disc") oc2 <- "alive_disc" else
        oc2 <- oc
    
    dt <- dt %>% 
      mutate(!!sym(paste0(oc, "_fixed_others")) := 
               !!sym(paste0(oc2, "_early_lung"))* stage_early_lung+
               !!sym(paste0(oc2, "_late_lung"))* stage_late_lung+
               !!sym(paste0(oc2, "_early_panc"))* stage_early_panc+
               !!sym(paste0(oc2, "_late_panc"))* stage_late_panc+
               !!sym(paste0(oc2, "_early_uter"))* stage_early_uter+
               !!sym(paste0(oc2, "_late_uter"))* stage_late_uter+
               !!sym(paste0(oc2, "_early_loGI"))* stage_early_loGI+
               !!sym(paste0(oc2, "_late_loGI"))* stage_late_loGI+
               !!sym(paste0(oc2, "_early_other"))* stage_early_other+
               !!sym(paste0(oc2, "_late_other"))* stage_late_other+
               !!sym(paste0(oc2, "_ben")) * non_cancer_rate
      ) %>% 
      mutate(
        # orignal 
        !!sym(paste0(oc, "_0")) := !!sym(paste0(oc2, "_early_ova"))*stage_early + 
          !!sym(paste0(oc2, "_late_ova"))*stage_late + 
          !!sym(paste0(oc, "_fixed_others")), 
        # update, with all non-ovarian cancer fixed, i.e. no stage shifting considered 
        !!sym(paste0(oc, "_1")) := !!sym(paste0(oc2, "_early_ova"))*stage_early_upd + 
          !!sym(paste0(oc2, "_late_ova"))*stage_late_upd + 
          !!sym(paste0(oc, "_fixed_others")),
        # update, with lung, panc, uter, loGI stage shifting considered
        # uo: update others
        !!sym(paste0(oc, "_uo_1")) := !!sym(paste0(oc2, "_early_ova"))*stage_early_upd + 
          !!sym(paste0(oc2, "_late_ova"))*stage_late_upd + 
          !!sym(paste0(oc2, "_early_lung"))*stage_early_lung_upd + 
          !!sym(paste0(oc2, "_late_lung"))*stage_late_lung_upd + 
          !!sym(paste0(oc2, "_early_panc"))*stage_early_panc_upd + 
          !!sym(paste0(oc2, "_late_panc"))*stage_late_panc_upd + 
          !!sym(paste0(oc2, "_early_uter"))*stage_early_uter_upd + 
          !!sym(paste0(oc2, "_late_uter"))*stage_late_uter_upd + 
          !!sym(paste0(oc2, "_early_loGI"))*stage_early_loGI_upd + 
          !!sym(paste0(oc2, "_late_loGI"))*stage_late_loGI_upd + 
          # for "other" cancers and benign, the same as fixed
          !!sym(paste0(oc2, "_early_other"))* stage_early_other+
          !!sym(paste0(oc2, "_late_other"))* stage_late_other+
          !!sym(paste0(oc2, "_ben")) * non_cancer_rate
      )
  }
  
  dt <- dt %>%
    mutate(
      cost_tot_0 = cost_0 + cost_pc_0, 
      cost_tot_1 = cost_1 + cost_pc_1,
      cost_tot_uo_1 = cost_uo_1 + cost_pc_1,
      
      cost_disc_tot_0 = cost_disc_0 + cost_pc_0, 
      cost_disc_tot_1 = cost_disc_1 + cost_pc_1,
      cost_disc_tot_uo_1 = cost_disc_uo_1 + cost_pc_1
    )

  # apply extra utility and cost
  dt <- dt %>% 
    # consider disutility of surgery for FP patients
    mutate(
      # Buys, Partridge et al. JAMA, 2011, 305(22) 2295-2303
      # of 3285 FP, 1080 underwent surgery
      # Oxley, Wei et al. Cancers, 2025, 16, 1358
      # adjusted disutility was -0.04 (-0.06, -0.01) within 1 year after Risk_reducing Salpingo-Oophorectomy 
      # FP_woUS_0: FP without US+ in primary
      # prev_us_woUS_FP_0: probability of FP without US+ in primary but US+ in secondary
      # so (FP_woUS_0 - prev_us_woUS_FP_0): FP without US+ in primary and secondary
      # FP_0_surg = (FP_0 - (FP_woUS_0 - prev_us_woUS_FP_0)) * non_cancer_rate * surg_prob,
      # FP_1_surg = (FP_1 - (FP_woUS_1 - prev_us_woUS_FP_1)) * non_cancer_rate * surg_prob,
      
      FP_0_surg = (FP_0 * non_cancer_rate_FP_0 - (FP_woUS_0 - prev_us_woUS_FP_0) - (FP_woCA125_0 - prev_ca125_woCA125_FP_0)) * surg_prob,
      FP_1_surg = (FP_1 * non_cancer_rate_FP_1 - (FP_woUS_1 - prev_us_woUS_FP_1) - (FP_woCA125_1 - prev_ca125_woCA125_FP_1)) * surg_prob,
      
      # # do not consider CA125 in secondary care
      # FP_0_surg = (FP_0 * non_cancer_rate_FP_0 - (FP_woUS_0 - prev_us_woUS_FP_0)) * surg_prob,
      # FP_1_surg = (FP_1 * non_cancer_rate_FP_1 - (FP_woUS_1 - prev_us_woUS_FP_1)) * surg_prob,
      
      # FP_0_disuti = FP_0_surg * (FP_data$surg_disu + alive_ben * FP_data$surg_ben),
      # FP_1_disuti = FP_1_surg * (FP_data$surg_disu + alive_ben * FP_data$surg_ben),
      # FP_0_disuti_disc = FP_0_surg * (FP_data$surg_disu + alive_disc_ben * FP_data$surg_ben),
      # FP_1_disuti_disc = FP_1_surg * (FP_data$surg_disu + alive_disc_ben * FP_data$surg_ben),
      
      # base case: consider -0.04 surgery harm in the first year
      FP_0_disuti = FP_0_surg * FP_data$surg_disu,
      FP_1_disuti = FP_1_surg * FP_data$surg_disu,

      FP_0_ben = FP_0_surg * (alive_ben * FP_data$surg_ben),
      FP_1_ben = FP_1_surg * (alive_ben * FP_data$surg_ben),
      FP_0_ben_disc = FP_0_surg * (alive_disc_ben * FP_data$surg_ben),
      FP_1_ben_disc = FP_1_surg * (alive_disc_ben * FP_data$surg_ben),
      
      # NHS tariff 21/22 HRG code: MB09F non-malignant gynaecological disorders with interventions
      # FP_0_cost = (FP_0 - (FP_woUS_0 - prev_us_woUS_FP_0)) * non_cancer_rate * 
      #   (surg_prob * FP_data$surg_cost + (1-surg_prob) * FP_data$no_surg_cost) +
      #   (FP_woUS_0 - prev_us_woUS_FP_0) * FP_data$no_surg_cost,
      # FP_1_cost = (FP_1 - (FP_woUS_1 - prev_us_woUS_FP_1)) * non_cancer_rate * 
      #   (surg_prob * FP_data$surg_cost + (1-surg_prob) * FP_data$no_surg_cost) +
      #   (FP_woUS_1 - prev_us_woUS_FP_1) * FP_data$no_surg_cost
      
      FP_0_cost = (FP_0 * non_cancer_rate_FP_0 - (FP_woUS_0 - prev_us_woUS_FP_0) - (FP_woCA125_0 - prev_ca125_woCA125_FP_0)) * 
        (surg_prob * FP_data$surg_cost + (1-surg_prob) * FP_data$no_surg_cost) +
        (FP_woUS_0 - prev_us_woUS_FP_0 + FP_woCA125_0 - prev_ca125_woCA125_FP_0) * FP_data$no_surg_cost,
      FP_1_cost = (FP_1 * non_cancer_rate_FP_1 - (FP_woUS_1 - prev_us_woUS_FP_1) - (FP_woCA125_1 - prev_ca125_woCA125_FP_1)) * 
        (surg_prob * FP_data$surg_cost + (1-surg_prob) * FP_data$no_surg_cost) +
        (FP_woUS_1 - prev_us_woUS_FP_1 + FP_woCA125_1 - prev_ca125_woCA125_FP_1) * FP_data$no_surg_cost
    )
  
  # consider CT cost if consider pathway effects on other cancers 
  # only apply to lung, lower GI, pancreatic, uterin cancers in FP
  dt <- dt %>% 
    mutate(
      ct_cost_0 = (uter_FP12_0 + loGI_FP12_0 + lung_FP12_0 + panc_FP12_0 + 
        uter_FP34_0 + loGI_FP34_0 + lung_FP34_0 + panc_FP34_0) * FP_data$ct_cost,
      
      ct_cost_1 = (uter_FP12_1 + loGI_FP12_1 + lung_FP12_1 + panc_FP12_1 + 
        uter_FP34_1 + loGI_FP34_1 + lung_FP34_1 + panc_FP34_1) * FP_data$ct_cost,
      
      ct_cost_dif_uo = ct_cost_1 - ct_cost_0
    )
  
  
  if (disuti){
    dt <- dt %>% 
      mutate(
        qaly_0 = qaly_0 + FP_0_disuti,
        qaly_1 = qaly_1 +  FP_1_disuti,
        qaly_uo_1 = qaly_uo_1 + FP_1_disuti,
        
        qaly_disc_0 = qaly_disc_0 + FP_0_disuti,
        qaly_disc_1 = qaly_disc_1 + FP_1_disuti,
        qaly_disc_uo_1 = qaly_disc_uo_1 + FP_1_disuti,
      )
  }
  
  if (uti_ben){
    dt <- dt %>% 
      mutate(
        qaly_0 = qaly_0 + FP_0_ben,
        qaly_1 = qaly_1 +  FP_1_ben,
        qaly_uo_1 = qaly_uo_1 + FP_1_ben,
        
        qaly_disc_0 = qaly_disc_0 + FP_0_ben_disc,
        qaly_disc_1 = qaly_disc_1 + FP_1_ben_disc,
        qaly_disc_uo_1 = qaly_disc_uo_1 + FP_1_ben_disc
      )
  }
  
  if (extra_cost){
    dt <- dt %>% 
      mutate(
        cost_tot_0 = cost_tot_0 + FP_0_cost,
        cost_tot_1 = cost_tot_1 + FP_1_cost,
        cost_tot_uo_1 = cost_tot_uo_1 + FP_1_cost + ct_cost_dif_uo,
        
        cost_disc_tot_0 = cost_disc_tot_0 + FP_0_cost,
        cost_disc_tot_1 = cost_disc_tot_1 + FP_1_cost,
        cost_disc_tot_uo_1 = cost_disc_tot_uo_1 + FP_1_cost + ct_cost_dif_uo
      )
  }
  
  dt <- dt %>% 
    mutate(
      ly_dif = ly_1 - ly_0,
      qaly_dif = qaly_1 - qaly_0,
      cost_tot_dif = cost_tot_1 - cost_tot_0,
      
      ly_disc_dif = ly_disc_1 - ly_disc_0,
      qaly_disc_dif = qaly_disc_1 - qaly_disc_0,
      cost_disc_tot_dif = cost_disc_tot_1 - cost_disc_tot_0,
      
      cd_dif = cd_1 - cd_0,
      ncd_dif = ncd_1 - ncd_0,
      ova_prev = stage_early + stage_late,
      # sub costs
      cost_dif = cost_1 - cost_0, # hosptial cost only
      cost_disc_dif = cost_disc_1 - cost_disc_0,
      cost_dif_fromFP = FP_1_cost - FP_0_cost,
      
      # update, with lung, panc, uter, loGI stage shifting considered 
      ly_uo_dif = ly_uo_1 - ly_0,
      qaly_uo_dif = qaly_uo_1 - qaly_0,
      cost_tot_uo_dif = cost_tot_uo_1 - cost_tot_0,
      
      ly_disc_uo_dif = ly_disc_uo_1 - ly_disc_0,
      qaly_disc_uo_dif = qaly_disc_uo_1 - qaly_disc_0,
      cost_disc_tot_uo_dif = cost_disc_tot_uo_1 - cost_disc_tot_0,
      
      cost_uo_dif = cost_uo_1 - cost_0, # hosptial cost only
      cost_disc_uo_dif = cost_disc_uo_1 - cost_disc_0,
      cd_uo_dif = cd_uo_1 - cd_0,
      ncd_uo_dif = ncd_uo_1 - ncd_0
    ) 
  
  dt <- dt%>%
    select(id, town5, ethn3, ethn4, age, age_group10, age_group_40_80, CA125,
           ova_prev, stage_early, stage_late, TP12_0, TP34_0, TP12_1, TP34_1, 
           stage_late_upd, stage_early_upd, stage_change, 
           FP_0, FP_dif, FP_true_dif, ref_0, ref_1, 
           ly_0, ly_1, ly_uo_1, qaly_0, qaly_1, qaly_uo_1,
           qaly_disc_0, qaly_disc_1, qaly_disc_uo_1,
           cost_disc_tot_0, cost_disc_tot_1, cost_disc_tot_uo_1,
           cost_0, cost_1, FP_0_cost, FP_1_cost, cost_uo_1,
           cost_tot_0, cost_tot_1, cost_tot_uo_1,
           cd_0, cd_1, cd_uo_1, ncd_0, ncd_1, ncd_uo_1,
           FP_0_disuti, FP_1_disuti, FP_0_surg, FP_1_surg,

           non_cancer_rate, non_cancer_rate_FP_0, non_cancer_rate_FP_1,
           ly_dif, ly_uo_dif, qaly_dif, qaly_uo_dif,
           cost_dif, cost_uo_dif, cost_tot_dif, cost_tot_uo_dif, cost_pc_0, cost_pc_dif,
           ly_disc_dif, ly_disc_uo_dif, qaly_disc_dif, qaly_disc_uo_dif,
           cost_disc_dif, cost_disc_uo_dif, cost_disc_tot_dif, cost_disc_tot_uo_dif,
           cd_dif, ncd_dif, cd_uo_dif, cost_dif_fromFP)

  dt <- dt %>% 
    mutate(age_group_40_70 = if_else(age_group_40_80 %in% c("70-79", ">=80"), 
                                     ">=70", age_group_40_80),
           age_group_40_70 = factor(age_group_40_70, 
                                    levels=c("<40", "40-49", "50-59","60-69", ">=70")),
           age_50 = if_else(age_group_40_80 %in% c("<40", "40-49"), 
                            "<50", ">=50"),
           age_50 = factor(age_50, levels=c("<50", ">=50"))
    ) %>% 
    mutate(town5 = factor(town5, levels=c("Q1", "Q2", "Q3", "Q4", "Q5"))) 
  
  
  return(dt)
}


# analytical results ------------------------------------------------------
# comp is the output of stage_shift_effect function
# dt is the output of link_markov function

# group="ethn3"

ana_rslt <- function(dt, age50filter=NULL, group=NULL){
  # age50filter can be "<50" or  ">=50"
  dtt <- dt
  if (!is.null(age50filter)) dtt <- dt %>% filter(age_50==age50filter)

  if (is.null(group)){
    rt <- dtt %>% 
      summarise(
        N = mean(n()),
        across(c("ova_prev","stage_late", "stage_early", "TP12_0", "TP34_0", 
                 "TP12_1", "TP34_1", "stage_late_upd", "stage_early_upd", 
                 "stage_change","FP_0", "FP_dif", "FP_true_dif", "ref_0", "ref_1", 
                 "ly_dif", "ly_uo_dif", "ly_disc_dif", "ly_disc_uo_dif", 
                 "qaly_dif", "qaly_uo_dif", "qaly_disc_dif", "qaly_disc_uo_dif",
                 "qaly_0", "qaly_1", "FP_0_disuti", "FP_1_disuti", "FP_0_surg", "FP_1_surg",
                 "qaly_disc_0", "qaly_disc_1", "qaly_disc_uo_1",
                 "cost_disc_tot_0", "cost_disc_tot_1", "cost_disc_tot_uo_1",
                 "cost_tot_dif","cost_tot_uo_dif", "cost_disc_tot_dif", "cost_disc_tot_uo_dif",
                 "cost_dif", "cost_uo_dif", "cost_pc_dif", "cost_dif_fromFP", 
                 "cost_tot_0", "cost_0", "cost_pc_0", "FP_0_cost"), 
               sum, .names = "sum_{col}")
      ) %>% 
      mutate(
        icer = sum_cost_disc_tot_dif/sum_qaly_disc_dif,
        icer_uo = sum_cost_disc_tot_uo_dif/sum_qaly_disc_uo_dif,
        icer_ly = sum_cost_disc_tot_dif/sum_ly_disc_dif,
        icer_ly_uo = sum_cost_disc_tot_uo_dif/sum_ly_disc_uo_dif,
        TP_pct_0 = (sum_TP12_0 + sum_TP34_0)/(sum_stage_early + sum_stage_late),
        TP12_pct_0 = sum_TP12_0/sum_stage_early,
        TP34_pct_0 = sum_TP34_0/sum_stage_late,
        TP_pct_1 = (sum_TP12_1 + sum_TP34_1)/(sum_stage_early + sum_stage_late),
        TP12_pct_1 = sum_TP12_1/sum_stage_early,
        TP34_pct_1 = sum_TP34_1/sum_stage_late,
        ref_pct_0 = sum_ref_0/N,
        ref_pct_1 = sum_ref_1/N,
        FP_pct_0 = sum_FP_0/sum_ref_0,
        FP_pct_1 = (sum_FP_0 + sum_FP_dif)/sum_ref_1,
        stage_change_pct=sum_stage_change/sum_stage_late
      ) 
  } else {
    # group <- deparse(substitute(group))
    rt <- dtt %>% 
      group_by(get(group)) %>% 
      summarise(
        N = mean(n()),
        across(c("ova_prev","stage_late", "stage_early", "TP12_0", "TP34_0", 
                 "TP12_1", "TP34_1", "stage_late_upd", "stage_early_upd", 
                 "stage_change","FP_0", "FP_dif", "FP_true_dif", "ref_0", "ref_1", 
                 "ly_dif", "ly_uo_dif", "ly_disc_dif", "ly_disc_uo_dif", 
                 "qaly_dif", "qaly_uo_dif", "qaly_disc_dif", "qaly_disc_uo_dif",
                 "qaly_0", "qaly_1", "FP_0_disuti", "FP_1_disuti", "FP_0_surg", "FP_1_surg",
                 "qaly_disc_0", "qaly_disc_1", "qaly_disc_uo_1",
                 "cost_disc_tot_0", "cost_disc_tot_1", "cost_disc_tot_uo_1",
                 "cost_tot_dif","cost_tot_uo_dif", "cost_disc_tot_dif", "cost_disc_tot_uo_dif",
                 "cost_dif", "cost_uo_dif", "cost_pc_dif", "cost_dif_fromFP", 
                 "cost_tot_0", "cost_0", "cost_pc_0", "FP_0_cost"), 
               sum, .names = "sum_{col}")
      ) %>%
      mutate(
        icer = sum_cost_disc_tot_dif/sum_qaly_disc_dif,
        icer_uo = sum_cost_disc_tot_uo_dif/sum_qaly_disc_uo_dif,
        icer_ly = sum_cost_disc_tot_dif/sum_ly_disc_dif,
        icer_ly_uo = sum_cost_disc_tot_uo_dif/sum_ly_disc_uo_dif,
        TP_pct_0 = (sum_TP12_0 + sum_TP34_0)/(sum_stage_early + sum_stage_late),
        TP12_pct_0 = sum_TP12_0/sum_stage_early,
        TP34_pct_0 = sum_TP34_0/sum_stage_late,
        TP_pct_1 = (sum_TP12_1 + sum_TP34_1)/(sum_stage_early + sum_stage_late),
        TP12_pct_1 = sum_TP12_1/sum_stage_early,
        TP34_pct_1 = sum_TP34_1/sum_stage_late,
        ref_pct_0 = sum_ref_0/N,
        ref_pct_1 = sum_ref_1/N,
        FP_pct_0 = sum_FP_0/sum_ref_0,
        FP_pct_1 = (sum_FP_0 + sum_FP_dif)/sum_ref_1,
        stage_change_pct=sum_stage_change/sum_stage_late
      )
  }
  
  return(rt)
  
}

# Map the result from simulation using the unique characteristics to the whole population
map2pop <- function(rst, id_dic, ipd2){
  
  ipd2 <- ipd2 %>% 
    # output result use the same psuedo ids as ipd2 here
    mutate(id=1:n()) %>% 
    relocate(id, .before = age_cent60)
  
  # use psuedo ids to link to ipd2 to add group id 
  rst2 <- rst %>% 
    left_join(ipd2[, c("id", "group")]) %>% 
    # get rid of pseudo ids
    select(-id)
  
  rt <- id_dic %>% 
    left_join(rst2)
  
  return(rt)
}









