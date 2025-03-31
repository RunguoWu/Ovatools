##################################################################
##         Master functions for decision tree PSA files         ##
##################################################################

# adapted from fn_master_decision_tree.R in HE_model folder

master_decision_tree <- function(i=1001, byAge = TRUE, accuracy_file, adj_pop = FALSE, 
                                 cost_us_cdc=FALSE, cost_us_op=FALSE, cost_us_input=NA,
                                 us_sens_adj=0, us_spec_adj=0){
  
  adj_pop_tag <- ""
  adj_addon <- ""
  if (adj_pop){
    adj_pop_tag <- "_adj_pop"
    adj_addon <- "_ca125"
  }
  # import prevalence and stage data
  prevalence <- readRDS(file.path(work_data, "prev/ova", paste0("ova_prev", adj_addon, i, ".rds")))
  ova_stage34 <- readRDS(file.path(work_data, "stage_rate/ova", paste0("ova_stage", adj_addon, i, ".rds")))
  
  pred_loGI <- readRDS(file.path(work_data, "prev/loGI", paste0("loGI_prev", adj_addon, i, ".rds")))
  loGI_stage34 <- readRDS(file.path(work_data, "stage_rate/loGI", paste0("loGI_stage", adj_addon, i, ".rds")))
  
  pred_uter <- readRDS(file.path(work_data, "prev/uter", paste0("uter_prev", adj_addon, i, ".rds")))
  uter_stage34 <- readRDS(file.path(work_data, "stage_rate/uter", paste0("uter_stage", adj_addon, i, ".rds")))
  
  pred_lung <- readRDS(file.path(work_data, "prev/lung", paste0("lung_prev", adj_addon, i, ".rds")))
  lung_stage34 <- readRDS(file.path(work_data, "stage_rate/lung", paste0("lung_stage", adj_addon, i, ".rds")))
  
  pred_panc <- readRDS(file.path(work_data, "prev/panc", paste0("panc_prev", adj_addon, i, ".rds")))
  panc_stage34 <- readRDS(file.path(work_data, "stage_rate/panc", paste0("panc_stage", adj_addon, i, ".rds")))
  
  pred_other <- readRDS(file.path(work_data, "prev/other", paste0("other_prev", adj_addon, i, ".rds")))
  other_stage34 <- readRDS(file.path(work_data, "stage_rate/other", paste0("other_stage", adj_addon, i, ".rds")))

  
  id_new <- cbind(id_dic2, prevalence, ova_stage34, pred_loGI, loGI_stage34,
                  pred_uter, uter_stage34, pred_lung, lung_stage34,
                  pred_panc, panc_stage34, pred_other, other_stage34)
  
  id_new <- id_new %>% 
    distinct(group, .keep_all = T) %>% 
    select(-id)
  
  IPD_i <- IPD %>% 
    left_join(id_new) # order of group is unchanged from ipd2
  
  # adjust prevalence to the population with ca125 
  # adj_pop_tag <- ""
  # if (adj_pop){
  #   IPD_i <- IPD_i %>% mutate(
  #     prevalence = prevalence * 1.3,
  #     pred_loGI = pred_loGI * 1.3,
  #     pred_lung = pred_lung * 1.3,
  #     pred_panc = pred_panc * 1.3,
  #     # uter is the same
  #     pred_other = pred_other * 1.2
  #   )
  #   adj_pop_tag <- "_adj_pop"
  # }
  
  cost_us_cdc_tag <- if (cost_us_cdc) "_cost_us_cdc" else  ""
  cost_us_op_tag <- if (cost_us_op) "_cost_us_op" else  ""
  cost_us_input_tag <- if(!is.na(cost_us_input)) paste0("_cost_us", cost_us_input) else ""
  
  us_sens_adj_tag <- ""
  if (us_sens_adj > 0 ) us_sens_adj_tag <- "_USsensUP" 
  if (us_sens_adj < 0 ) us_sens_adj_tag <- "_USsensDOWN" 
  
  us_spec_adj_tag <- ""
  if (us_spec_adj > 0 ) us_spec_adj_tag <- "_USspecUP" 
  if (us_spec_adj < 0 ) us_spec_adj_tag <- "_USspecDOWN" 
  
  # accuracy data
  if (byAge) {
    
    IPD_i_u49 <- IPD_i %>% filter(age_cent60 < -1)
    IPD_i_o50 <- IPD_i %>% filter(age_cent60 >= -1)
    
    accuracy_tot <- readRDS(file.path(work_data, paste0("accuracy_psa_", accuracy_file, ".rds")))[[i]]
    
    accuracy_tot_u49 <- accuracy_tot[["u49"]]
    accuracy_tot_o50 <- accuracy_tot[["o50"]]
    
    # CA125, current
    u49 <- decision_tree(ca125_first=1, us_first=0, ovatools=0, dt=IPD_i_u49,
                         accuracy_tot=accuracy_tot_u49, 
                         cost_us_cdc=cost_us_cdc, cost_us_op=cost_us_op, cost_us_input=cost_us_input,
                         us_sens_adj=us_sens_adj, us_spec_adj=us_spec_adj)
    o50 <- decision_tree(ca125_first=1, us_first=0, ovatools=0, dt=IPD_i_o50,
                         accuracy_tot=accuracy_tot_o50, 
                         cost_us_cdc=cost_us_cdc, cost_us_op=cost_us_op, cost_us_input=cost_us_input,
                         us_sens_adj=us_sens_adj, us_spec_adj=us_spec_adj)
    x <- rbind(u49, o50) %>% select(-id)
    decitree100 <- IPD_ids %>% left_join(x)
    saveRDS(decitree100, file.path(output, "decitree100",
                                   paste0("decitree100_", accuracy_file, i, adj_pop_tag, 
                                          cost_us_cdc_tag, cost_us_op_tag, cost_us_input_tag,
                                          us_sens_adj_tag, us_spec_adj_tag, ".RDS")))

    # CA125, concurrent
    u49 <- decision_tree(ca125_first=1, us_first=1, ovatools=0, dt=IPD_i_u49,
                         accuracy_tot=accuracy_tot_u49, 
                         cost_us_cdc=cost_us_cdc, cost_us_op=cost_us_op, cost_us_input=cost_us_input,
                         us_sens_adj=us_sens_adj, us_spec_adj=us_spec_adj)
    o50 <- decision_tree(ca125_first=1, us_first=1, ovatools=0, dt=IPD_i_o50,
                         accuracy_tot=accuracy_tot_o50, 
                         cost_us_cdc=cost_us_cdc, cost_us_op=cost_us_op, cost_us_input=cost_us_input,
                         us_sens_adj=us_sens_adj, us_spec_adj=us_spec_adj)
    x <- rbind(u49, o50) %>% select(-id)
    decitree110 <- IPD_ids %>% left_join(x)
    saveRDS(decitree110, file.path(output, "psa/dectr", "decitree110",
                                   paste0("decitree110_", accuracy_file, i, adj_pop_tag, 
                                          cost_us_cdc_tag, cost_us_op_tag, cost_us_input_tag,
                                          us_sens_adj_tag, us_spec_adj_tag, ".RDS")))

    # Ovatools sequential
    u49 <- decision_tree(ca125_first=1, us_first=0, ovatools=1, dt=IPD_i_u49,
                         accuracy_tot=accuracy_tot_u49, 
                         cost_us_cdc=cost_us_cdc, cost_us_op=cost_us_op, cost_us_input=cost_us_input,
                         us_sens_adj=us_sens_adj, us_spec_adj=us_spec_adj)
    o50 <- decision_tree(ca125_first=1, us_first=0, ovatools=1, dt=IPD_i_o50,
                         accuracy_tot=accuracy_tot_o50, 
                         cost_us_cdc=cost_us_cdc, cost_us_op=cost_us_op, cost_us_input=cost_us_input,
                         us_sens_adj=us_sens_adj, us_spec_adj=us_spec_adj)
    x <- rbind(u49, o50) %>% select(-id)
    decitree101 <- IPD_ids %>% left_join(x)
    saveRDS(decitree101, file.path(output, "psa/dectr", "decitree101",
                                   paste0("decitree101_", accuracy_file, i, adj_pop_tag, 
                                          cost_us_cdc_tag, cost_us_op_tag, cost_us_input_tag,
                                          us_sens_adj_tag, us_spec_adj_tag, ".RDS")))

    # Ovatools concurrent
    u49 <- decision_tree(ca125_first=1, us_first=1, ovatools=2, dt=IPD_i_u49,
                         accuracy_tot=accuracy_tot_u49, 
                         cost_us_cdc=cost_us_cdc, cost_us_op=cost_us_op, cost_us_input=cost_us_input,
                         us_sens_adj=us_sens_adj, us_spec_adj=us_spec_adj)
    o50 <- decision_tree(ca125_first=1, us_first=1, ovatools=2, dt=IPD_i_o50,
                         accuracy_tot=accuracy_tot_o50, 
                         cost_us_cdc=cost_us_cdc, cost_us_op=cost_us_op, cost_us_input=cost_us_input,
                         us_sens_adj=us_sens_adj, us_spec_adj=us_spec_adj)
    x <- rbind(u49, o50) %>% select(-id)
    decitree112 <- IPD_ids %>% left_join(x)
    saveRDS(decitree112, file.path(output, "psa/dectr", "decitree112",
                                   paste0("decitree112_", accuracy_file, i, adj_pop_tag, 
                                          cost_us_cdc_tag, cost_us_op_tag, cost_us_input_tag,
                                          us_sens_adj_tag, us_spec_adj_tag, ".RDS")))

    # CA125, varying threshold, sequential
    accuracy_tot_var <- readRDS(file.path(work_data, "psa", "accuracy_psa_CA125varyingThreshold_niceUS.rds"))[[i]]
    # accuracy_tot_var <- readRDS(file.path(work_data, "psa", "accuracy_psa_CA125varyingThreshold_2stepUS.rds"))[[i]]

    accuracy_tot_u49_var <- accuracy_tot_var[["u49"]]
    accuracy_tot_o50_var <- accuracy_tot_var[["o50"]]

    # use the 101 pathway, as Ovatools sequential. 
    # accuracy has been replaced with CA125 varying threshold file
    u49 <- decision_tree(ca125_first=1, us_first=0, ovatools=1, dt=IPD_i_u49,
                         accuracy_tot=accuracy_tot_u49_var, 
                         cost_us_cdc=cost_us_cdc, cost_us_op=cost_us_op, cost_us_input=cost_us_input,
                         us_sens_adj=us_sens_adj, us_spec_adj=us_spec_adj)
    o50 <- decision_tree(ca125_first=1, us_first=0, ovatools=1, dt=IPD_i_o50,
                         accuracy_tot=accuracy_tot_o50_var, 
                         cost_us_cdc=cost_us_cdc, cost_us_op=cost_us_op, cost_us_input=cost_us_input,
                         us_sens_adj=us_sens_adj, us_spec_adj=us_spec_adj)
    x <- rbind(u49, o50) %>% select(-id)
    decitree101var <- IPD_ids %>% left_join(x)

    saveRDS(decitree101var, file.path(
      output, "psa/dectr", "decitree101var",
      paste0("decitree101var_", "CA125varyingThreshold", i, adj_pop_tag, 
             cost_us_cdc_tag, cost_us_op_tag, cost_us_input_tag,
             us_sens_adj_tag, us_spec_adj_tag, ".RDS")))
      # paste0("decitree101var_", "CA125varyingThreshold_2stepUS", i, adj_pop_tag, cost_us_cdc_tag, ".RDS")))
    
    
    # CA125, varying threshold, concurrent
    accuracy_tot_var <- readRDS(file.path(work_data, "psa", "accuracy_psa_CA125varyingThreshold_niceUS.rds"))[[i]]
    # accuracy_tot_var <- readRDS(file.path(work_data, "psa", "accuracy_psa_CA125varyingThreshold_2stepUS.rds"))[[i]]
    
    accuracy_tot_u49_var <- accuracy_tot_var[["u49"]]
    accuracy_tot_o50_var <- accuracy_tot_var[["o50"]]
    
    # use the 112 pathway, as Ovatools concurrent 
    # accuracy has been replaced with CA125 varying threshold 3%
    u49 <- decision_tree(ca125_first=1, us_first=1, ovatools=2, dt=IPD_i_u49,
                         accuracy_tot=accuracy_tot_u49_var, 
                         cost_us_cdc=cost_us_cdc, cost_us_op=cost_us_op, cost_us_input=cost_us_input,
                         us_sens_adj=us_sens_adj, us_spec_adj=us_spec_adj)
    o50 <- decision_tree(ca125_first=1, us_first=1, ovatools=2, dt=IPD_i_o50,
                         accuracy_tot=accuracy_tot_o50_var, 
                         cost_us_cdc=cost_us_cdc, cost_us_op=cost_us_op, cost_us_input=cost_us_input,
                         us_sens_adj=us_sens_adj, us_spec_adj=us_spec_adj)
    x <- rbind(u49, o50) %>% select(-id)
    decitree112var <- IPD_ids %>% left_join(x)
    
    saveRDS(decitree112var, file.path(
      output, "psa/dectr", "decitree112var",
      paste0("decitree112var_", "CA125varyingThreshold", i, adj_pop_tag, 
             cost_us_cdc_tag, cost_us_op_tag, cost_us_input_tag,
             us_sens_adj_tag, us_spec_adj_tag, ".RDS")))

  } else {
    accuracy_tot <- readRDS(file.path(work_data, "psa", paste0("accuracy_psa_", accuracy_file, ".rds")))[[i]]
    
    run the model
    decitree100 <- decision_tree(ca125_first=1, us_first=0, ovatools=0, dt=IPD_i,
                                 accuracy_tot=accuracy_tot,
                                 cost_us_cdc=cost_us_cdc, cost_us_op=cost_us_op,
                                 us_sens_adj=us_sens_adj, us_spec_adj=us_spec_adj)
    saveRDS(decitree100, file.path(output, "psa/dectr", "decitree100",
                                   paste0("decitree100_", accuracy_file, i, adj_pop_tag,
                                          cost_us_cdc_tag, cost_us_op_tag,
                                          us_sens_adj_tag, us_spec_adj_tag, ".RDS")))

    decitree110 <- decision_tree(ca125_first=1, us_first=1, ovatools=0, dt=IPD_i,
                                 accuracy_tot=accuracy_tot,
                                 cost_us_cdc=cost_us_cdc, cost_us_op=cost_us_op,
                                 us_sens_adj=us_sens_adj, us_spec_adj=us_spec_adj)
    saveRDS(decitree110, file.path(output, "psa/dectr", "decitree110",
                                   paste0("decitree110_", accuracy_file, i,
                                          adj_pop_tag, cost_us_cdc_tag,
                                          us_sens_adj_tag, us_spec_adj_tag, ".RDS")))

    decitree101 <- decision_tree(ca125_first=1, us_first=0, ovatools=1, dt=IPD_i,
                                 accuracy_tot=accuracy_tot,
                                 cost_us_cdc=cost_us_cdc, cost_us_op=cost_us_op,
                                 us_sens_adj=us_sens_adj, us_spec_adj=us_spec_adj)
    saveRDS(decitree101, file.path(output, "psa/dectr", "decitree101",
                                   paste0("decitree101_", accuracy_file, i, adj_pop_tag,
                                          cost_us_cdc_tag, cost_us_op_tag,
                                          us_sens_adj_tag, us_spec_adj_tag, ".RDS")))

    decitree112 <- decision_tree(ca125_first=1, us_first=1, ovatools=2, dt=IPD_i,
                                 accuracy_tot=accuracy_tot,
                                 cost_us_cdc=cost_us_cdc, cost_us_op=cost_us_op,
                                 us_sens_adj=us_sens_adj, us_spec_adj=us_spec_adj)
    saveRDS(decitree112, file.path(output, "psa/dectr", "decitree112",
                                   paste0("decitree112_", accuracy_file, i, adj_pop_tag,
                                          cost_us_cdc_tag, cost_us_op_tag,
                                          us_sens_adj_tag, us_spec_adj_tag, ".RDS")))
  }
}



# Current and Ovatools sequential only ------------------------------------

# for the ICER test only

master_decision_tree2 <- function(i=1001, j, byAge = TRUE, accuracy_file, adj_pop = FALSE, 
                                 cost_us_cdc=FALSE, cost_us_op=FALSE, cost_us_input=NA,
                                 us_sens_adj=0, us_spec_adj=0){
  
  adj_pop_tag <- ""
  adj_addon <- ""
  if (adj_pop){
    adj_pop_tag <- "_adj_pop"
    adj_addon <- "_ca125"
  }
  # import prevalence and stage data
  prevalence <- readRDS(file.path(work_data, "psa/prev/ova", paste0("ova_prev", adj_addon, i, ".rds")))
  ova_stage34 <- readRDS(file.path(work_data, "psa/stage_rate/ova", paste0("ova_stage", adj_addon, i, ".rds")))
  
  pred_loGI <- readRDS(file.path(work_data, "psa/prev/loGI", paste0("loGI_prev", adj_addon, i, ".rds")))
  loGI_stage34 <- readRDS(file.path(work_data, "psa/stage_rate/loGI", paste0("loGI_stage", adj_addon, i, ".rds")))
  
  pred_uter <- readRDS(file.path(work_data, "psa/prev/uter", paste0("uter_prev", adj_addon, i, ".rds")))
  uter_stage34 <- readRDS(file.path(work_data, "psa/stage_rate/uter", paste0("uter_stage", adj_addon, i, ".rds")))
  
  pred_lung <- readRDS(file.path(work_data, "psa/prev/lung", paste0("lung_prev", adj_addon, i, ".rds")))
  lung_stage34 <- readRDS(file.path(work_data, "psa/stage_rate/lung", paste0("lung_stage", adj_addon, i, ".rds")))
  
  pred_panc <- readRDS(file.path(work_data, "psa/prev/panc", paste0("panc_prev", adj_addon, i, ".rds")))
  panc_stage34 <- readRDS(file.path(work_data, "psa/stage_rate/panc", paste0("panc_stage", adj_addon, i, ".rds")))
  
  pred_other <- readRDS(file.path(work_data, "psa/prev/other", paste0("other_prev", adj_addon, i, ".rds")))
  other_stage34 <- readRDS(file.path(work_data, "psa/stage_rate/other", paste0("other_stage", adj_addon, i, ".rds")))
  
  
  id_new <- cbind(id_dic2, prevalence, ova_stage34, pred_loGI, loGI_stage34,
                  pred_uter, uter_stage34, pred_lung, lung_stage34,
                  pred_panc, panc_stage34, pred_other, other_stage34)
  
  id_new <- id_new %>% 
    distinct(group, .keep_all = T) %>% 
    select(-id)
  
  IPD_i <- IPD %>% 
    left_join(id_new) # order of group is unchanged from ipd2

  cost_us_cdc_tag <- if (cost_us_cdc) "_cost_us_cdc" else  ""
  cost_us_op_tag <- if (cost_us_op) "_cost_us_op" else  ""
  cost_us_input_tag <- if(!is.na(cost_us_input)) paste0("_cost_us", cost_us_input) else ""
  
  us_sens_adj_tag <- ""
  if (us_sens_adj > 0 ) us_sens_adj_tag <- "_USsensUP" 
  if (us_sens_adj < 0 ) us_sens_adj_tag <- "_USsensDOWN" 
  
  us_spec_adj_tag <- ""
  if (us_spec_adj > 0 ) us_spec_adj_tag <- "_USspecUP" 
  if (us_spec_adj < 0 ) us_spec_adj_tag <- "_USspecDOWN" 
  
  # accuracy data
  
  IPD_i_u49 <- IPD_i %>% filter(age_cent60 < -1)
  IPD_i_o50 <- IPD_i %>% filter(age_cent60 >= -1)
  
  accuracy_tot <- readRDS(file.path(work_data, "psa", paste0("accuracy_psa_", accuracy_file, ".rds")))[[j]]
  
  accuracy_tot_u49 <- accuracy_tot[["u49"]]
  accuracy_tot_o50 <- accuracy_tot[["o50"]]
  
  # CA125, current
  u49 <- decision_tree(ca125_first=1, us_first=0, ovatools=0, dt=IPD_i_u49,
                       accuracy_tot=accuracy_tot_u49, 
                       cost_us_cdc=cost_us_cdc, cost_us_op=cost_us_op, cost_us_input=cost_us_input,
                       us_sens_adj=us_sens_adj, us_spec_adj=us_spec_adj)
  o50 <- decision_tree(ca125_first=1, us_first=0, ovatools=0, dt=IPD_i_o50,
                       accuracy_tot=accuracy_tot_o50, 
                       cost_us_cdc=cost_us_cdc, cost_us_op=cost_us_op, cost_us_input=cost_us_input,
                       us_sens_adj=us_sens_adj, us_spec_adj=us_spec_adj)
  x <- rbind(u49, o50) %>% select(-id)
  decitree100 <- IPD_ids %>% left_join(x)
  saveRDS(decitree100, file.path(output, "psa/dectr", "decitree100",
                                 paste0("decitree100_", accuracy_file, j, adj_pop_tag, 
                                        cost_us_cdc_tag, cost_us_op_tag, cost_us_input_tag,
                                        us_sens_adj_tag, us_spec_adj_tag, ".RDS")))

  # Ovatools sequential
  u49 <- decision_tree(ca125_first=1, us_first=0, ovatools=1, dt=IPD_i_u49,
                       accuracy_tot=accuracy_tot_u49, 
                       cost_us_cdc=cost_us_cdc, cost_us_op=cost_us_op, cost_us_input=cost_us_input,
                       us_sens_adj=us_sens_adj, us_spec_adj=us_spec_adj)
  o50 <- decision_tree(ca125_first=1, us_first=0, ovatools=1, dt=IPD_i_o50,
                       accuracy_tot=accuracy_tot_o50, 
                       cost_us_cdc=cost_us_cdc, cost_us_op=cost_us_op, cost_us_input=cost_us_input,
                       us_sens_adj=us_sens_adj, us_spec_adj=us_spec_adj)
  x <- rbind(u49, o50) %>% select(-id)
  decitree101 <- IPD_ids %>% left_join(x)
  saveRDS(decitree101, file.path(output, "psa/dectr", "decitree101",
                                 paste0("decitree101_", accuracy_file, j, adj_pop_tag, 
                                        cost_us_cdc_tag, cost_us_op_tag, cost_us_input_tag,
                                        us_sens_adj_tag, us_spec_adj_tag, ".RDS")))
}
  