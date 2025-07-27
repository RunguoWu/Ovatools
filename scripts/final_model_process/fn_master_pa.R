#################################################################
##              Master function for post-analysis              ##
#################################################################

# adapte from the same name file in post_analysis folder

fn_pa_master <- function(i=1001, dect0="decitree100", dect1="decitree101",  
                         accuracy_file, adj_pop=FALSE, 
                         cost_us_cdc=FALSE, cost_us_op=FALSE, cost_us_input=NA,
                         disuti = TRUE,
                         uti_ben = TRUE,
                         extra_cost = TRUE,
                         complication = FALSE,
                         twoStepUS = FALSE,
                         Markov_snr = "", 
                         surg_prob_adj = 0, 
                         surg_disu_adj = 0,
                         shift_late_sens = "", # "upper" or "lower"
                         shift_early_sens="",
                         us_sens_adj_tag="", us_spec_adj_tag="" # "_USsensUP", "_USsensDOWN", "_USspecUP", "_USspecDOWN"
# Markov_snr = "NostageQoL", "stageQoL006", "cd_prob15yr", "disc015"
                         ){
  # this function keep more details
  # the following post_ana4psa only keep essentials
  adj_pop_tag <- if (adj_pop) "_adj_pop" else ""
  
  disuti_tag <- if (!disuti) "_NOdisuti" else ""
  
  extra_cost_tag <- if (!extra_cost) "_NOextra_cost" else ""
  
  complication_tag <- if (complication) "_comp" else ""
  
  uti_ben_tag <- if (!uti_ben) "_NOuti_ben" else ""
  
  cost_us_cdc_tag <- if (cost_us_cdc) "_cost_us_cdc" else  ""
  cost_us_op_tag <- if (cost_us_op) "_cost_us_op" else  ""
  cost_us_input_tag <- if(!is.na(cost_us_input)) paste0("_cost_us", cost_us_input) else ""
  
  surg_prob_adj_tag <- ""
  if (surg_prob_adj > 0 ) surg_prob_adj_tag <- "_surg_probUP" 
  if (surg_prob_adj < 0 ) surg_prob_adj_tag <- "_surg_probDOWN" 
  
  surg_disu_adj_tag <- if (surg_disu_adj !=0) paste0("surg_disu_adj", surg_disu_adj) else ""
  
  vs_name <- gsub("decitree", "", (paste0("pa", dect0, "v", dect1)))
  
  ca125varying_name <- if (twoStepUS) "CA125varyingThreshold_2stepUS" else "CA125varyingThreshold"
  
  accuracy_file0 <- if (dect0 == "decitree112var" | dect0 == "decitree101var") ca125varying_name else accuracy_file
  
  decitree0 <- readRDS(file.path(
    output,"psa/dectr", dect0,
    paste0(dect0, "_", accuracy_file0, i, adj_pop_tag, cost_us_cdc_tag, cost_us_op_tag, cost_us_input_tag,
           us_sens_adj_tag, us_spec_adj_tag, ".RDS")))
  
  accuracy_file1 <- if (dect1 == "decitree112var" | dect1 == "decitree101var") ca125varying_name else accuracy_file
  
  decitree1 <- readRDS(file.path(
    output,"psa/dectr", dect1,
    paste0(dect1, "_", accuracy_file1, i, adj_pop_tag, cost_us_cdc_tag, cost_us_op_tag, cost_us_input_tag,
           us_sens_adj_tag, us_spec_adj_tag, ".RDS")))
  
  stage_shift <- stage_shift_psa[[i]]
  FP_data <- FP_data_psa[[i]]
  
  # map them to the whole population
  rst_list <- c("decitree0", "decitree1")
  for (rst in rst_list) {
    to_update <- get(rst)
    updated <- map2pop(to_update, id_dic3, ipd2)
    updated <- updated %>% select(-group)
    assign(rst, updated)
  }

  # import Markov results
  if (Markov_snr==""){
    benign <- read.csv(file.path(output, "psa", "benign_0", paste0("psa_benign_stage0_useMort_TRUE_", i, "_sum.csv")))
    ova_early <- read.csv(file.path(output, "psa", "ova_0", paste0("psa_ova_stage0_useMort_TRUE_", i, "_sum.csv")))
    ova_late <- read.csv(file.path(output, "psa", "ova_1", paste0("psa_ova_stage1_useMort_TRUE_", i, "_sum.csv")))
    lung_early <- read.csv(file.path(output, "psa", "lung_0", paste0("psa_lung_stage0_useMort_TRUE_", i, "_sum.csv")))
    lung_late <- read.csv(file.path(output, "psa", "lung_1", paste0("psa_lung_stage1_useMort_TRUE_", i, "_sum.csv")))
    panc_early <- read.csv(file.path(output, "psa", "panc_0", paste0("psa_panc_stage0_useMort_TRUE_", i, "_sum.csv")))
    panc_late <- read.csv(file.path(output, "psa", "panc_1", paste0("psa_panc_stage1_useMort_TRUE_", i, "_sum.csv")))
    uter_early <- read.csv(file.path(output, "psa", "uter_0", paste0("psa_uter_stage0_useMort_TRUE_", i, "_sum.csv")))
    uter_late <- read.csv(file.path(output, "psa", "uter_1", paste0("psa_uter_stage1_useMort_TRUE_", i, "_sum.csv")))
    loGI_early <- read.csv(file.path(output, "psa", "loGI_0", paste0("psa_loGI_stage0_useMort_TRUE_", i, "_sum.csv")))
    loGI_late <- read.csv(file.path(output, "psa", "loGI_1", paste0("psa_loGI_stage1_useMort_TRUE_", i, "_sum.csv")))
    other_early <- read.csv(file.path(output, "psa", "other_0", paste0("psa_other_stage0_useMort_TRUE_", i, "_sum.csv")))
    other_late <- read.csv(file.path(output, "psa", "other_1", paste0("psa_other_stage1_useMort_TRUE_", i, "_sum.csv")))
  } else {
    # sensitivity analysis
    benign <- read.csv(file.path(output, paste0(Markov_snr, "_benign_stage0_useMort_TRUE_sum.csv")))
    ova_early <- read.csv(file.path(output, paste0(Markov_snr, "_ova_stage0_useMort_TRUE_sum.csv")))
    ova_late <- read.csv(file.path(output, paste0(Markov_snr, "_ova_stage1_useMort_TRUE_sum.csv")))
    lung_early <- read.csv(file.path(output, paste0(Markov_snr, "_lung_stage0_useMort_TRUE_sum.csv")))
    lung_late <- read.csv(file.path(output, paste0(Markov_snr, "_lung_stage1_useMort_TRUE_sum.csv")))
    panc_early <- read.csv(file.path(output, paste0(Markov_snr, "_panc_stage0_useMort_TRUE_sum.csv")))
    panc_late <- read.csv(file.path(output, paste0(Markov_snr, "_panc_stage1_useMort_TRUE_sum.csv")))
    uter_early <- read.csv(file.path(output, paste0(Markov_snr, "_uter_stage0_useMort_TRUE_sum.csv")))
    uter_late <- read.csv(file.path(output, paste0(Markov_snr, "_uter_stage1_useMort_TRUE_sum.csv")))
    loGI_early <- read.csv(file.path(output, paste0(Markov_snr, "_loGI_stage0_useMort_TRUE_sum.csv")))
    loGI_late <- read.csv(file.path(output, paste0(Markov_snr, "_loGI_stage1_useMort_TRUE_sum.csv")))
    other_early <- read.csv(file.path(output, paste0(Markov_snr, "_other_stage0_useMort_TRUE_sum.csv")))
    other_late <- read.csv(file.path(output, paste0(Markov_snr, "_other_stage1_useMort_TRUE_sum.csv")))
    Markov_snr <- paste0("_", Markov_snr)
  }
    
  # map them to the whole population
  rst_list <- c("benign", "ova_early", "ova_late", "lung_early", "lung_late", 
                "panc_early", "panc_late", "uter_early", "uter_late", 
                "loGI_early", "loGI_late", "other_early", "other_late")
  
  for (rst in rst_list) {
    to_update <- get(rst)
    updated <- map2pop(to_update, id_dic3, ipd2)
    updated <- updated %>% select(-group, -epatid)
    assign(rst, updated)
  }
  
  comp <- stage_shift_effect(decitree0, decitree1, with_real_dist = F, stage_shift,
                                    shift_other = T, 
                             shift_late_sens = shift_late_sens, 
                             shift_early_sens = shift_early_sens)
  
  if (shift_late_sens !="") shift_late_sens <- paste0("_shift_late_", shift_late_sens)
  if (shift_early_sens !="") shift_early_sens <- paste0("_shift_early_", shift_early_sens)
  
  dt <- link_markov(he_ova = he_ova, comp = comp,
                    benign = benign,
                    ova_early = ova_early, ova_late = ova_late,
                    lung_early = lung_early, lung_late = lung_late,
                    panc_early = panc_early, panc_late = panc_late,
                    uter_early = uter_early, uter_late = uter_late,
                    loGI_early = loGI_early, loGI_late = loGI_late,
                    other_early = other_early, other_late = other_late,
                    FP_data = FP_data,
                    disuti = disuti,
                    uti_ben = uti_ben,
                    extra_cost = extra_cost,
                    complication = complication,
                    surg_prob_adj = surg_prob_adj, 
                    surg_disu_adj = surg_disu_adj
  )

  # only look at those with CA125 records
  if (adj_pop) dt <- dt %>% filter(!is.na(CA125))
  
  rt1 <- ana_rslt(dt) %>% bind_cols(rownames(.), .)
  rt2 <- ana_rslt(dt, group="age_50") %>% bind_cols(rownames(.), .)
  rt3 <- ana_rslt(dt, age50filter=">=50", group="ethn4") %>% bind_cols(rownames(.), .)
  rt4 <- ana_rslt(dt, age50filter=">=50", group="town5") %>% bind_cols(rownames(.), .)
  
  rt <- bind_rows(rt1, rt2,rt3, rt4) %>%
    relocate("get(group)") %>%
    rename(group = "get(group)") %>% 
    mutate(group = if_else(is.na(group), "Total", group)) %>% 
    select(-"...1") %>%
    select("group", "N", "icer", "icer_uo", "icer_ly", "icer_ly_uo", "TP_pct_0", "TP12_pct_0", # "sum_ova_prev","sum_stage_late",
           "TP34_pct_0", "TP_pct_1", "TP12_pct_1", "TP34_pct_1", "ref_pct_0",
           "ref_pct_1", "FP_pct_0", "FP_pct_1", "stage_change_pct",
           "sum_qaly_disc_0", "sum_qaly_disc_1", "sum_qaly_disc_uo_1",
           "sum_cost_disc_tot_0", "sum_cost_disc_tot_1", "sum_cost_disc_tot_uo_1",
           "sum_ly_dif", "sum_ly_uo_dif", "sum_ly_disc_dif", "sum_ly_disc_uo_dif",
           "sum_qaly_dif", "sum_qaly_uo_dif", "sum_qaly_disc_dif", "sum_qaly_disc_uo_dif",
           "sum_cost_tot_dif","sum_cost_tot_uo_dif", "sum_cost_disc_tot_dif", "sum_cost_disc_tot_uo_dif",
           "sum_cost_dif", "sum_cost_uo_dif", "sum_cost_pc_dif", "sum_cost_dif_fromFP"
    )
  write.csv(rt, file.path(output, "psa", "CEA", vs_name,
                          paste0(vs_name, "_", i, "_utiben_", accuracy_file,
                                 adj_pop_tag, Markov_snr,
                                 cost_us_cdc_tag, cost_us_op_tag, cost_us_input_tag, disuti_tag,
                                 uti_ben_tag, extra_cost_tag, complication_tag, shift_late_sens,
                                 shift_early_sens, us_sens_adj_tag, us_spec_adj_tag,
                                 surg_prob_adj_tag, surg_disu_adj_tag, ".csv")))
  
}



psa_summary <- function(versus = "100v101", utiben = FALSE, accuracy_file, adj_pop=FALSE){
  
  disu <- if (utiben) "_utiben" else ""
  adj_pop_tag <- if (adj_pop) "_adj_pop" else ""
  # deterministic analysis
  i <- 1001
  
  folder_name <- paste0("pa", versus)
  
  file_name <- paste0(folder_name, "_", i, disu, "_", accuracy_file, adj_pop_tag, ".csv")
  
  det <- read.csv(file.path(output, "psa", "CEA", folder_name, file_name))
  
  det <- det %>%
    mutate(ly_disc_dif_per1000 = sum_ly_disc_dif/N * 1000,
           qaly_disc_dif_per1000 = sum_qaly_disc_dif/N * 1000,
           ly_disc_uo_dif_per1000 = sum_ly_disc_uo_dif/N * 1000,
           qaly_disc_uo_dif_per1000 = sum_qaly_disc_uo_dif/N * 1000,
           cost_disc_tot_dif_per1000 = sum_cost_disc_tot_dif/N * 1000,
           cost_disc_tot_uo_dif_per1000 = sum_cost_disc_tot_uo_dif/N * 1000
    )
  
  col_names <- colnames(det)[(-1:-3)]
  for(nam in col_names){
    x <- det[, nam]
    assign(nam, x)
  }
  
  groups <- det$group
  
  # PSA
  for (j in 1:1000) {
    
    file_name <- paste0(folder_name, "_", j, disu, "_", accuracy_file, adj_pop_tag, ".csv")
    psa <- read.csv(file.path(output, "psa", "CEA",folder_name, file_name))
    
    psa <- psa %>%
    mutate(ly_disc_dif_per1000 = sum_ly_disc_dif/N * 1000,
           qaly_disc_dif_per1000 = sum_qaly_disc_dif/N * 1000,
           ly_disc_uo_dif_per1000 = sum_ly_disc_uo_dif/N * 1000,
           qaly_disc_uo_dif_per1000 = sum_qaly_disc_uo_dif/N * 1000,
           cost_disc_tot_dif_per1000 = sum_cost_disc_tot_dif/N * 1000,
           cost_disc_tot_uo_dif_per1000 = sum_cost_disc_tot_uo_dif/N * 1000
    )
    
    for(nam in col_names){
      x <- psa[, nam]
      assign(nam, cbind(get(nam), x))
    }
    
  }
  
  # summary
  rt_all <- c()
  
  for(nam in col_names){
    det <- get(nam)[, 1]
    psa <- get(nam)[,-1]
    
    prob <- rowQuantiles(psa, probs = c(0.025, 0.975))
    
    rt <- cbind(rep(nam, length(det)), groups, det, prob)
    # rt <- rbind(c(nam, NA, NA), rt)
    
    rt_all <- rbind(rt_all, rt)
  }
  
  return(rt_all)
}


# to extract detection rate, 

output1 <- function(x, col, decimal=2, type="pct"){ 
  # type = pct for percentage; NA for LY and QALY C for cost
  
  scale_up <- 1
  scale_down <- 1
  pct_symb <- ""
  if (type == "pct") {
    scale_up <- 100
    # pct_symb <- "%"
  } else if (type == "cost"){
    scale_down <- 1 # 1/1000000 # per million
  }
  
  x <- as.data.frame(x)
  
  groups <- x[x$V1==col, ]$groups
  
  det <- round(as.numeric(x[x$V1==col,]$det)*scale_down, decimal)*scale_up
  lb <- round(as.numeric(x[x$V1==col,]$`2.5%`)*scale_down, decimal)*scale_up
  ub <- round(as.numeric(x[x$V1==col,]$`97.5%`)*scale_down, decimal)*scale_up
  
  rt <- cbind(groups, paste0(det, pct_symb, " (", lb, " - ", ub, ")"))
  
  return(rt)
}

# to quickly load scenarios

snr_load <- function(snr_tag, adj_pop_tag = "", accuracy_file="byAgeOnly_niceUS"){
  
  # snr_tag <- "_NOuti_ben"
  rt <- c()
  name_x <- paste0("icer_101", snr_tag)
  x <- read.csv(
    file.path(output, "psa", "CEA", "pa100v101", 
              paste0("pa100v101_1001_utiben_", accuracy_file, adj_pop_tag, snr_tag, ".csv"))
  ) %>% filter(group==">=50") %>% select(icer)%>% as.numeric() %>% round()
  rt[name_x]=x

  # name_x <- paste0("icer_112", snr_tag)
  # x <- read.csv(
  #   file.path(output, "psa", "CEA", "pa100v112",  
  #             paste0("pa100v112_1001_utiben_", accuracy_file, adj_pop_tag, snr_tag, ".csv"))
  # ) %>% filter(group==">=50") %>% select(icer)%>% as.numeric() %>% round()
  # rt[name_x]=x
  # 
  # name_x <- paste0("icer_101var", snr_tag)
  # x <- read.csv(
  #   file.path(output, "psa", "CEA", "pa100v101var", 
  #             paste0("pa100v101var_1001_utiben_", accuracy_file, adj_pop_tag, snr_tag, ".csv"))
  # ) %>% filter(group==">=50") %>% select(icer)%>% as.numeric() %>% round()
  # rt[name_x]=x
  # 
  # name_x <- paste0("icer_110", snr_tag)
  # x <- read.csv(
  #   file.path(output, "psa", "CEA", "pa100v110",  
  #             paste0("pa100v110_1001_utiben_", accuracy_file, adj_pop_tag, snr_tag, ".csv"))
  # ) %>% filter(group==">=50") %>% select(icer)%>% as.numeric() %>% round()
  # rt[name_x]=x
  
  name_x <- paste0("icer_110", snr_tag)
  x <- read.csv(
    file.path(output, "psa", "CEA", "pa101v110",  
              paste0("pa101v110_1001_utiben_", accuracy_file, adj_pop_tag, snr_tag, ".csv"))
  ) %>% filter(group==">=50") %>% select(icer)%>% as.numeric() %>% round()
  rt[name_x]=x
  
  return(rt)  
}

# For ICER test only

fn_pa_master2 <- function(i=1001, j, dect0="decitree100", dect1="decitree101",  
                         accuracy_file, adj_pop=FALSE, 
                         cost_us_cdc=FALSE, cost_us_op=FALSE, cost_us_input=NA,
                         disuti = TRUE,
                         uti_ben = TRUE,
                         extra_cost = TRUE,
                         twoStepUS = FALSE,
                         Markov_snr = "", 
                         surg_prob_adj = 0, 
                         shift_late_sens = "", # "upper" or "lower"
                         shift_early_sens="",
                         us_sens_adj_tag="", us_spec_adj_tag="" # "_USsensUP", "_USsensDOWN", "_USspecUP", "_USspecDOWN"
                         # Markov_snr = "NostageQoL", "stageQoL006", "cd_prob15yr", "disc015"
){
  # this function keep more details
  # the following post_ana4psa only keep essentials
  adj_pop_tag <- if (adj_pop) "_adj_pop" else ""
  
  disuti_tag <- if (!disuti) "_NOdisuti" else ""
  
  extra_cost_tag <- if (!extra_cost) "_NOextra_cost" else ""
  
  uti_ben_tag <- if (!uti_ben) "_NOuti_ben" else ""
  
  cost_us_cdc_tag <- if (cost_us_cdc) "_cost_us_cdc" else  ""
  cost_us_op_tag <- if (cost_us_op) "_cost_us_op" else  ""
  cost_us_input_tag <- if(!is.na(cost_us_input)) paste0("_cost_us", cost_us_input) else ""
  
  surg_prob_adj_tag <- ""
  if (surg_prob_adj > 0 ) surg_prob_adj_tag <- "_surg_probUP" 
  if (surg_prob_adj < 0 ) surg_prob_adj_tag <- "_surg_probDOWN" 
  
  vs_name <- gsub("decitree", "", (paste0("pa", dect0, "v", dect1)))
  
  ca125varying_name <- if (twoStepUS) "CA125varyingThreshold_2stepUS" else "CA125varyingThreshold"
  
  accuracy_file0 <- if (dect0 == "decitree100var" | dect0 == "decitree101var") ca125varying_name else accuracy_file
  
  decitree0 <- readRDS(file.path(
    output,"psa/dectr", dect0,
    paste0(dect0, "_", accuracy_file0, j, adj_pop_tag, cost_us_cdc_tag, cost_us_op_tag, cost_us_input_tag,
           us_sens_adj_tag, us_spec_adj_tag, ".RDS")))
  
  accuracy_file1 <- if (dect1 == "decitree100var" | dect1 == "decitree101var") ca125varying_name else accuracy_file
  
  decitree1 <- readRDS(file.path(
    output,"psa/dectr", dect1,
    paste0(dect1, "_", accuracy_file1, j, adj_pop_tag, cost_us_cdc_tag, cost_us_op_tag, cost_us_input_tag,
           us_sens_adj_tag, us_spec_adj_tag, ".RDS")))
  
  stage_shift <- stage_shift_psa[[i]]
  FP_data <- FP_data_psa[[i]]
  
  # map them to the whole population
  rst_list <- c("decitree0", "decitree1")
  for (rst in rst_list) {
    to_update <- get(rst)
    updated <- map2pop(to_update, id_dic3, ipd2)
    updated <- updated %>% select(-group)
    assign(rst, updated)
  }
  
  # import Markov results
  if (Markov_snr==""){
    benign <- read.csv(file.path(output, "psa", "benign_0", paste0("psa_benign_stage0_useMort_TRUE_", i, "_sum.csv")))
    ova_early <- read.csv(file.path(output, "psa", "ova_0", paste0("psa_ova_stage0_useMort_TRUE_", i, "_sum.csv")))
    ova_late <- read.csv(file.path(output, "psa", "ova_1", paste0("psa_ova_stage1_useMort_TRUE_", i, "_sum.csv")))
    lung_early <- read.csv(file.path(output, "psa", "lung_0", paste0("psa_lung_stage0_useMort_TRUE_", i, "_sum.csv")))
    lung_late <- read.csv(file.path(output, "psa", "lung_1", paste0("psa_lung_stage1_useMort_TRUE_", i, "_sum.csv")))
    panc_early <- read.csv(file.path(output, "psa", "panc_0", paste0("psa_panc_stage0_useMort_TRUE_", i, "_sum.csv")))
    panc_late <- read.csv(file.path(output, "psa", "panc_1", paste0("psa_panc_stage1_useMort_TRUE_", i, "_sum.csv")))
    uter_early <- read.csv(file.path(output, "psa", "uter_0", paste0("psa_uter_stage0_useMort_TRUE_", i, "_sum.csv")))
    uter_late <- read.csv(file.path(output, "psa", "uter_1", paste0("psa_uter_stage1_useMort_TRUE_", i, "_sum.csv")))
    loGI_early <- read.csv(file.path(output, "psa", "loGI_0", paste0("psa_loGI_stage0_useMort_TRUE_", i, "_sum.csv")))
    loGI_late <- read.csv(file.path(output, "psa", "loGI_1", paste0("psa_loGI_stage1_useMort_TRUE_", i, "_sum.csv")))
    other_early <- read.csv(file.path(output, "psa", "other_0", paste0("psa_other_stage0_useMort_TRUE_", i, "_sum.csv")))
    other_late <- read.csv(file.path(output, "psa", "other_1", paste0("psa_other_stage1_useMort_TRUE_", i, "_sum.csv")))
  } else {
    # sensitivity analysis
    benign <- read.csv(file.path(output, paste0(Markov_snr, "_benign_stage0_useMort_TRUE_sum.csv")))
    ova_early <- read.csv(file.path(output, paste0(Markov_snr, "_ova_stage0_useMort_TRUE_sum.csv")))
    ova_late <- read.csv(file.path(output, paste0(Markov_snr, "_ova_stage1_useMort_TRUE_sum.csv")))
    lung_early <- read.csv(file.path(output, paste0(Markov_snr, "_lung_stage0_useMort_TRUE_sum.csv")))
    lung_late <- read.csv(file.path(output, paste0(Markov_snr, "_lung_stage1_useMort_TRUE_sum.csv")))
    panc_early <- read.csv(file.path(output, paste0(Markov_snr, "_panc_stage0_useMort_TRUE_sum.csv")))
    panc_late <- read.csv(file.path(output, paste0(Markov_snr, "_panc_stage1_useMort_TRUE_sum.csv")))
    uter_early <- read.csv(file.path(output, paste0(Markov_snr, "_uter_stage0_useMort_TRUE_sum.csv")))
    uter_late <- read.csv(file.path(output, paste0(Markov_snr, "_uter_stage1_useMort_TRUE_sum.csv")))
    loGI_early <- read.csv(file.path(output, paste0(Markov_snr, "_loGI_stage0_useMort_TRUE_sum.csv")))
    loGI_late <- read.csv(file.path(output, paste0(Markov_snr, "_loGI_stage1_useMort_TRUE_sum.csv")))
    other_early <- read.csv(file.path(output, paste0(Markov_snr, "_other_stage0_useMort_TRUE_sum.csv")))
    other_late <- read.csv(file.path(output, paste0(Markov_snr, "_other_stage1_useMort_TRUE_sum.csv")))
    Markov_snr <- paste0("_", Markov_snr)
  }
  
  # map them to the whole population
  rst_list <- c("benign", "ova_early", "ova_late", "lung_early", "lung_late", 
                "panc_early", "panc_late", "uter_early", "uter_late", 
                "loGI_early", "loGI_late", "other_early", "other_late")
  
  for (rst in rst_list) {
    to_update <- get(rst)
    updated <- map2pop(to_update, id_dic3, ipd2)
    updated <- updated %>% select(-group, -epatid)
    assign(rst, updated)
  }
  
  comp <- stage_shift_effect(decitree0, decitree1, with_real_dist = F, stage_shift,
                             shift_other = T, 
                             shift_late_sens = shift_late_sens, 
                             shift_early_sens = shift_early_sens)
  
  if (shift_late_sens !="") shift_late_sens <- paste0("_shift_late_", shift_late_sens)
  if (shift_early_sens !="") shift_early_sens <- paste0("_shift_early_", shift_early_sens)
  
  dt <- link_markov(he_ova = he_ova, comp = comp,
                    benign = benign,
                    ova_early = ova_early, ova_late = ova_late,
                    lung_early = lung_early, lung_late = lung_late,
                    panc_early = panc_early, panc_late = panc_late,
                    uter_early = uter_early, uter_late = uter_late,
                    loGI_early = loGI_early, loGI_late = loGI_late,
                    other_early = other_early, other_late = other_late,
                    FP_data = FP_data,
                    disuti = disuti,
                    uti_ben = uti_ben,
                    extra_cost = extra_cost,
                    surg_prob_adj = surg_prob_adj
  )
  
  # only look at those with CA125 records
  if (adj_pop) dt <- dt %>% filter(!is.na(CA125))
  
  rt1 <- ana_rslt(dt) %>% bind_cols(rownames(.), .)
  rt2 <- ana_rslt(dt, group="age_50") %>% bind_cols(rownames(.), .)
  rt3 <- ana_rslt(dt, age50filter=">=50", group="ethn4") %>% bind_cols(rownames(.), .)
  rt4 <- ana_rslt(dt, age50filter=">=50", group="town5") %>% bind_cols(rownames(.), .)
  
  rt <- bind_rows(rt1, rt2,rt3, rt4) %>%
    relocate("get(group)") %>%
    rename(group = "get(group)") %>% 
    mutate(group = if_else(is.na(group), "Total", group)) %>% 
    select(-"...1") %>%
    select("group", "N", "icer", "icer_uo", "icer_ly", "icer_ly_uo", "TP_pct_0", "TP12_pct_0", # "sum_ova_prev","sum_stage_late",
           "TP34_pct_0", "TP_pct_1", "TP12_pct_1", "TP34_pct_1", "ref_pct_0",
           "ref_pct_1", "FP_pct_0", "FP_pct_1", "stage_change_pct",
           "sum_qaly_disc_0", "sum_qaly_disc_1", "sum_qaly_disc_uo_1",
           "sum_cost_disc_tot_0", "sum_cost_disc_tot_1", "sum_cost_disc_tot_uo_1",
           "sum_ly_dif", "sum_ly_uo_dif", "sum_ly_disc_dif", "sum_ly_disc_uo_dif",
           "sum_qaly_dif", "sum_qaly_uo_dif", "sum_qaly_disc_dif", "sum_qaly_disc_uo_dif",
           "sum_cost_tot_dif","sum_cost_tot_uo_dif", "sum_cost_disc_tot_dif", "sum_cost_disc_tot_uo_dif",
           "sum_cost_dif", "sum_cost_uo_dif", "sum_cost_pc_dif", "sum_cost_dif_fromFP"
    )
  write.csv(rt, file.path(output, "psa", "CEA", vs_name, 
                          paste0(vs_name, "_", j, "_utiben_", accuracy_file, 
                                 adj_pop_tag, Markov_snr, 
                                 cost_us_cdc_tag, cost_us_op_tag, cost_us_input_tag, disuti_tag,
                                 uti_ben_tag, extra_cost_tag, shift_late_sens, 
                                 shift_early_sens, us_sens_adj_tag, us_spec_adj_tag,
                                 surg_prob_adj_tag, ".csv")))
  
}
