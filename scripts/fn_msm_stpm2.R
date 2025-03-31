##################################################################
##                       Define functions                       ##
##################################################################

# From survival probability to transition probability ---------------------

surv_tp_cd <- function(surv_probs, j) {
  
  if (j==1) p_cd <- 1 - surv_probs[j] else
    
    p_cd <- ((1 - surv_probs[j]) - (1 - surv_probs[j-1]))/surv_probs[j-1]
  
  return(p_cd)
}

# Update by cycle ---------------------------------------------------------

gen_next_cycle <- function(
    j, vP_0, vP_1, b_qol, lam_b_i, cf_t, surv_probs
) { 
  
  ### time-updated characteristics ------------
  
  CurrAge <- vP_0["age_cent60"]
  vP_1["cycle"] <- j
  if(j==1) {
    vP_1["age_cent60"] <- CurrAge 
  } else {
    vP_1["age_cent60"] <- CurrAge + 1/10
  }
  
  # record flexpoint as parameter
  if (CurrAge < 1) {
    vP_0["age_cent60_1"] <- CurrAge
    vP_0["age_cent60_2"] <- 0
  } else {
    vP_0["age_cent60_1"] <- 1
    vP_0["age_cent60_2"] <- CurrAge - 1
  }
  
  ### predict cd, ncd -----------
  if (c_type!="benign" & (!use_mort | j<=8)) {
    p_cd <- surv_tp_cd(surv_probs, j) 
    
    p_ncd <- mort_prob[near(mort_prob[, "age_cent60"], CurrAge), c("ncd")]
    
    vP_1["ncd"] <- vP_0["alive"] * p_ncd
    
    vP_1["cd"] <- vP_0["alive"] * (1-p_ncd) * p_cd  
    
  } else{
    p_cd <- mort_prob[near(mort_prob[, "age_cent60"], CurrAge), c("cd")]
    
    p_ncd <- mort_prob[near(mort_prob[, "age_cent60"], CurrAge), c("ncd")]
    
    vP_1["ncd"] <- vP_0["alive"] * p_ncd
    
    vP_1["cd"] <- vP_0["alive"] * p_cd  
  }
  
  vP_1["alive"] <- max(vP_0["alive"] - vP_1["cd"] - vP_1["ncd"], 0)
  
  
  ### predict QoL -----------
  cf_t_qol <- cf_t[["qol"]]
  
  if (j==1) {
    # ids with baseline cancer have 0% chance to have the same cancer in 1 yr in decision tree
    # e.g. if c_type = "lung", ids with lung_2=1 can exist here,
    # but their result will multiply 0 from the decision tree result 
    # so we proceed with ids with lung_1=lung_2=1, but their result will not be included in further analysis
    if (c_type %in% c("lung", "loGI")) b_qol[paste0(c_type, "_1")] <- 1 else
      if (c_type %in% c("ova", "uter")) b_qol["gynae_1"] <- 1 else
        if (c_type == "panc") b_qol["upGI_1"] <- 1
  } else {
    b_qol["gynae_1"] <- b_qol["lung_1"] <- b_qol["upGI_1"] <- b_qol["loGI_1"] <- 0
    if (c_type %in% c("lung", "loGI")) b_qol[paste0(c_type, "_2")] <- 1 else
      if (c_type %in% c("ova", "uter")) b_qol["gynae_2"] <- 1 else
        if (c_type == "panc") b_qol["upGI_2"] <- 1
  }
  
  vP_0_names <- intersect(names(cf_t_qol), names(vP_0))
  b_qol_names <- intersect(names(cf_t_qol), names(b_qol))
  
  qol <- lam_b_i["qol"] + vP_0[vP_0_names] %*% cf_t_qol[vP_0_names] + 
    b_qol[b_qol_names] %*% cf_t_qol[b_qol_names]
  
  vP_1["qaly"] <- qol * vP_1["alive"]
  
  ### TODO predict cost -----------
  # the simulation start from diagnosis
  # cost start from diagnosis
  fcts <- rep(0 , 5)
  if (j >= 5){
    fcts[5] <- 1
  } else if (j >= 1 & j <=4){
    fcts[j] <- 1
  }
  
  vP_cost <- c(age_cent60 = unname(CurrAge),
               admi_year_fct1 = fcts[1],
               admi_year_fct2 = fcts[2],
               admi_year_fct3 = fcts[3],
               admi_year_fct4 = fcts[4],
               admi_year_fct5 = fcts[5],
               death_c = unname(vP_1["cd"]),
               death_nc = unname(vP_1["ncd"]),
               "stage_late:admi_year_fct1" = stage_late * fcts[1],
               "stage_late:admi_year_fct2" = stage_late * fcts[2],
               "stage_late:admi_year_fct3" = stage_late * fcts[3],
               "stage_late:admi_year_fct4" = stage_late * fcts[4],
               "stage_late:admi_year_fct5" = stage_late * fcts[5]
  )
  
  # Part one
  cost_p1_xb <- lam_b_i["cost_p1"] + vP_cost %*% cf_t[["cost_p1"]]
  cost_p1_prob <- 1/(1 + exp(-cost_p1_xb))
  # Part two
  # all cost model use log link
  cost_p2_cost <- exp(lam_b_i["cost_p2"] + vP_cost %*% cf_t[["cost_p2"]])
  
  vP_1["cost"] <- cost_p1_prob * cost_p2_cost * vP_1["alive"]
  
  vP_1["alive_disc"] <- vP_1["alive"]/(1 + 0.035)^(j-1)
  vP_1["qaly_disc"] <- vP_1["qaly"]/(1 + 0.035)^(j-1)
  vP_1["cost_disc"] <- vP_1["cost"]/(1 + 0.035)^(j-1)
  
  return(vP_1)
}


##################################################################
##                         Master Model                         ##
##################################################################

master <- function(
    IPD, mort_prob, age_stop, n_cores, output_folder, output_name, cf,
    c_type, stage_late, use_mort, save_cycle
){
  
  # load data ---------------------------------------------------------------
  
  IPD <- cbind(IPD, stage_late=stage_late, Intercept=1)    
  
  if(c_type=="benign") {
    # IPD[ ,"ben_gynae_1"] <- 1# include baseline and incident, based on estimate of he_ova
    # do not assume all non-cancer patients have benign gynae.
    
    IPD[ ,"stage_late"] <- 0 # ignore the initial input of stage_late
  }
  if(c_type=="other") {
    IPD[, "cancers_other_1"] <- 1 # for QoL calculation; it is not time-varing
  }
  
  stage_name <- if(stage_late==1) "late" else "early"
  
  cf_temp <- list(qol = cf[["qol"]], cost = cf[["cost"]])
  
  # extract the right CD model by c_type and stage -------------------------
    if (c_type == "benign") surv_probs_all <- NULL else
      surv_probs_all <- cf$surv_prob_list[[paste0(c_type, "_", stage_late)]]
  
  lam_b <- list()
  # QoL
  # consider baseline other cancer and benign gynae
  cf_b <- cf_temp[["qol"]][["cf_b"]]
  lam_b[["qol"]] <- IPD[, names(cf_b)] %*% cf_b
  
  cf_t <- list()
  cf_t[["qol"]] <- cf_temp[["qol"]][["cf_t"]]
  
  # cost
  cf_b_p1 <- cf_temp[["cost"]][[c_type]][["p1"]][["cf_b"]]
  cf_b_p2 <- cf_temp[["cost"]][[c_type]][["p2"]][["cf_b"]]
  lam_b[["cost_p1"]] <- IPD[, names(cf_b_p1)] %*% cf_b_p1
  lam_b[["cost_p2"]] <- IPD[, names(cf_b_p2)] %*% cf_b_p2
  
  cf_t[["cost_p1"]] <- cf_temp[["cost"]][[c_type]][["p1"]][["cf_t"]]
  cf_t[["cost_p2"]] <- cf_temp[["cost"]][[c_type]][["p2"]][["cf_t"]]
  
  # Start loop --------------------------------------------------------------
  
  ids <- 1:nrow(IPD) 
  
  cl <- makeSOCKcluster(n_cores)
  clusterExport(cl, c("surv_tp_cd", "gen_next_cycle", "IPD", "age_stop", 
                      "mort_prob", "c_type", "stage_late", "use_mort"
  ))
  clusterEvalQ(cl, library(tidyverse))
  registerDoSNOW(cl)
  
  retval <- parLapply(cl = cl, ids, function(i){
    
    # initial individual vector
    vP_0 <- c(IPD[i,]["id"], 
              cycle=0, 
              IPD[i,]["age_cent60"], 
              alive = 1,
              cd = 0, 
              ncd = 0,
              qaly = 0, 
              cost = 0, 
              alive_disc = 0,
              qaly_disc = 0,
              cost_disc = 0
    )
    
    vP_1 <- vP_0 # at t1
    
    b_qol <- c(IPD[i,]["gynae_2"],
               IPD[i,]["upGI_2"],
               IPD[i,]["lung_2"],
               IPD[i,]["loGI_2"],
               gynae_1=0,
               lung_1=0,
               upGI_1=0,
               loGI_1=0
    )
    
    J <- eval(age_stop)
    
    lam_b_i <- sapply(lam_b, "[[", i)
    
    # initialise output dataframe
    output_i <- matrix(nrow = J, ncol = length(vP_0))
    colnames(output_i) <- names(vP_0)
    
    # extract survival probability 
    if (c_type != "benign") {
      surv_probs <- surv_probs_all[i, ]
    } else surv_probs <- NULL
    
    ### update cycle -----
    for (j in 1:J) {
      
      alpha <- gen_next_cycle(
        j=j, vP_0=vP_0, vP_1=vP_1, b_qol=b_qol, lam_b_i=lam_b_i, cf_t=cf_t, 
        surv_probs = surv_probs
      )
      
      output_i[j, ] <- alpha
      
      vP_0 <- vP_1 <- alpha
      
    }
    
    return(output_i)
    
  })
  
  stopCluster(cl)
  
  
  # After loops -------------------------------------------------------------
  
  # cycle level output
  df_output <- do.call("rbind", retval)
  
  if (save_cycle) {
    write.csv(df_output, file.path(output_folder, paste0(output_name,".csv")), row.names = F)
  }
  
  # aggregate output
  df_output2 <- aggregate(df_output[, c("alive", "cd", "ncd", "qaly", "cost", 
                                        "alive_disc", "qaly_disc", "cost_disc")],
                          list(df_output[, "id"]), sum)
  # force the first columne named id
  colnames(df_output2)[1] <- "id"
  
  IPD <- cbind(IPD, age = IPD[, "age_cent60"]*10+60)
  
  # df_output2 <- cbind(IPD[, c("id", "age")], 
  #                     df_output2[, c("alive", "cd", "ncd", "qaly", "cost")]
  # )
  
  # merge by id is much safer!
  df_output2 <- merge(as.data.frame(IPD[, c("id", "age")]),df_output2, by = "id")
  
  write.csv(df_output2, file.path(output_folder, paste0(output_name,"_sum.csv")), row.names = F)
  return(df_output)
  
}









