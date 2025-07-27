##################################################################
##                 Master for the decision tree                 ##
##                            for PSA                           ##
##################################################################

# adapted from master_decision_tree_psa.R in HE_model folder
rm(list = ls())

library(tidyverse)
library(haven)
library(foreign)
library(dtplyr)
library(data.table)
library(foreach)
library(snowfall)
library(doSNOW)

source("//Qmcr/dfs/WIPH/Primary Care/OVATOOLS/r_project/scripts/file_path.R")
source(file.path(scripts,"final_model_process", "/fn_decision_tree_split_stage.R"))
source(file.path(scripts,"final_model_process", "fn_master_decision_tree.R"))

# import ipd data
IPD <- readRDS(file.path(work_data, "ipd2_20241126.RDS")) %>% 
  mutate(id=1:n()) %>% 
  relocate(id, .before = age_cent60)
id_dic2 <- readRDS(file.path(work_data, "id_dic2_20241126.RDS"))

IPD_ids <- IPD %>% select(id, group)

N=1001
ncores <- 20

sfInit(cpus = ncores, parallel = T)
cl <- sfGetCluster()
# clusterExport(cl=cl, list = list("fn_x", "fn_ptr", "fn_npv", "fn_ppv",
#                                  "fn_prevalence_fl_ca125", "fn_prevalence_fl_ovt",
#                                  "fn_prevalence_fl_us"))
clusterExport(cl=cl, list = list())
clusterEvalQ(cl, library(tidyverse))
clusterEvalQ(cl, library(dtplyr))
clusterEvalQ(cl, library(data.table))

registerDoSNOW(cl)
ptm <- proc.time()
# i=1001

foreach(i = 1:N)%dopar%{
  
  master_decision_tree(i, accuracy_file = "byAgeOnly_niceUS")
  
  master_decision_tree(i, accuracy_file = "byAgeOnly_niceUS", adj_pop = TRUE)
  
}

sfStop()

print(proc.time() - ptm)
print(Sys.time())


# scenario analyses -------------------------------------------------------
adj_pop = FALSE # in the study population or in the CA125 record population

#> Varying GP US cost----
# Reduce US cost to CDC cost
master_decision_tree(accuracy_file = "byAgeOnly_niceUS", adj_pop=adj_pop, cost_us_cdc=TRUE)
# increase US cost to outpatient cost
master_decision_tree(accuracy_file = "byAgeOnly_niceUS", adj_pop=adj_pop, cost_us_op=TRUE)

#> Using by age and stage accuracy for CA125 and Ovatools----
master_decision_tree(accuracy_file = "byAgeStage_niceUS", adj_pop=adj_pop)

#> Varying USS accuracy----
# sensitivity plus 0.05
master_decision_tree(accuracy_file = "byAgeOnly_niceUS", adj_pop=adj_pop, us_sens_adj=0.05)
# sensitivity minus 0.05
master_decision_tree(accuracy_file = "byAgeOnly_niceUS", adj_pop=adj_pop, us_sens_adj=-0.05)
# specificity plus 0.05
master_decision_tree(accuracy_file = "byAgeOnly_niceUS", adj_pop=adj_pop, us_spec_adj=0.05)
# specificity minus 0.05
master_decision_tree(accuracy_file = "byAgeOnly_niceUS", adj_pop=adj_pop, us_spec_adj=-0.05)











