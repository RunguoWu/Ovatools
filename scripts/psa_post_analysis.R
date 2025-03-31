##################################################################################
##  Combining simulation data from decision tree, Markov model and other steps  ##
##################################################################################

# adapte from master_post_analysis_psa.R in post_analysis folder

rm(list = ls())

library(tidyverse)
library(foreach)
library(doSNOW)
library(snowfall)
source("//Qmcr/dfs/WIPH/Primary Care/OVATOOLS/r_project/scripts/file_path.R")
source(file.path(scripts, "final_model_process", "fn_post_analysis.R"))
source(file.path(scripts, "final_model_process", "fn_master_pa.R"))

he_ova <- readRDS(file.path(work_data, "he_ova_step4_withSymptom_prescr_NCRASupdated.RDS"))
stage_shift_psa <- readRDS(file.path(work_data, "psa", "stage_shift_psa.rds"))
FP_data_psa <- readRDS(file.path(work_data, "psa", "FP_data_psa.rds"))

# load in full population and unique characteristics connectors
ipd2 <- readRDS(file.path(work_data, "ipd2_20241126.rds"))
# ipd2_20241126.rds and ipd2_20240906.RDS only have different var "group"
# as the size and order are not changed, Markov model results are unchanged from ipd2_20240906.RDS are unchanged
# because Markov model results does not include variable "group"
id_dic3 <- readRDS(file.path(work_data, "id_dic3_20241126.rds"))

ncores = 20
sfInit(cpus = ncores, parallel = T)
cl <- sfGetCluster()
clusterExport(cl=cl, list=list())
clusterEvalQ(cl, library(tidyverse))
registerDoSNOW(cl)

foreach (i = 1:1000) %dopar%{

  fn_pa_master(i, dect0="decitree100", dect1="decitree101",
               accuracy_file= "byAgeOnly_niceUS")
  fn_pa_master(i, dect0="decitree100", dect1="decitree110",
               accuracy_file= "byAgeOnly_niceUS")
  fn_pa_master(i, dect0="decitree100", dect1="decitree112",
               accuracy_file= "byAgeOnly_niceUS")
  fn_pa_master(i, dect0="decitree100", dect1="decitree101var",
               accuracy_file= "byAgeOnly_niceUS")
  fn_pa_master(i, dect0="decitree100", dect1="decitree112var",
               accuracy_file= "byAgeOnly_niceUS")

  # fn_pa_master(i, dect0="decitree100", dect1="decitree101",
  #              accuracy_file= "byAgeOnly_niceUS", adj_pop=TRUE)
  # fn_pa_master(i, dect0="decitree100", dect1="decitree110",
  #              accuracy_file= "byAgeOnly_niceUS", adj_pop=TRUE)
  # fn_pa_master(i, dect0="decitree100", dect1="decitree112",
  #              accuracy_file= "byAgeOnly_niceUS", adj_pop=TRUE)
  # fn_pa_master(i, dect0="decitree100", dect1="decitree101var",
  #              accuracy_file= "byAgeOnly_niceUS", adj_pop=TRUE)
  # fn_pa_master(i, dect0="decitree100", dect1="decitree112var",
  #              accuracy_file= "byAgeOnly_niceUS", adj_pop=TRUE)

}

sfStop()


# Scenario Analyses -------------------------------------------------------
adj_pop <- FALSE # in the study population or in the CA125 record population

#> base case----
# fn_pa_master(dect0="decitree100", dect1="decitree101",
#              accuracy_file= "byAgeOnly_niceUS", adj_pop=adj_pop)

fn_pa_master(dect0="decitree101", dect1="decitree112",
             accuracy_file= "byAgeOnly_niceUS", adj_pop=adj_pop)
fn_pa_master(dect0="decitree112", dect1="decitree110",
             accuracy_file= "byAgeOnly_niceUS", adj_pop=adj_pop)
fn_pa_master(dect0="decitree101", dect1="decitree110",
             accuracy_file= "byAgeOnly_niceUS", adj_pop=adj_pop)

fn_pa_master(dect0="decitree101var", dect1="decitree112var",
             accuracy_file= "byAgeOnly_niceUS", adj_pop=adj_pop)
fn_pa_master(dect0="decitree112var", dect1="decitree110",
             accuracy_file= "byAgeOnly_niceUS", adj_pop=adj_pop)
fn_pa_master(dect0="decitree101var", dect1="decitree110",
             accuracy_file= "byAgeOnly_niceUS", adj_pop=adj_pop)

#> Stage shift----
# lower effect, i.e. upper bound
fn_pa_master(dect0="decitree100", dect1="decitree101",
             accuracy_file= "byAgeOnly_niceUS", adj_pop=adj_pop, shift_late_sens= "upper")
# fn_pa_master(dect0="decitree100", dect1="decitree112",
#              accuracy_file= "byAgeOnly_niceUS", adj_pop=adj_pop, shift_late_sens= "upper")
# fn_pa_master(dect0="decitree100", dect1="decitree110",
#              accuracy_file= "byAgeOnly_niceUS", adj_pop=adj_pop, shift_late_sens= "upper")
# fn_pa_master(dect0="decitree100", dect1="decitree101var",
#              accuracy_file= "byAgeOnly_niceUS", adj_pop=adj_pop, shift_late_sens= "upper")
fn_pa_master(dect0="decitree101", dect1="decitree110",
             accuracy_file= "byAgeOnly_niceUS", adj_pop=adj_pop, shift_late_sens= "upper")


# higher effect, i.e. lower bound
fn_pa_master(dect0="decitree100", dect1="decitree101",
             accuracy_file= "byAgeOnly_niceUS", adj_pop=adj_pop, shift_late_sens= "lower")
# fn_pa_master(dect0="decitree100", dect1="decitree112",
#              accuracy_file= "byAgeOnly_niceUS", adj_pop=adj_pop, shift_late_sens= "lower")
# fn_pa_master(dect0="decitree100", dect1="decitree110",
#              accuracy_file= "byAgeOnly_niceUS", adj_pop=adj_pop, shift_late_sens= "lower")
# fn_pa_master(dect0="decitree100", dect1="decitree101var",
#              accuracy_file= "byAgeOnly_niceUS", adj_pop=adj_pop, shift_late_sens= "lower")
fn_pa_master(dect0="decitree101", dect1="decitree110",
             accuracy_file= "byAgeOnly_niceUS", adj_pop=adj_pop, shift_late_sens= "lower")

#> FA related utility loss and gain----

#> no short-term disutility
fn_pa_master(dect0="decitree100", dect1="decitree101",
             accuracy_file= "byAgeOnly_niceUS", adj_pop=adj_pop, disuti = FALSE)
# fn_pa_master(dect0="decitree100", dect1="decitree112",
#              accuracy_file= "byAgeOnly_niceUS", adj_pop=adj_pop, disuti = FALSE)
# fn_pa_master(dect0="decitree100", dect1="decitree110",
#              accuracy_file= "byAgeOnly_niceUS", adj_pop=adj_pop, disuti = FALSE)
# fn_pa_master(dect0="decitree100", dect1="decitree101var",
#              accuracy_file= "byAgeOnly_niceUS", adj_pop=adj_pop, disuti = FALSE)
fn_pa_master(dect0="decitree101", dect1="decitree110",
             accuracy_file= "byAgeOnly_niceUS", adj_pop=adj_pop, disuti = FALSE)

#> no long-term benefit
fn_pa_master(dect0="decitree100", dect1="decitree101",
             accuracy_file= "byAgeOnly_niceUS", adj_pop=adj_pop, uti_ben = FALSE)
# fn_pa_master(dect0="decitree100", dect1="decitree112",
#              accuracy_file= "byAgeOnly_niceUS", adj_pop=adj_pop, uti_ben = FALSE)
# fn_pa_master(dect0="decitree100", dect1="decitree110",
#              accuracy_file= "byAgeOnly_niceUS", adj_pop=adj_pop, uti_ben = FALSE)
# fn_pa_master(dect0="decitree100", dect1="decitree101var",
#              accuracy_file= "byAgeOnly_niceUS", adj_pop=adj_pop, uti_ben = FALSE)
fn_pa_master(dect0="decitree101", dect1="decitree110",
             accuracy_file= "byAgeOnly_niceUS", adj_pop=adj_pop, uti_ben = FALSE)

#> Stage related QoL change in a long term----
#> no QoL loss due to late stage
fn_pa_master(dect0="decitree100", dect1="decitree101",
             accuracy_file= "byAgeOnly_niceUS", adj_pop=adj_pop, Markov_snr = "NostageQoL")
# fn_pa_master(dect0="decitree100", dect1="decitree112",
#              accuracy_file= "byAgeOnly_niceUS", adj_pop=adj_pop, Markov_snr = "NostageQoL")
# fn_pa_master(dect0="decitree100", dect1="decitree110",
#              accuracy_file= "byAgeOnly_niceUS", adj_pop=adj_pop, Markov_snr = "NostageQoL")
# fn_pa_master(dect0="decitree100", dect1="decitree101var",
#              accuracy_file= "byAgeOnly_niceUS", adj_pop=adj_pop, Markov_snr = "NostageQoL")
fn_pa_master(dect0="decitree101", dect1="decitree110",
             accuracy_file= "byAgeOnly_niceUS", adj_pop=adj_pop, Markov_snr = "NostageQoL")

#> higher QoL loss due to late stage
fn_pa_master(dect0="decitree100", dect1="decitree101",
             accuracy_file= "byAgeOnly_niceUS", adj_pop=adj_pop, Markov_snr = "stageQoL006")
# fn_pa_master(dect0="decitree100", dect1="decitree112",
#              accuracy_file= "byAgeOnly_niceUS", adj_pop=adj_pop, Markov_snr = "stageQoL006")
# fn_pa_master(dect0="decitree100", dect1="decitree110",
#              accuracy_file= "byAgeOnly_niceUS", adj_pop=adj_pop, Markov_snr = "stageQoL006")
# fn_pa_master(dect0="decitree100", dect1="decitree101var",
#              accuracy_file= "byAgeOnly_niceUS", adj_pop=adj_pop, Markov_snr = "stageQoL006")
fn_pa_master(dect0="decitree101", dect1="decitree110",
             accuracy_file= "byAgeOnly_niceUS", adj_pop=adj_pop, Markov_snr = "stageQoL006")

#> Extending use of model predicted cancer mortality to 15 years----
fn_pa_master(dect0="decitree100", dect1="decitree101",
             accuracy_file= "byAgeOnly_niceUS", adj_pop=adj_pop, Markov_snr = "cd_prob15yr")
# fn_pa_master(dect0="decitree100", dect1="decitree112",
#              accuracy_file= "byAgeOnly_niceUS", adj_pop=adj_pop, Markov_snr = "cd_prob15yr")
# fn_pa_master(dect0="decitree100", dect1="decitree110",
#              accuracy_file= "byAgeOnly_niceUS", adj_pop=adj_pop, Markov_snr = "cd_prob15yr")
# fn_pa_master(dect0="decitree100", dect1="decitree101var",
#              accuracy_file= "byAgeOnly_niceUS", adj_pop=adj_pop, Markov_snr = "cd_prob15yr")
fn_pa_master(dect0="decitree101", dect1="decitree110",
             accuracy_file= "byAgeOnly_niceUS", adj_pop=adj_pop, Markov_snr = "cd_prob15yr")

#> Discount 0.015----
fn_pa_master(dect0="decitree100", dect1="decitree101",
             accuracy_file= "byAgeOnly_niceUS", adj_pop=adj_pop, Markov_snr = "disc015")
# fn_pa_master(dect0="decitree100", dect1="decitree112",
#              accuracy_file= "byAgeOnly_niceUS", adj_pop=adj_pop, Markov_snr = "disc015")
# fn_pa_master(dect0="decitree100", dect1="decitree110",
#              accuracy_file= "byAgeOnly_niceUS", adj_pop=adj_pop, Markov_snr = "disc015")
# fn_pa_master(dect0="decitree100", dect1="decitree101var",
#              accuracy_file= "byAgeOnly_niceUS", adj_pop=adj_pop, Markov_snr = "disc015")
fn_pa_master(dect0="decitree101", dect1="decitree110",
             accuracy_file= "byAgeOnly_niceUS", adj_pop=adj_pop, Markov_snr = "disc015")

#> Varying GP US cost----
# Reduce US cost to CDC cost
fn_pa_master(dect0="decitree100", dect1="decitree101",
             accuracy_file= "byAgeOnly_niceUS", adj_pop=adj_pop, cost_us_cdc=TRUE)
# fn_pa_master(dect0="decitree100", dect1="decitree112",
#              accuracy_file= "byAgeOnly_niceUS", adj_pop=adj_pop, cost_us_cdc=TRUE)
# fn_pa_master(dect0="decitree100", dect1="decitree110",
#              accuracy_file= "byAgeOnly_niceUS", adj_pop=adj_pop, cost_us_cdc=TRUE)
# fn_pa_master(dect0="decitree100", dect1="decitree101var",
#              accuracy_file= "byAgeOnly_niceUS", adj_pop=adj_pop, cost_us_cdc=TRUE)
fn_pa_master(dect0="decitree101", dect1="decitree110",
             accuracy_file= "byAgeOnly_niceUS", adj_pop=adj_pop, cost_us_cdc=TRUE)
# increase US cost to outpatient cost
fn_pa_master(dect0="decitree100", dect1="decitree101",
             accuracy_file= "byAgeOnly_niceUS", adj_pop=adj_pop, cost_us_op=TRUE)
# fn_pa_master(dect0="decitree100", dect1="decitree112",
#              accuracy_file= "byAgeOnly_niceUS", adj_pop=adj_pop, cost_us_op=TRUE)
# fn_pa_master(dect0="decitree100", dect1="decitree110",
#              accuracy_file= "byAgeOnly_niceUS", adj_pop=adj_pop, cost_us_op=TRUE)
# fn_pa_master(dect0="decitree100", dect1="decitree101var",
#              accuracy_file= "byAgeOnly_niceUS", adj_pop=adj_pop, cost_us_op=TRUE)
fn_pa_master(dect0="decitree101", dect1="decitree110",
             accuracy_file= "byAgeOnly_niceUS", adj_pop=adj_pop, cost_us_op=TRUE)

#> by age and stage accuracy----
fn_pa_master(dect0="decitree100", dect1="decitree101",
             accuracy_file= "byAgeStage_niceUS", adj_pop=adj_pop)
# fn_pa_master(dect0="decitree100", dect1="decitree112",
#              accuracy_file= "byAgeStage_niceUS", adj_pop=adj_pop)
# fn_pa_master(dect0="decitree100", dect1="decitree110",
#              accuracy_file= "byAgeStage_niceUS", adj_pop=adj_pop)
# fn_pa_master(dect0="decitree100", dect1="decitree101var",
#              accuracy_file= "byAgeStage_niceUS", adj_pop=adj_pop)
fn_pa_master(dect0="decitree101", dect1="decitree110",
             accuracy_file= "byAgeStage_niceUS", adj_pop=adj_pop)

#> Varying USS accuracy----
# sensitivity plus 0.05
fn_pa_master(dect0="decitree100", dect1="decitree101",
             accuracy_file= "byAgeOnly_niceUS", adj_pop=adj_pop, us_sens_adj_tag="_USsensUP")
# fn_pa_master(dect0="decitree100", dect1="decitree112",
#              accuracy_file= "byAgeOnly_niceUS", adj_pop=adj_pop, us_sens_adj_tag="_USsensUP")
# fn_pa_master(dect0="decitree100", dect1="decitree110",
#              accuracy_file= "byAgeOnly_niceUS", adj_pop=adj_pop, us_sens_adj_tag="_USsensUP")
# fn_pa_master(dect0="decitree100", dect1="decitree101var",
#              accuracy_file= "byAgeOnly_niceUS", adj_pop=adj_pop, us_sens_adj_tag="_USsensUP")
fn_pa_master(dect0="decitree101", dect1="decitree110",
             accuracy_file= "byAgeOnly_niceUS", adj_pop=adj_pop, us_sens_adj_tag="_USsensUP")

# sensitivity minus 0.05
fn_pa_master(dect0="decitree100", dect1="decitree101",
             accuracy_file= "byAgeOnly_niceUS", adj_pop=adj_pop, us_sens_adj_tag="_USsensDOWN")
# fn_pa_master(dect0="decitree100", dect1="decitree112",
#              accuracy_file= "byAgeOnly_niceUS", adj_pop=adj_pop, us_sens_adj_tag="_USsensDOWN")
# fn_pa_master(dect0="decitree100", dect1="decitree110",
#              accuracy_file= "byAgeOnly_niceUS", adj_pop=adj_pop, us_sens_adj_tag="_USsensDOWN")
# fn_pa_master(dect0="decitree100", dect1="decitree101var",
#              accuracy_file= "byAgeOnly_niceUS", adj_pop=adj_pop, us_sens_adj_tag="_USsensDOWN")
fn_pa_master(dect0="decitree101", dect1="decitree110",
             accuracy_file= "byAgeOnly_niceUS", adj_pop=adj_pop, us_sens_adj_tag="_USsensDOWN")

# specificity plus 0.05
fn_pa_master(dect0="decitree100", dect1="decitree101",
             accuracy_file= "byAgeOnly_niceUS", adj_pop=adj_pop, us_spec_adj_tag="_USspecUP")
# fn_pa_master(dect0="decitree100", dect1="decitree112",
#              accuracy_file= "byAgeOnly_niceUS", adj_pop=adj_pop, us_spec_adj_tag="_USspecUP")
# fn_pa_master(dect0="decitree100", dect1="decitree110",
#              accuracy_file= "byAgeOnly_niceUS", adj_pop=adj_pop, us_spec_adj_tag="_USspecUP")
# fn_pa_master(dect0="decitree100", dect1="decitree101var",
#              accuracy_file= "byAgeOnly_niceUS", adj_pop=adj_pop, us_spec_adj_tag="_USspecUP")
fn_pa_master(dect0="decitree101", dect1="decitree110",
             accuracy_file= "byAgeOnly_niceUS", adj_pop=adj_pop, us_spec_adj_tag="_USspecUP")

# specificity minus 0.05
fn_pa_master(dect0="decitree100", dect1="decitree101",
             accuracy_file= "byAgeOnly_niceUS", adj_pop=adj_pop, us_spec_adj_tag="_USspecDOWN")
# fn_pa_master(dect0="decitree100", dect1="decitree112",
#              accuracy_file= "byAgeOnly_niceUS", adj_pop=adj_pop, us_spec_adj_tag="_USspecDOWN")
# fn_pa_master(dect0="decitree100", dect1="decitree110",
#              accuracy_file= "byAgeOnly_niceUS", adj_pop=adj_pop, us_spec_adj_tag="_USspecDOWN")
# fn_pa_master(dect0="decitree100", dect1="decitree101var",
#              accuracy_file= "byAgeOnly_niceUS", adj_pop=adj_pop, us_spec_adj_tag="_USspecDOWN")
fn_pa_master(dect0="decitree101", dect1="decitree110",
             accuracy_file= "byAgeOnly_niceUS", adj_pop=adj_pop, us_spec_adj_tag="_USspecDOWN")


#> Sugery rate varying----
# higher rate
fn_pa_master(dect0="decitree100", dect1="decitree101",
             accuracy_file= "byAgeOnly_niceUS", adj_pop=adj_pop, surg_prob_adj = 0.15)
# fn_pa_master(dect0="decitree100", dect1="decitree112",
#              accuracy_file= "byAgeOnly_niceUS", adj_pop=adj_pop, surg_prob_adj = 0.15)
# fn_pa_master(dect0="decitree100", dect1="decitree110",
#              accuracy_file= "byAgeOnly_niceUS", adj_pop=adj_pop, surg_prob_adj = 0.15)
# fn_pa_master(dect0="decitree100", dect1="decitree101var",
#              accuracy_file= "byAgeOnly_niceUS", adj_pop=adj_pop, surg_prob_adj = 0.15)
fn_pa_master(dect0="decitree101", dect1="decitree110",
             accuracy_file= "byAgeOnly_niceUS", adj_pop=adj_pop, surg_prob_adj = 0.15)

# lower rate
fn_pa_master(dect0="decitree100", dect1="decitree101",
             accuracy_file= "byAgeOnly_niceUS", adj_pop=adj_pop, surg_prob_adj = -0.15)
# fn_pa_master(dect0="decitree100", dect1="decitree112",
#              accuracy_file= "byAgeOnly_niceUS", adj_pop=adj_pop, surg_prob_adj = -0.15)
# fn_pa_master(dect0="decitree100", dect1="decitree110",
#              accuracy_file= "byAgeOnly_niceUS", adj_pop=adj_pop, surg_prob_adj = -0.15)
# fn_pa_master(dect0="decitree100", dect1="decitree101var",
#              accuracy_file= "byAgeOnly_niceUS", adj_pop=adj_pop, surg_prob_adj = -0.15)
fn_pa_master(dect0="decitree101", dect1="decitree110",
             accuracy_file= "byAgeOnly_niceUS", adj_pop=adj_pop, surg_prob_adj = -0.15)


comp <- "pa100v110"
adj <- TRUE
adj_tag <- if (adj) "_adj_pop" else ""

dt <- read.csv(file.path(output, "psa/CEA", comp, 
                         paste0(comp, "_1001_utiben_byAgeOnly_niceUS", adj_tag, ".csv" )))

dt$sum_qaly_disc_dif[dt$group==">=50"]/dt$N[dt$group==">=50"]*1000

dt$sum_qaly_disc_dif[dt$group=="<50"]/dt$N[dt$group=="<50"]*1000

