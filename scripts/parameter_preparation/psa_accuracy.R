##################################################################
##               Generate PSA accuracy parameters               ##
##################################################################
rm(list = ls())
library(bannerCommenter)
library(tidyverse)
library(matrixStats)
source("//Qmcr/dfs/WIPH/Primary Care/OVATOOLS/r_project/scripts/file_path.R")

source(file.path(scripts, "fn_data_clean.R"))


# Data --------------------------------------------------------------------

## CA125 and Ovatools ----

ca125 <- list(
  ##> By age accuracy----
  byAge = list(
    # age 49- 
    sens_ca125_u49 = c(mean = 0.753, lb = 0.700, ub = 0.800),
    spec_ca125_u49 = c(mean = 0.925, lb = 0.923, ub = 0.926),
    
    sens_ovt1_u49 = c(mean = 0.619, lb = 0.561, ub = 0.674),
    spec_ovt1_u49 = c(mean = 0.969, lb = 0.968, ub = 0.970),
    
    sens_ovt3_u49 = c(mean = 0.452, lb = 0.394, ub = 0.510),
    spec_ovt3_u49 = c(mean = 0.994, lb = 0.993, ub = 0.994),
    
    # age 50+
    sens_ca125_o50 = c(mean = 0.865, lb = 0.848, ub = 0.880),
    spec_ca125_o50 = c(mean = 0.943, lb = 0.942, ub = 0.944),
    
    sens_ovt1_o50 = c(mean = 0.911, lb = 0.897, ub = 0.924),
    spec_ovt1_o50 = c(mean = 0.891, lb = 0.889, ub = 0.892),
    
    sens_ovt3_o50 = c(mean = 0.831, lb = 0.813, ub = 0.848),
    spec_ovt3_o50 = c(mean = 0.965, lb = 0.964, ub = 0.966)
  ),
  ##> By stage accuracy----
  byStage = list(
    # early
    sens_ca125_early = c(mean = 0.667, lb = 0.627, ub = 0.706),
    spec_ca125_early = c(mean = 0.936, lb = 0.935, ub = 0.937),
    
    sens_ovt1_early = c(mean = 0.707, lb = 0.668, ub = 0.744),
    spec_ovt1_early = c(mean = 0.925, lb = 0.924, ub = 0.925),
    
    sens_ovt3_early = c(mean = 0.517, lb = 0.476, ub = 0.559),
    spec_ovt3_early = c(mean = 0.978, lb = 0.977, ub = 0.978),
    
    # late
    sens_ca125_late = c(mean = 0.933, lb = 0.918, ub = 0.947),
    spec_ca125_late = c(mean = 0.936, lb = 0.935, ub = 0.937),
    
    sens_ovt1_late = c(mean = 0.946, lb = 0.932, ub = 0.958),
    spec_ovt1_late = c(mean = 0.925, lb = 0.924, ub = 0.925),
    
    sens_ovt3_late = c(mean = 0.896, lb = 0.877, ub = 0.912),
    spec_ovt3_late = c(mean = 0.978, lb = 0.977, ub = 0.978)
  ),
  ##> By age and stage accuracy----
  byAge_stage = list(
    # age 49- 
    # early
    sens_ca125_u49_early = c(mean = 0.633, lb = 0.547, ub = 0.713),
    spec_ca125_u49_early = c(mean = 0.925, lb = 0.924, ub = 0.927),
    
    sens_ovt1_u49_early = c(mean = 0.439, lb = 0.355, ub = 0.525),
    spec_ovt1_u49_early = c(mean = 0.970, lb = 0.969, ub = 0.971),
    
    sens_ovt3_u49_early = c(mean = 0.245, lb = 0.176, ub = 0.325),
    spec_ovt3_u49_early = c(mean = 0.994, lb = 0.993, ub = 0.994),
    
    # late
    sens_ca125_u49_late = c(mean = 0.927, lb = 0.867, ub = 0.966),
    spec_ca125_u49_late = c(mean = 0.925, lb = 0.924, ub = 0.927),
    
    sens_ovt1_u49_late = c(mean = 0.847, lb = 0.771, ub = 0.905),
    spec_ovt1_u49_late = c(mean = 0.970, lb = 0.969, ub = 0.971),
    
    sens_ovt3_u49_late = c(mean = 0.694, lb = 0.604, ub = 0.773),
    spec_ovt3_u49_late = c(mean = 0.994, lb = 0.993, ub = 0.994),
    
    # age 50+
    # early
    sens_ca125_o50_early = c(mean = 0.678, lb = 0.632, ub = 0.721),
    spec_ca125_o50_early = c(mean = 0.944, lb = 0.943, ub = 0.945),
    
    sens_ovt1_o50_early = c(mean = 0.791, lb = 0.750, ub = 0.828),
    spec_ovt1_o50_early = c(mean = 0.892, lb = 0.890, ub = 0.893),
    
    sens_ovt3_o50_early = c(mean = 0.603, lb = 0.556, ub = 0.649),
    spec_ovt3_o50_early = c(mean = 0.966, lb = 0.965, ub = 0.967),
    
    # late
    sens_ca125_o50_late = c(mean = 0.934, lb = 0.918, ub = 0.948),
    spec_ca125_o50_late = c(mean = 0.944, lb = 0.943, ub = 0.945),
    
    sens_ovt1_o50_late = c(mean = 0.957, lb = 0.944, ub = 0.968),
    spec_ovt1_o50_late = c(mean = 0.892, lb = 0.890, ub = 0.893),
    
    sens_ovt3_o50_late = c(mean = 0.918, lb = 0.900, ub = 0.933),
    spec_ovt3_o50_late = c(mean = 0.966, lb = 0.965, ub = 0.967)
  ),
  ##> CA125 threshold varying, by age groups----
  byAge50 = list( # different for ca125 with different threshold by age
    # equavalent to ovt 1%
    
    # age 49-, >=46
    sens_ca125_u49_var = c(mean = 0.635, lb = 0.578, ub = 0.690),
    spec_ca125_u49_var = c(mean = 0.969, lb = 0.968, ub = 0.969),

    # age 49-, >=123, ovt3%
    sens_ca125_u49_var_ovt3 = c(mean = 0.435, lb = 0.378, ub = 0.493),
    spec_ca125_u49_var_ovt3 = c(mean = 0.993, lb = 0.993, ub = 0.994),
    
    # age 50+, threshold varies as below
    # this is overall accuracy
    sens_ca125_o50_var = c(mean = 0.910, lb = 0.896, ub = 0.923),
    spec_ca125_o50_var = c(mean = 0.885, lb = 0.883, ub = 0.886),
    # age 50+, overall, ovt3%
    sens_ca125_o50_var_ovt3 = c(mean = 0.828, lb = 0.810, ub = 0.846),
    spec_ca125_o50_var_ovt3 = c(mean = 0.965, lb = 0.964, ub = 0.966),
    
    # age 50-59, >=26
    sens_ca125_5059 = c(mean = 0.805, lb = 0.763, ub = 0.843),
    spec_ca125_5059 = c(mean = 0.957, lb = 0.937, ub = 0.977), # to correct
    
    sens_ca125_5059_var = c(mean = 0.848, lb = 0.809, ub = 0.882),
    spec_ca125_5059_var = c(mean = 0.916, lb = 0.914, ub = 0.918),
    
    # age 60-69, >=22
    sens_ca125_6069 = c(mean = 0.869, lb = 0.839, ub = 0.895),
    spec_ca125_6069 = c(mean = 0.959, lb = 0.958, ub = 0.961),
    
    sens_ca125_6069_var = c(mean = 0.924, lb = 0.900, ub = 0.944),
    spec_ca125_6069_var = c(mean = 0.893, lb = 0.890, ub = 0.895),
    
    # age 70-79, >=22
    sens_ca125_7079 = c(mean = 0.877, lb = 0.846, ub = 0.903),
    spec_ca125_7079 = c(mean = 0.936, lb = 0.934, ub = 0.938),
    
    sens_ca125_7079_var = c(mean = 0.935, lb = 0.910, ub = 0.954),
    spec_ca125_7079_var = c(mean = 0.847, lb = 0.844, ub = 0.851),
    
    # age 80-89, >=26
    sens_ca125_8089 = c(mean = 0.914, lb = 0.877, ub = 0.943),
    spec_ca125_8089 = c(mean = 0.873, lb = 0.869, ub = 0.878),
    
    sens_ca125_8089_var = c(mean = 0.922, lb = 0.881, ub = 0.951),
    spec_ca125_8089_var = c(mean = 0.818, lb = 0.812, ub = 0.823)
  )
)

# ***----
## USS ----

##> Landolfo 2023----
# recommended by Sudha
uss <- list(
  uss_2step= list(
    ###> Overall----
    overall = list(
      sens_us = c(mean = 0.911, lb = 0.844, ub = 0.951),
      spec_us = c(mean = 0.845, lb = 0.801, ub = 0.881),
      sens_us_ca125 = c(mean = 0.910, lb = 0.845, ub = 0.950),
      spec_us_ca125 = c(mean = 0.856, lb = 0.814, ub = 0.891)
    ),
    ###> By age accuracy----
    byAge = list(
      # age 49- 
      sens_us_u49 = c(mean = 0.875, lb = 0.793, ub = 0.928),
      spec_us_u49 = c(mean = 0.887, lb = 0.851, ub = 0.915),
      sens_us_ca125_u49 = c(mean = 0.869, lb = 0.784, ub = 0.924),
      spec_us_ca125_u49 = c(mean = 0.894, lb = 0.860, ub = 0.921),
      
      # age 50+
      sens_us_o50 = c(mean = 0.927, lb = 0.856, ub = 0.965),
      spec_us_o50 = c(mean = 0.767, lb = 0.702, ub = 0.822),
      sens_us_ca125_o50 = c(mean = 0.920, lb = 0.848, ub = 0.960),
      spec_us_ca125_o50 = c(mean = 0.780, lb = 0.720, ub = 0.830)
    )
  ),
  ##>Sudha's ROCKeTS data----
  uss_sudha = list(
    ###> By age accuracy----
    byAge = list(
      # age 49- 
      sens_us_u49 = c(mean = 0.891, lb = 0.891*(1-0.05), ub = 0.891*(1+0.05)),
      spec_us_u49 = c(mean = 0.751, lb = 0.751*(1-0.05), ub = 0.751*(1+0.05)),
      
      # age 50+
      sens_us_o50 = c(mean = 0.961, lb = 0.922, ub = 0.984),
      spec_us_o50 = c(mean = 0.585, lb = 0.547, ub = 0.621)
    ),
    
    ###> By age and stage accuracy----
    byAge_stage = list(
      
    )
  ),
  
  ##>NICE guidelines data----
  uss_nice = list(
    ###> Overall----
    overall = list(
      sens_us = c(mean = 0.85, lb = 0.83, ub = 0.87),
      spec_us = c(mean = 0.83, lb = 0.81, ub = 0.85)
    )
  )
)


## Others -----
# pathway involving US+
set.seed(1234)
fp_uter_us <- sample_beta2(4+1, 373+46) 
fp_loGI_us <- sample_beta2(1, 373+46)     
fp_uter_ovt3 <- sample_beta2(48, 384)*0.1235       
fp_loGI_ovt3 <- sample_beta2(58, 384)*0.1235
fp_lung_us <- 0
fp_panc_us <- 0
fp_lung_ovt3 <- sample_beta2(49, 384)*0.1235
fp_panc_ovt3 <- sample_beta2(46, 384)*0.1235

# Costs
cost_us=204
cost_us_cdc = 64 # price in directly accessed community diagnostic centre
cost_ca125=10
cost_gp1= 3.99*10 # face-to-face cons 10 min
cost_gp2= 3.99*5 # telephone cons after test 
cost_nurse=9
cost_ct=117# CT of one area, without contrast, 19 years and over
cost_us_op = 209

# ******----

# Monte Carlos sample -----------------------------------------------------

## CA125 accuracy sample----
ca125_sample <- list(
  ##> By age accuracy----
  byAge = list(
    # age 49- 
    sens_ca125_u49 = sample_beta(ca125$byAge$sens_ca125_u49["mean"], 
                                 ca125$byAge$sens_ca125_u49["lb"], 
                                 ca125$byAge$sens_ca125_u49["ub"]),
    
    spec_ca125_u49 = sample_beta(ca125$byAge$spec_ca125_u49["mean"], 
                                 ca125$byAge$spec_ca125_u49["lb"], 
                                 ca125$byAge$spec_ca125_u49["ub"]),
    
    sens_ovt1_u49 = sample_beta(ca125$byAge$sens_ovt1_u49["mean"], 
                                ca125$byAge$sens_ovt1_u49["lb"], 
                                ca125$byAge$sens_ovt1_u49["ub"]),
    
    spec_ovt1_u49 = sample_beta(ca125$byAge$spec_ovt1_u49["mean"], 
                                ca125$byAge$spec_ovt1_u49["lb"], 
                                ca125$byAge$spec_ovt1_u49["ub"]),
    
    sens_ovt3_u49 = sample_beta(ca125$byAge$sens_ovt3_u49["mean"], 
                                ca125$byAge$sens_ovt3_u49["lb"], 
                                ca125$byAge$sens_ovt3_u49["ub"]),
    
    spec_ovt3_u49 = sample_beta(ca125$byAge$spec_ovt3_u49["mean"], 
                                ca125$byAge$spec_ovt3_u49["lb"], 
                                ca125$byAge$spec_ovt3_u49["ub"]),
    
    # age 50+
    sens_ca125_o50 = sample_beta(ca125$byAge$sens_ca125_o50["mean"], 
                                 ca125$byAge$sens_ca125_o50["lb"], 
                                 ca125$byAge$sens_ca125_o50["ub"]),
    
    spec_ca125_o50 = sample_beta(ca125$byAge$spec_ca125_o50["mean"], 
                                 ca125$byAge$spec_ca125_o50["lb"], 
                                 ca125$byAge$spec_ca125_o50["ub"]),
    
    sens_ovt1_o50 = sample_beta(ca125$byAge$sens_ovt1_o50["mean"], 
                                ca125$byAge$sens_ovt1_o50["lb"], 
                                ca125$byAge$sens_ovt1_o50["ub"]),
    
    spec_ovt1_o50 = sample_beta(ca125$byAge$spec_ovt1_o50["mean"], 
                                ca125$byAge$spec_ovt1_o50["lb"], 
                                ca125$byAge$spec_ovt1_o50["ub"]),
    
    sens_ovt3_o50 = sample_beta(ca125$byAge$sens_ovt3_o50["mean"], 
                                ca125$byAge$sens_ovt3_o50["lb"], 
                                ca125$byAge$sens_ovt3_o50["ub"]),
    
    spec_ovt3_o50 = sample_beta(ca125$byAge$spec_ovt3_o50["mean"], 
                                ca125$byAge$spec_ovt3_o50["lb"], 
                                ca125$byAge$spec_ovt3_o50["ub"])
  ),
  ##> By stage accuracy----
  byStage = list(
    # early
    sens_ca125_early = sample_beta(ca125$byStage$sens_ca125_early["mean"], 
                                   ca125$byStage$sens_ca125_early["lb"], 
                                   ca125$byStage$sens_ca125_early["ub"]),
    
    spec_ca125_early = sample_beta(ca125$byStage$spec_ca125_early["mean"], 
                                   ca125$byStage$spec_ca125_early["lb"], 
                                   ca125$byStage$spec_ca125_early["ub"]),
    
    sens_ovt1_early = sample_beta(ca125$byStage$sens_ovt1_early["mean"], 
                                  ca125$byStage$sens_ovt1_early["lb"], 
                                  ca125$byStage$sens_ovt1_early["ub"]),
    
    spec_ovt1_early = sample_beta(ca125$byStage$spec_ovt1_early["mean"], 
                                  ca125$byStage$spec_ovt1_early["lb"], 
                                  ca125$byStage$spec_ovt1_early["ub"]),
    
    sens_ovt3_early = sample_beta(ca125$byStage$sens_ovt3_early["mean"], 
                                  ca125$byStage$sens_ovt3_early["lb"], 
                                  ca125$byStage$sens_ovt3_early["ub"]),
    
    spec_ovt3_early = sample_beta(ca125$byStage$spec_ovt3_early["mean"], 
                                  ca125$byStage$spec_ovt3_early["lb"], 
                                  ca125$byStage$spec_ovt3_early["ub"]),
    
    # late
    sens_ca125_late = sample_beta(ca125$byStage$sens_ca125_late["mean"], 
                                  ca125$byStage$sens_ca125_late["lb"], 
                                  ca125$byStage$sens_ca125_late["ub"]),
    
    spec_ca125_late = sample_beta(ca125$byStage$spec_ca125_late["mean"], 
                                  ca125$byStage$spec_ca125_late["lb"], 
                                  ca125$byStage$spec_ca125_late["ub"]),
    
    sens_ovt1_late = sample_beta(ca125$byStage$sens_ovt1_late["mean"], 
                                 ca125$byStage$sens_ovt1_late["lb"], 
                                 ca125$byStage$sens_ovt1_late["ub"]),
    
    spec_ovt1_late = sample_beta(ca125$byStage$spec_ovt1_late["mean"], 
                                 ca125$byStage$spec_ovt1_late["lb"], 
                                 ca125$byStage$spec_ovt1_late["ub"]),
    
    sens_ovt3_late = sample_beta(ca125$byStage$sens_ovt3_late["mean"], 
                                 ca125$byStage$sens_ovt3_late["lb"], 
                                 ca125$byStage$sens_ovt3_late["ub"]),
    
    spec_ovt3_late = sample_beta(ca125$byStage$spec_ovt3_late["mean"], 
                                 ca125$byStage$spec_ovt3_late["lb"], 
                                 ca125$byStage$spec_ovt3_late["ub"])
  ),
  ##> By age and stage accuracy----
  byAge_stage = list(
    # age 49- 
    # early
    sens_ca125_u49_early = sample_beta(ca125$byAge_stage$sens_ca125_u49_early["mean"], 
                                       ca125$byAge_stage$sens_ca125_u49_early["lb"], 
                                       ca125$byAge_stage$sens_ca125_u49_early["ub"]),
    
    spec_ca125_u49_early = sample_beta(ca125$byAge_stage$spec_ca125_u49_early["mean"], 
                                       ca125$byAge_stage$spec_ca125_u49_early["lb"], 
                                       ca125$byAge_stage$spec_ca125_u49_early["ub"]),
    
    sens_ovt1_u49_early = sample_beta(ca125$byAge_stage$sens_ovt1_u49_early["mean"], 
                                      ca125$byAge_stage$sens_ovt1_u49_early["lb"], 
                                      ca125$byAge_stage$sens_ovt1_u49_early["ub"]),
    
    spec_ovt1_u49_early = sample_beta(ca125$byAge_stage$spec_ovt1_u49_early["mean"], 
                                      ca125$byAge_stage$spec_ovt1_u49_early["lb"], 
                                      ca125$byAge_stage$spec_ovt1_u49_early["ub"]),
    
    sens_ovt3_u49_early = sample_beta(ca125$byAge_stage$sens_ovt3_u49_early["mean"], 
                                      ca125$byAge_stage$sens_ovt3_u49_early["lb"], 
                                      ca125$byAge_stage$sens_ovt3_u49_early["ub"]),
    
    spec_ovt3_u49_early = sample_beta(ca125$byAge_stage$spec_ovt3_u49_early["mean"], 
                                      ca125$byAge_stage$spec_ovt3_u49_early["lb"], 
                                      ca125$byAge_stage$spec_ovt3_u49_early["ub"]),
    # late
    sens_ca125_u49_late = sample_beta(ca125$byAge_stage$sens_ca125_u49_late["mean"], 
                                      ca125$byAge_stage$sens_ca125_u49_late["lb"], 
                                      ca125$byAge_stage$sens_ca125_u49_late["ub"]),
    
    spec_ca125_u49_late = sample_beta(ca125$byAge_stage$spec_ca125_u49_late["mean"], 
                                      ca125$byAge_stage$spec_ca125_u49_late["lb"], 
                                      ca125$byAge_stage$spec_ca125_u49_late["ub"]),
    
    sens_ovt1_u49_late = sample_beta(ca125$byAge_stage$sens_ovt1_u49_late["mean"], 
                                     ca125$byAge_stage$sens_ovt1_u49_late["lb"], 
                                     ca125$byAge_stage$sens_ovt1_u49_late["ub"]),
    
    spec_ovt1_u49_late = sample_beta(ca125$byAge_stage$spec_ovt1_u49_late["mean"], 
                                     ca125$byAge_stage$spec_ovt1_u49_late["lb"], 
                                     ca125$byAge_stage$spec_ovt1_u49_late["ub"]),
    
    sens_ovt3_u49_late = sample_beta(ca125$byAge_stage$sens_ovt3_u49_late["mean"], 
                                     ca125$byAge_stage$sens_ovt3_u49_late["lb"], 
                                     ca125$byAge_stage$sens_ovt3_u49_late["ub"]),
    
    spec_ovt3_u49_late = sample_beta(ca125$byAge_stage$spec_ovt3_u49_late["mean"], 
                                     ca125$byAge_stage$spec_ovt3_u49_late["lb"], 
                                     ca125$byAge_stage$spec_ovt3_u49_late["ub"]),
    
    # age 50+
    # early
    sens_ca125_o50_early = sample_beta(ca125$byAge_stage$sens_ca125_o50_early["mean"], 
                                       ca125$byAge_stage$sens_ca125_o50_early["lb"], 
                                       ca125$byAge_stage$sens_ca125_o50_early["ub"]),
    
    spec_ca125_o50_early = sample_beta(ca125$byAge_stage$spec_ca125_o50_early["mean"], 
                                       ca125$byAge_stage$spec_ca125_o50_early["lb"], 
                                       ca125$byAge_stage$spec_ca125_o50_early["ub"]),
    
    sens_ovt1_o50_early = sample_beta(ca125$byAge_stage$sens_ovt1_o50_early["mean"], 
                                      ca125$byAge_stage$sens_ovt1_o50_early["lb"], 
                                      ca125$byAge_stage$sens_ovt1_o50_early["ub"]),
    
    spec_ovt1_o50_early = sample_beta(ca125$byAge_stage$spec_ovt1_o50_early["mean"], 
                                      ca125$byAge_stage$spec_ovt1_o50_early["lb"], 
                                      ca125$byAge_stage$spec_ovt1_o50_early["ub"]),
    
    sens_ovt3_o50_early = sample_beta(ca125$byAge_stage$sens_ovt3_o50_early["mean"], 
                                      ca125$byAge_stage$sens_ovt3_o50_early["lb"], 
                                      ca125$byAge_stage$sens_ovt3_o50_early["ub"]),
    
    spec_ovt3_o50_early = sample_beta(ca125$byAge_stage$spec_ovt3_o50_early["mean"], 
                                      ca125$byAge_stage$spec_ovt3_o50_early["lb"], 
                                      ca125$byAge_stage$spec_ovt3_o50_early["ub"]),
    # late
    sens_ca125_o50_late = sample_beta(ca125$byAge_stage$sens_ca125_o50_late["mean"], 
                                      ca125$byAge_stage$sens_ca125_o50_late["lb"], 
                                      ca125$byAge_stage$sens_ca125_o50_late["ub"]),
    
    spec_ca125_o50_late = sample_beta(ca125$byAge_stage$spec_ca125_o50_late["mean"], 
                                      ca125$byAge_stage$spec_ca125_o50_late["lb"], 
                                      ca125$byAge_stage$spec_ca125_o50_late["ub"]),
    
    sens_ovt1_o50_late = sample_beta(ca125$byAge_stage$sens_ovt1_o50_late["mean"], 
                                     ca125$byAge_stage$sens_ovt1_o50_late["lb"], 
                                     ca125$byAge_stage$sens_ovt1_o50_late["ub"]),
    
    spec_ovt1_o50_late = sample_beta(ca125$byAge_stage$spec_ovt1_o50_late["mean"], 
                                     ca125$byAge_stage$spec_ovt1_o50_late["lb"], 
                                     ca125$byAge_stage$spec_ovt1_o50_late["ub"]),
    
    sens_ovt3_o50_late = sample_beta(ca125$byAge_stage$sens_ovt3_o50_late["mean"], 
                                     ca125$byAge_stage$sens_ovt3_o50_late["lb"], 
                                     ca125$byAge_stage$sens_ovt3_o50_late["ub"]),
    
    spec_ovt3_o50_late = sample_beta(ca125$byAge_stage$spec_ovt3_o50_late["mean"], 
                                     ca125$byAge_stage$spec_ovt3_o50_late["lb"], 
                                     ca125$byAge_stage$spec_ovt3_o50_late["ub"])
  ),
  ##> By age50, with varying threshold for CA125, corresponding to 1% of Ovatools----
  byAge50 = list(
    # age 49- 
    sens_ca125_u49_var = sample_beta(ca125$byAge50$sens_ca125_u49_var["mean"], 
                                 ca125$byAge50$sens_ca125_u49_var["lb"], 
                                 ca125$byAge50$sens_ca125_u49_var["ub"]),
    
    spec_ca125_u49_var = sample_beta(ca125$byAge50$spec_ca125_u49_var["mean"], 
                                 ca125$byAge50$spec_ca125_u49_var["lb"], 
                                 ca125$byAge50$spec_ca125_u49_var["ub"]),
    
    sens_ca125_u49_var_ovt3 = sample_beta(ca125$byAge50$sens_ca125_u49_var_ovt3["mean"], 
                                     ca125$byAge50$sens_ca125_u49_var_ovt3["lb"], 
                                     ca125$byAge50$sens_ca125_u49_var_ovt3["ub"]),
    
    spec_ca125_u49_var_ovt3 = sample_beta(ca125$byAge50$spec_ca125_u49_var_ovt3["mean"], 
                                     ca125$byAge50$spec_ca125_u49_var_ovt3["lb"], 
                                     ca125$byAge50$spec_ca125_u49_var_ovt3["ub"]),
    
    # age 50+, threshold varies as below
    # this is overall accuracy
    sens_ca125_o50_var = sample_beta(ca125$byAge50$sens_ca125_o50_var["mean"], 
                                     ca125$byAge50$sens_ca125_o50_var["lb"], 
                                     ca125$byAge50$sens_ca125_o50_var["ub"]),
    
    spec_ca125_o50_var = sample_beta(ca125$byAge50$spec_ca125_o50_var["mean"], 
                                     ca125$byAge50$spec_ca125_o50_var["lb"], 
                                     ca125$byAge50$spec_ca125_o50_var["ub"]),
    
    sens_ca125_o50_var_ovt3 = sample_beta(ca125$byAge50$sens_ca125_o50_var_ovt3["mean"], 
                                     ca125$byAge50$sens_ca125_o50_var_ovt3["lb"], 
                                     ca125$byAge50$sens_ca125_o50_var_ovt3["ub"]),
    
    spec_ca125_o50_var_ovt3 = sample_beta(ca125$byAge50$spec_ca125_o50_var_ovt3["mean"], 
                                     ca125$byAge50$spec_ca125_o50_var_ovt3["lb"], 
                                     ca125$byAge50$spec_ca125_o50_var_ovt3["ub"]),
    
    # age 50-59
    sens_ca125_5059 = sample_beta(ca125$byAge50$sens_ca125_5059["mean"], 
                                      ca125$byAge50$sens_ca125_5059["lb"], 
                                      ca125$byAge50$sens_ca125_5059["ub"]),
    
    spec_ca125_5059 = sample_beta(ca125$byAge50$spec_ca125_5059["mean"], 
                                      ca125$byAge50$spec_ca125_5059["lb"], 
                                      ca125$byAge50$spec_ca125_5059["ub"]),
    
    sens_ca125_5059_var = sample_beta(ca125$byAge50$sens_ca125_5059_var["mean"], 
                                  ca125$byAge50$sens_ca125_5059_var["lb"], 
                                  ca125$byAge50$sens_ca125_5059_var["ub"]),
    
    spec_ca125_5059_var = sample_beta(ca125$byAge50$spec_ca125_5059_var["mean"], 
                                  ca125$byAge50$spec_ca125_5059_var["lb"], 
                                  ca125$byAge50$spec_ca125_5059_var["ub"]),
    
    # 60-69
    sens_ca125_6069 = sample_beta(ca125$byAge50$sens_ca125_6069["mean"], 
                                      ca125$byAge50$sens_ca125_6069["lb"], 
                                      ca125$byAge50$sens_ca125_6069["ub"]),
    
    spec_ca125_6069 = sample_beta(ca125$byAge50$spec_ca125_6069["mean"], 
                                      ca125$byAge50$spec_ca125_6069["lb"], 
                                      ca125$byAge50$spec_ca125_6069["ub"]),
    
    
    sens_ca125_6069_var = sample_beta(ca125$byAge50$sens_ca125_6069_var["mean"], 
                                  ca125$byAge50$sens_ca125_6069_var["lb"], 
                                  ca125$byAge50$sens_ca125_6069_var["ub"]),
    
    spec_ca125_6069_var = sample_beta(ca125$byAge50$spec_ca125_6069_var["mean"], 
                                  ca125$byAge50$spec_ca125_6069_var["lb"], 
                                  ca125$byAge50$spec_ca125_6069_var["ub"]),
    # 70-79
    sens_ca125_7079 = sample_beta(ca125$byAge50$sens_ca125_7079["mean"], 
                                      ca125$byAge50$sens_ca125_7079["lb"], 
                                      ca125$byAge50$sens_ca125_7079["ub"]),
    
    spec_ca125_7079 = sample_beta(ca125$byAge50$spec_ca125_7079["mean"], 
                                      ca125$byAge50$spec_ca125_7079["lb"], 
                                      ca125$byAge50$spec_ca125_7079["ub"]),
    
    sens_ca125_7079_var = sample_beta(ca125$byAge50$sens_ca125_7079_var["mean"], 
                                  ca125$byAge50$sens_ca125_7079_var["lb"], 
                                  ca125$byAge50$sens_ca125_7079_var["ub"]),
    
    spec_ca125_7079_var = sample_beta(ca125$byAge50$spec_ca125_7079_var["mean"], 
                                  ca125$byAge50$spec_ca125_7079_var["lb"], 
                                  ca125$byAge50$spec_ca125_7079_var["ub"]),
    # 80-89
    sens_ca125_8089 = sample_beta(ca125$byAge50$sens_ca125_8089["mean"], 
                                      ca125$byAge50$sens_ca125_8089["lb"], 
                                      ca125$byAge50$sens_ca125_8089["ub"]),
    
    spec_ca125_8089 = sample_beta(ca125$byAge50$spec_ca125_8089["mean"], 
                                      ca125$byAge50$spec_ca125_8089["lb"], 
                                      ca125$byAge50$spec_ca125_8089["ub"]), 
    
    sens_ca125_8089_var = sample_beta(ca125$byAge50$sens_ca125_8089_var["mean"], 
                                  ca125$byAge50$sens_ca125_8089_var["lb"], 
                                  ca125$byAge50$sens_ca125_8089_var["ub"]),
    
    spec_ca125_8089_var = sample_beta(ca125$byAge50$spec_ca125_8089_var["mean"], 
                                  ca125$byAge50$spec_ca125_8089_var["lb"], 
                                  ca125$byAge50$spec_ca125_8089_var["ub"])
  )
)
# ***----
## USS accuracy sample----
uss_sample <- list(
  
  ##> 2 step accuracy----
  uss_2step= list(
    ###> Overall----
    overall = list(
      sens_us = sample_beta(uss$uss_2step$overall$sens_us["mean"],
                            uss$uss_2step$overall$sens_us["lb"],
                            uss$uss_2step$overall$sens_us["ub"]),
      
      spec_us = sample_beta(uss$uss_2step$overall$spec_us["mean"],
                            uss$uss_2step$overall$spec_us["lb"],
                            uss$uss_2step$overall$spec_us["ub"]),
      
      sens_us_ca125 = sample_beta(uss$uss_2step$overall$sens_us_ca125["mean"],
                                  uss$uss_2step$overall$sens_us_ca125["lb"],
                                  uss$uss_2step$overall$sens_us_ca125["ub"]),
      
      spec_us_ca125 = sample_beta(uss$uss_2step$overall$spec_us_ca125["mean"],
                                  uss$uss_2step$overall$spec_us_ca125["lb"],
                                  uss$uss_2step$overall$spec_us_ca125["ub"])
    ),
    ###> By age accuracy----
    byAge = list(
      # age 49- 
      sens_us_u49 = sample_beta(uss$uss_2step$byAge$sens_us_u49["mean"],
                                uss$uss_2step$byAge$sens_us_u49["lb"],
                                uss$uss_2step$byAge$sens_us_u49["ub"]),
      
      spec_us_u49 = sample_beta(uss$uss_2step$byAge$spec_us_u49["mean"],
                                uss$uss_2step$byAge$spec_us_u49["lb"],
                                uss$uss_2step$byAge$spec_us_u49["ub"]),
      
      sens_us_ca125_u49 = sample_beta(uss$uss_2step$byAge$sens_us_ca125_u49["mean"],
                                      uss$uss_2step$byAge$sens_us_ca125_u49["lb"],
                                      uss$uss_2step$byAge$sens_us_ca125_u49["ub"]),
      
      spec_us_ca125_u49 = sample_beta(uss$uss_2step$byAge$spec_us_ca125_u49["mean"],
                                      uss$uss_2step$byAge$spec_us_ca125_u49["lb"],
                                      uss$uss_2step$byAge$spec_us_ca125_u49["ub"]),
      
      # age 50+
      sens_us_o50 = sample_beta(uss$uss_2step$byAge$sens_us_o50["mean"],
                                uss$uss_2step$byAge$sens_us_o50["lb"],
                                uss$uss_2step$byAge$sens_us_o50["ub"]),
      
      spec_us_o50 = sample_beta(uss$uss_2step$byAge$spec_us_o50["mean"],
                                uss$uss_2step$byAge$spec_us_o50["lb"],
                                uss$uss_2step$byAge$spec_us_o50["ub"]),
      
      sens_us_ca125_o50 = sample_beta(uss$uss_2step$byAge$sens_us_ca125_o50["mean"],
                                      uss$uss_2step$byAge$sens_us_ca125_o50["lb"],
                                      uss$uss_2step$byAge$sens_us_ca125_o50["ub"]),
      
      spec_us_ca125_o50 = sample_beta(uss$uss_2step$byAge$spec_us_ca125_o50["mean"],
                                      uss$uss_2step$byAge$spec_us_ca125_o50["lb"],
                                      uss$uss_2step$byAge$spec_us_ca125_o50["ub"])
    )
  ),
  ##> Sudha's ROCKeTS data----
  uss_sudha = list(
    ###> By age accuracy----
    byAge = list(
      # age 49- 
      sens_us_u49 = sample_beta(uss$uss_sudha$byAge$sens_us_u49["mean"],
                                uss$uss_sudha$byAge$sens_us_u49["lb"],
                                uss$uss_sudha$byAge$sens_us_u49["ub"]),
      
      spec_us_u49 = sample_beta(uss$uss_sudha$byAge$spec_us_u49["mean"],
                                uss$uss_sudha$byAge$spec_us_u49["lb"],
                                uss$uss_sudha$byAge$spec_us_u49["ub"]),
      
      # age 50+
      sens_us_o50 = sample_beta(uss$uss_sudha$byAge$sens_us_o50["mean"],
                                uss$uss_sudha$byAge$sens_us_o50["lb"],
                                uss$uss_sudha$byAge$sens_us_o50["ub"]),
      
      spec_us_o50 = sample_beta(uss$uss_sudha$byAge$spec_us_o50["mean"],
                                uss$uss_sudha$byAge$spec_us_o50["lb"],
                                uss$uss_sudha$byAge$spec_us_o50["ub"])
    ),
    
    ###> By age and stage accuracy----
    byAge_stage = list(
      
    )
  ),
  
  ##> NICE guidelines data----
  uss_nice = list(
    ###> Overall----
    overall = list(
      sens_us = sample_beta(uss$uss_nice$overall$sens_us["mean"],
                            uss$uss_nice$overall$sens_us["lb"],
                            uss$uss_nice$overall$sens_us["ub"]),
      
      spec_us = sample_beta(uss$uss_nice$overall$spec_us["mean"],
                            uss$uss_nice$overall$spec_us["lb"],
                            uss$uss_nice$overall$spec_us["ub"])
    )
  )
)

# ******----
# Create file -------------------------------------------------------------

## Use NICE US data----
##> By age CA125/Ovatools accuracy----

N <- 1001
accuracy <- list()
for (i in 1:N) {
  
  accuracy[[i]] <- list(
    
    # 49 or under
    u49 = list(
      accuracy_early = c(
        sens_ca125 = ca125_sample$byAge$sens_ca125_u49[i],
        spec_ca125 = ca125_sample$byAge$spec_ca125_u49[i],
        sens_ovt1 = ca125_sample$byAge$sens_ovt1_u49[i],
        spec_ovt1 = ca125_sample$byAge$spec_ovt1_u49[i],
        sens_ovt3 = ca125_sample$byAge$sens_ovt3_u49[i],
        spec_ovt3 = ca125_sample$byAge$spec_ovt3_u49[i],
        sens_ovt3_fl_us = ca125_sample$byAge$sens_ovt3_u49[i],
        spec_ovt3_fl_us = ca125_sample$byAge$spec_ovt3_u49[i],
        
        sens_us = uss_sample$uss_nice$overall$sens_us[i],
        spec_us = uss_sample$uss_nice$overall$spec_us[i],
        sens_us_fl_ca125 = uss_sample$uss_nice$overall$sens_us[i],
        spec_us_fl_ca125 = uss_sample$uss_nice$overall$spec_us[i],
        sens_us_fl_ovt = uss_sample$uss_nice$overall$sens_us[i],
        spec_us_fl_ovt = uss_sample$uss_nice$overall$spec_us[i],
        
        fp_uter_us = fp_uter_us[i], 
        fp_loGI_us = fp_loGI_us[i], 
        fp_uter_ovt3 = fp_uter_ovt3[i],       
        fp_loGI_ovt3 = fp_loGI_ovt3[i], 
        fp_lung_us=fp_lung_us,
        fp_panc_us=fp_panc_us,
        fp_lung_ovt3 = fp_lung_ovt3[i], 
        fp_panc_ovt3 = fp_panc_ovt3[i], 
        
        # no uncertainty for service item costs
        cost_us = cost_us,
        cost_us_cdc = cost_us_cdc,
        cost_us_op = cost_us_op,
        cost_ca125 = cost_ca125,
        cost_gp1 = cost_gp1, # face-to-face cons
        cost_gp2 = cost_gp2, # telephone cons after test 
        cost_nurse =cost_nurse,
        cost_ct = cost_ct
      )
    ),
    
    # 51 or over
    o50 = list(
      accuracy_early = c(
        sens_ca125 = ca125_sample$byAge$sens_ca125_o50[i],
        spec_ca125 = ca125_sample$byAge$spec_ca125_o50[i],
        sens_ovt1 = ca125_sample$byAge$sens_ovt1_o50[i],
        spec_ovt1 = ca125_sample$byAge$spec_ovt1_o50[i],
        sens_ovt3 = ca125_sample$byAge$sens_ovt3_o50[i],
        spec_ovt3 = ca125_sample$byAge$spec_ovt3_o50[i],
        sens_ovt3_fl_us = ca125_sample$byAge$sens_ovt3_o50[i],
        spec_ovt3_fl_us = ca125_sample$byAge$spec_ovt3_o50[i],
        
        sens_us = uss_sample$uss_nice$overall$sens_us[i],
        spec_us = uss_sample$uss_nice$overall$spec_us[i],
        sens_us_fl_ca125 = uss_sample$uss_nice$overall$sens_us[i],
        spec_us_fl_ca125 = uss_sample$uss_nice$overall$spec_us[i],
        sens_us_fl_ovt = uss_sample$uss_nice$overall$sens_us[i],
        spec_us_fl_ovt = uss_sample$uss_nice$overall$spec_us[i],
        
        fp_uter_us = fp_uter_us[i], 
        fp_loGI_us = fp_loGI_us[i], 
        fp_uter_ovt3 = fp_uter_ovt3[i],       
        fp_loGI_ovt3 = fp_loGI_ovt3[i], 
        fp_lung_us=fp_lung_us,
        fp_panc_us=fp_panc_us,
        fp_lung_ovt3 = fp_lung_ovt3[i], 
        fp_panc_ovt3 = fp_panc_ovt3[i], 
        
        # no uncertainty for service item costs
        cost_us = cost_us,
        cost_us_cdc = cost_us_cdc,
        cost_us_op = cost_us_op,
        cost_ca125 = cost_ca125,
        cost_gp1 = cost_gp1, # face-to-face cons
        cost_gp2 = cost_gp2, # telephone cons after test 
        cost_nurse =cost_nurse,
        cost_ct = cost_ct
      )
    )
  )
}

for(i in 1:N){
  
  accuracy[[i]]$u49$accuracy_late <- accuracy[[i]]$u49$accuracy_early
  accuracy[[i]]$o50$accuracy_late <- accuracy[[i]]$o50$accuracy_early
  
}

acc_u49 <- lapply(accuracy, function(x) x[[1]])
acc_o50 <- lapply(accuracy, function(x) x[[2]])
# 
test_psa3(accuracy = acc_u49)
test_psa3(accuracy = acc_o50)

saveRDS(accuracy, file.path(work_data, "psa", "accuracy_psa_byAgeOnly_niceUS.rds"))


##> By age and stage CA125/Ovatools data----

N <- 1001
accuracy <- list()
for (i in 1:N) {
  
  accuracy[[i]] <- list(
    
    # 49 or under
    u49 = list(
      accuracy_early = c(
        sens_ca125 = ca125_sample$byAge_stage$sens_ca125_u49_early[i],
        spec_ca125 = ca125_sample$byAge_stage$spec_ca125_u49_early[i],
        sens_ovt1 = ca125_sample$byAge_stage$sens_ovt1_u49_early[i],
        spec_ovt1 = ca125_sample$byAge_stage$spec_ovt1_u49_early[i],
        sens_ovt3 = ca125_sample$byAge_stage$sens_ovt3_u49_early[i],
        spec_ovt3 = ca125_sample$byAge_stage$spec_ovt3_u49_early[i],
        
        sens_ovt3_fl_us = ca125_sample$byAge_stage$sens_ovt3_u49_early[i],
        spec_ovt3_fl_us = ca125_sample$byAge_stage$spec_ovt3_u49_early[i],
        
        sens_us = uss_sample$uss_nice$overall$sens_us[i],
        spec_us = uss_sample$uss_nice$overall$spec_us[i],
        sens_us_fl_ca125 = uss_sample$uss_nice$overall$sens_us[i],
        spec_us_fl_ca125 = uss_sample$uss_nice$overall$spec_us[i],
        sens_us_fl_ovt = uss_sample$uss_nice$overall$sens_us[i],
        spec_us_fl_ovt = uss_sample$uss_nice$overall$spec_us[i],
        
        fp_uter_us = fp_uter_us[i], 
        fp_loGI_us = fp_loGI_us[i], 
        fp_uter_ovt3 = fp_uter_ovt3[i],       
        fp_loGI_ovt3 = fp_loGI_ovt3[i], 
        fp_lung_us=fp_lung_us,
        fp_panc_us=fp_panc_us,
        fp_lung_ovt3 = fp_lung_ovt3[i], 
        fp_panc_ovt3 = fp_panc_ovt3[i], 
        
        # no uncertainty for service item costs
        cost_us = cost_us,
        cost_us_cdc = cost_us_cdc,
        cost_us_op = cost_us_op,
        cost_ca125 = cost_ca125,
        cost_gp1 = cost_gp1, # face-to-face cons
        cost_gp2 = cost_gp2, # telephone cons after test 
        cost_nurse =cost_nurse,
        cost_ct = cost_ct
      ),
      accuracy_late = c(
        sens_ca125 = ca125_sample$byAge_stage$sens_ca125_u49_late[i],
        spec_ca125 = ca125_sample$byAge_stage$spec_ca125_u49_late[i],
        sens_ovt1 = ca125_sample$byAge_stage$sens_ovt1_u49_late[i],
        spec_ovt1 = ca125_sample$byAge_stage$spec_ovt1_u49_late[i],
        sens_ovt3 = ca125_sample$byAge_stage$sens_ovt3_u49_late[i],
        spec_ovt3 = ca125_sample$byAge_stage$spec_ovt3_u49_late[i],
        
        sens_ovt3_fl_us = ca125_sample$byAge_stage$sens_ovt3_u49_late[i],
        spec_ovt3_fl_us = ca125_sample$byAge_stage$spec_ovt3_u49_late[i],
        
        sens_us = uss_sample$uss_nice$overall$sens_us[i],
        spec_us = uss_sample$uss_nice$overall$spec_us[i],
        sens_us_fl_ca125 = uss_sample$uss_nice$overall$sens_us[i],
        spec_us_fl_ca125 = uss_sample$uss_nice$overall$spec_us[i],
        sens_us_fl_ovt = uss_sample$uss_nice$overall$sens_us[i],
        spec_us_fl_ovt = uss_sample$uss_nice$overall$spec_us[i],
        
        fp_uter_us = fp_uter_us[i], 
        fp_loGI_us = fp_loGI_us[i], 
        fp_uter_ovt3 = fp_uter_ovt3[i],       
        fp_loGI_ovt3 = fp_loGI_ovt3[i], 
        fp_lung_us=fp_lung_us,
        fp_panc_us=fp_panc_us,
        fp_lung_ovt3 = fp_lung_ovt3[i], 
        fp_panc_ovt3 = fp_panc_ovt3[i], 
        
        # no uncertainty for service item costs
        cost_us = cost_us,
        cost_us_cdc = cost_us_cdc,
        cost_us_op = cost_us_op,
        cost_ca125 = cost_ca125,
        cost_gp1 = cost_gp1, # face-to-face cons
        cost_gp2 = cost_gp2, # telephone cons after test 
        cost_nurse =cost_nurse,
        cost_ct = cost_ct
      )
    ),
    
    # 51 or over
    o50 = list(
      accuracy_early = c(
        sens_ca125 = ca125_sample$byAge_stage$sens_ca125_o50_early[i],
        spec_ca125 = ca125_sample$byAge_stage$spec_ca125_o50_early[i],
        sens_ovt1 = ca125_sample$byAge_stage$sens_ovt1_o50_early[i],
        spec_ovt1 = ca125_sample$byAge_stage$spec_ovt1_o50_early[i],
        sens_ovt3 = ca125_sample$byAge_stage$sens_ovt3_o50_early[i],
        spec_ovt3 = ca125_sample$byAge_stage$spec_ovt3_o50_early[i],
        sens_ovt3_fl_us = ca125_sample$byAge_stage$sens_ovt3_o50_early[i],
        spec_ovt3_fl_us = ca125_sample$byAge_stage$spec_ovt3_o50_early[i],
        
        sens_us = uss_sample$uss_nice$overall$sens_us[i],
        spec_us = uss_sample$uss_nice$overall$spec_us[i],
        sens_us_fl_ca125 = uss_sample$uss_nice$overall$sens_us[i],
        spec_us_fl_ca125 = uss_sample$uss_nice$overall$spec_us[i],
        sens_us_fl_ovt = uss_sample$uss_nice$overall$sens_us[i],
        spec_us_fl_ovt = uss_sample$uss_nice$overall$spec_us[i],
        
        fp_uter_us = fp_uter_us[i], 
        fp_loGI_us = fp_loGI_us[i], 
        fp_uter_ovt3 = fp_uter_ovt3[i],       
        fp_loGI_ovt3 = fp_loGI_ovt3[i], 
        fp_lung_us=fp_lung_us,
        fp_panc_us=fp_panc_us,
        fp_lung_ovt3 = fp_lung_ovt3[i], 
        fp_panc_ovt3 = fp_panc_ovt3[i], 
        
        # no uncertainty for service item costs
        cost_us = cost_us,
        cost_us_cdc = cost_us_cdc,
        cost_us_op = cost_us_op,
        cost_ca125 = cost_ca125,
        cost_gp1 = cost_gp1, # face-to-face cons
        cost_gp2 = cost_gp2, # telephone cons after test 
        cost_nurse =cost_nurse,
        cost_ct = cost_ct
      ), 
      accuracy_late = c(
        sens_ca125 = ca125_sample$byAge_stage$sens_ca125_o50_late[i],
        spec_ca125 = ca125_sample$byAge_stage$spec_ca125_o50_late[i],
        sens_ovt1 = ca125_sample$byAge_stage$sens_ovt1_o50_late[i],
        spec_ovt1 = ca125_sample$byAge_stage$spec_ovt1_o50_late[i],
        sens_ovt3 = ca125_sample$byAge_stage$sens_ovt3_o50_late[i],
        spec_ovt3 = ca125_sample$byAge_stage$spec_ovt3_o50_late[i],
        sens_ovt3_fl_us = ca125_sample$byAge_stage$sens_ovt3_o50_late[i],
        spec_ovt3_fl_us = ca125_sample$byAge_stage$spec_ovt3_o50_late[i],
        
        sens_us = uss_sample$uss_nice$overall$sens_us[i],
        spec_us = uss_sample$uss_nice$overall$spec_us[i],
        sens_us_fl_ca125 = uss_sample$uss_nice$overall$sens_us[i],
        spec_us_fl_ca125 = uss_sample$uss_nice$overall$spec_us[i],
        sens_us_fl_ovt = uss_sample$uss_nice$overall$sens_us[i],
        spec_us_fl_ovt = uss_sample$uss_nice$overall$spec_us[i],
        
        fp_uter_us = fp_uter_us[i], 
        fp_loGI_us = fp_loGI_us[i], 
        fp_uter_ovt3 = fp_uter_ovt3[i],       
        fp_loGI_ovt3 = fp_loGI_ovt3[i], 
        fp_lung_us=fp_lung_us,
        fp_panc_us=fp_panc_us,
        fp_lung_ovt3 = fp_lung_ovt3[i], 
        fp_panc_ovt3 = fp_panc_ovt3[i], 
        
        # no uncertainty for service item costs
        cost_us = cost_us,
        cost_us_cdc = cost_us_cdc,
        cost_us_op = cost_us_op,
        cost_ca125 = cost_ca125,
        cost_gp1 = cost_gp1, # face-to-face cons
        cost_gp2 = cost_gp2, # telephone cons after test 
        cost_nurse =cost_nurse,
        cost_ct = cost_ct
      )   
    )
  )
}

acc_u49 <- lapply(accuracy, function(x) x[[1]])
acc_o50 <- lapply(accuracy, function(x) x[[2]])
# 
test_psa3(accuracy = acc_u49)
test_psa3(accuracy = acc_o50)

saveRDS(accuracy, file.path(work_data, "psa", "accuracy_psa_byAgeStage_niceUS.rds"))

##> CA125 varying threshold----

N <- 1001
accuracy <- list()
for (i in 1:N) {
  
  accuracy[[i]] <- list(
    
    # 49 or under
    u49 = list(
      accuracy_early = c(
        sens_ca125 = ca125_sample$byAge$sens_ca125_u49[i],
        spec_ca125 = ca125_sample$byAge$spec_ca125_u49[i],

        sens_us = uss_sample$uss_nice$overall$sens_us[i],
        spec_us = uss_sample$uss_nice$overall$spec_us[i],
        sens_us_fl_ca125 = uss_sample$uss_nice$overall$sens_us[i],
        spec_us_fl_ca125 = uss_sample$uss_nice$overall$spec_us[i],
        sens_us_fl_ovt = uss_sample$uss_nice$overall$sens_us[i],
        spec_us_fl_ovt = uss_sample$uss_nice$overall$spec_us[i],

        # in a new ca125 threshold modified pathway
        # used ca125 accuracy under modified threshold for age groups
        # use the pathway of decision tree 101  
        # so just replace sens and spec with ca125 modified accuracy corresponding to 1% and 3%
        sens_ovt1 = ca125_sample$byAge50$sens_ca125_u49_var[i],
        spec_ovt1 = ca125_sample$byAge50$spec_ca125_u49_var[i],
        
        sens_ovt3 = ca125_sample$byAge50$sens_ca125_u49_var_ovt3[i],
        spec_ovt3 = ca125_sample$byAge50$spec_ca125_u49_var_ovt3[i],
        sens_ovt3_fl_us = ca125_sample$byAge50$sens_ca125_u49_var_ovt3[i],
        spec_ovt3_fl_us = ca125_sample$byAge50$spec_ca125_u49_var_ovt3[i],
        
        fp_uter_us = fp_uter_us[i], 
        fp_loGI_us = fp_loGI_us[i], 
        fp_uter_ovt3 = fp_uter_ovt3[i],       
        fp_loGI_ovt3 = fp_loGI_ovt3[i], 
        fp_lung_us=fp_lung_us,
        fp_panc_us=fp_panc_us,
        fp_lung_ovt3 = fp_lung_ovt3[i], 
        fp_panc_ovt3 = fp_panc_ovt3[i], 
        
        # no uncertainty for service item costs
        cost_us = cost_us,
        cost_us_cdc = cost_us_cdc,
        cost_us_op = cost_us_op,
        cost_ca125 = cost_ca125,
        cost_gp1 = cost_gp1, # face-to-face cons
        cost_gp2 = cost_gp2, # telephone cons after test 
        cost_nurse =cost_nurse,
        cost_ct = cost_ct
      )
    ),
    
    # 51 or over
    o50 = list(
      accuracy_early = c(
        sens_ca125 = ca125_sample$byAge$sens_ca125_o50[i],
        spec_ca125 = ca125_sample$byAge$spec_ca125_o50[i],

        sens_us = uss_sample$uss_nice$overall$sens_us[i],
        spec_us = uss_sample$uss_nice$overall$spec_us[i],
        sens_us_fl_ca125 = uss_sample$uss_nice$overall$sens_us[i],
        spec_us_fl_ca125 = uss_sample$uss_nice$overall$spec_us[i],
        sens_us_fl_ovt = uss_sample$uss_nice$overall$sens_us[i],
        spec_us_fl_ovt = uss_sample$uss_nice$overall$spec_us[i],

        # in a new ca125 threshold modified pathway
        # used ca125 accuracy under modified threshold for age groups
        # use the pathway of decision tree 101  
        # so just replace sens and spec with ca125 modified accuracy corresponding to 1% and 3%
        sens_ovt1 = ca125_sample$byAge50$sens_ca125_o50_var[i],
        spec_ovt1 = ca125_sample$byAge50$spec_ca125_o50_var[i],
        
        sens_ovt3 = ca125_sample$byAge50$sens_ca125_o50_var_ovt3[i],
        spec_ovt3 = ca125_sample$byAge50$spec_ca125_o50_var_ovt3[i],
        sens_ovt3_fl_us = ca125_sample$byAge50$sens_ca125_o50_var_ovt3[i],
        spec_ovt3_fl_us = ca125_sample$byAge50$spec_ca125_o50_var_ovt3[i],
        
        fp_uter_us = fp_uter_us[i], 
        fp_loGI_us = fp_loGI_us[i], 
        fp_uter_ovt3 = fp_uter_ovt3[i],       
        fp_loGI_ovt3 = fp_loGI_ovt3[i], 
        fp_lung_us=fp_lung_us,
        fp_panc_us=fp_panc_us,
        fp_lung_ovt3 = fp_lung_ovt3[i], 
        fp_panc_ovt3 = fp_panc_ovt3[i], 
        
        # no uncertainty for service item costs
        cost_us = cost_us,
        cost_us_cdc = cost_us_cdc,
        cost_us_op = cost_us_op,
        cost_ca125 = cost_ca125,
        cost_gp1 = cost_gp1, # face-to-face cons
        cost_gp2 = cost_gp2, # telephone cons after test 
        cost_nurse =cost_nurse,
        cost_ct = cost_ct
      )
    )
  )
}

for(i in 1:N){
  
  accuracy[[i]]$u49$accuracy_late <- accuracy[[i]]$u49$accuracy_early
  accuracy[[i]]$o50$accuracy_late <- accuracy[[i]]$o50$accuracy_early
  
}

acc_u49 <- lapply(accuracy, function(x) x[[1]])
acc_o50 <- lapply(accuracy, function(x) x[[2]])
# 
test_psa3(accuracy = acc_u49)
test_psa3(accuracy = acc_o50)

saveRDS(accuracy, file.path(work_data, "psa", "accuracy_psa_CA125varyingThreshold_niceUS.rds"))

## Use 2 step US accuracy----
##> By age accuracy for CA125 and US----

N <- 1001
accuracy <- list()
for (i in 1:N) {
  
  accuracy[[i]] <- list(
    
    # 49 or under
    u49 = list(
      accuracy_early = c(
        sens_ca125 = ca125_sample$byAge$sens_ca125_u49[i],
        spec_ca125 = ca125_sample$byAge$spec_ca125_u49[i],
        sens_ovt1 = ca125_sample$byAge$sens_ovt1_u49[i],
        spec_ovt1 = ca125_sample$byAge$spec_ovt1_u49[i],
        sens_ovt3 = ca125_sample$byAge$sens_ovt3_u49[i],
        spec_ovt3 = ca125_sample$byAge$spec_ovt3_u49[i],
        sens_ovt3_fl_us = ca125_sample$byAge$sens_ovt3_u49[i],
        spec_ovt3_fl_us = ca125_sample$byAge$spec_ovt3_u49[i],
        
        sens_us = uss_sample$uss_2step$byAge$sens_us_u49[i],
        spec_us = uss_sample$uss_2step$byAge$spec_us_u49[i],
        sens_us_fl_ca125 = uss_sample$uss_2step$byAge$sens_us_ca125_u49[i],
        spec_us_fl_ca125 = uss_sample$uss_2step$byAge$spec_us_ca125_u49[i],
        sens_us_fl_ovt = uss_sample$uss_2step$byAge$sens_us_ca125_u49[i],
        spec_us_fl_ovt = uss_sample$uss_2step$byAge$spec_us_ca125_u49[i],
        
        fp_uter_us = fp_uter_us[i], 
        fp_loGI_us = fp_loGI_us[i], 
        fp_uter_ovt3 = fp_uter_ovt3[i],       
        fp_loGI_ovt3 = fp_loGI_ovt3[i], 
        fp_lung_us=fp_lung_us,
        fp_panc_us=fp_panc_us,
        fp_lung_ovt3 = fp_lung_ovt3[i], 
        fp_panc_ovt3 = fp_panc_ovt3[i], 
        
        # no uncertainty for service item costs
        cost_us = cost_us,
        cost_us_cdc = cost_us_cdc,
        cost_us_op = cost_us_op,
        cost_ca125 = cost_ca125,
        cost_gp1 = cost_gp1, # face-to-face cons
        cost_gp2 = cost_gp2, # telephone cons after test 
        cost_nurse =cost_nurse,
        cost_ct = cost_ct
      )
    ),
    
    # 51 or over
    o50 = list(
      accuracy_early = c(
        sens_ca125 = ca125_sample$byAge$sens_ca125_o50[i],
        spec_ca125 = ca125_sample$byAge$spec_ca125_o50[i],
        sens_ovt1 = ca125_sample$byAge$sens_ovt1_o50[i],
        spec_ovt1 = ca125_sample$byAge$spec_ovt1_o50[i],
        sens_ovt3 = ca125_sample$byAge$sens_ovt3_o50[i],
        spec_ovt3 = ca125_sample$byAge$spec_ovt3_o50[i],
        sens_ovt3_fl_us = ca125_sample$byAge$sens_ovt3_o50[i],
        spec_ovt3_fl_us = ca125_sample$byAge$spec_ovt3_o50[i],
        
        sens_us = uss_sample$uss_2step$byAge$sens_us_o50[i],
        spec_us = uss_sample$uss_2step$byAge$spec_us_o50[i],
        sens_us_fl_ca125 = uss_sample$uss_2step$byAge$sens_us_ca125_o50[i],
        spec_us_fl_ca125 = uss_sample$uss_2step$byAge$spec_us_ca125_o50[i],
        sens_us_fl_ovt = uss_sample$uss_2step$byAge$sens_us_ca125_o50[i],
        spec_us_fl_ovt = uss_sample$uss_2step$byAge$spec_us_ca125_o50[i],
        
        fp_uter_us = fp_uter_us[i], 
        fp_loGI_us = fp_loGI_us[i], 
        fp_uter_ovt3 = fp_uter_ovt3[i],       
        fp_loGI_ovt3 = fp_loGI_ovt3[i], 
        fp_lung_us=fp_lung_us,
        fp_panc_us=fp_panc_us,
        fp_lung_ovt3 = fp_lung_ovt3[i], 
        fp_panc_ovt3 = fp_panc_ovt3[i], 
        
        # no uncertainty for service item costs
        cost_us = cost_us,
        cost_us_cdc = cost_us_cdc,
        cost_us_op = cost_us_op,
        cost_ca125 = cost_ca125,
        cost_gp1 = cost_gp1, # face-to-face cons
        cost_gp2 = cost_gp2, # telephone cons after test 
        cost_nurse =cost_nurse,
        cost_ct = cost_ct
      )
    )
  )
}

for(i in 1:N){
  
  accuracy[[i]]$u49$accuracy_late <- accuracy[[i]]$u49$accuracy_early
  accuracy[[i]]$o50$accuracy_late <- accuracy[[i]]$o50$accuracy_early
  
}

acc_u49 <- lapply(accuracy, function(x) x[[1]])
acc_o50 <- lapply(accuracy, function(x) x[[2]])
# 
test_psa3(accuracy = acc_u49)
test_psa3(accuracy = acc_o50)

saveRDS(accuracy, file.path(work_data, "psa", "accuracy_psa_byAgeOnly_2stepUS.rds"))

##> By stage for CA125 and overall for US----

N <- 1001
accuracy <- list()
for (i in 1:N) {
  
  accuracy[[i]] <- list(
    
    accuracy_early = c(
      sens_ca125 = ca125_sample$byStage$sens_ca125_early[i],
      spec_ca125 = ca125_sample$byStage$spec_ca125_early[i],
      sens_ovt1 = ca125_sample$byStage$sens_ovt1_early[i],
      spec_ovt1 = ca125_sample$byStage$spec_ovt1_early[i],
      sens_ovt3 = ca125_sample$byStage$sens_ovt3_early[i],
      spec_ovt3 = ca125_sample$byStage$spec_ovt3_early[i],
      sens_ovt3_fl_us = ca125_sample$byStage$sens_ovt3_early[i],
      spec_ovt3_fl_us = ca125_sample$byStage$spec_ovt3_early[i],
      
      sens_us = uss_sample$uss_2step$overall$sens_us[i],
      spec_us = uss_sample$uss_2step$overall$spec_us[i],
      sens_us_fl_ca125 = uss_sample$uss_2step$overall$sens_us_ca125[i],
      spec_us_fl_ca125 = uss_sample$uss_2step$overall$spec_us_ca125[i],
      sens_us_fl_ovt = uss_sample$uss_2step$overall$sens_us_ca125[i],
      spec_us_fl_ovt = uss_sample$uss_2step$overall$spec_us_ca125[i],
      
      
      fp_uter_us = fp_uter_us[i], 
      fp_loGI_us = fp_loGI_us[i], 
      fp_uter_ovt3 = fp_uter_ovt3[i],       
      fp_loGI_ovt3 = fp_loGI_ovt3[i], 
      fp_lung_us=fp_lung_us,
      fp_panc_us=fp_panc_us,
      fp_lung_ovt3 = fp_lung_ovt3[i], 
      fp_panc_ovt3 = fp_panc_ovt3[i], 
      
      # no uncertainty for service item costs
      cost_us = cost_us,
      cost_us_cdc = cost_us_cdc,
      cost_us_op = cost_us_op,
      cost_ca125 = cost_ca125,
      cost_gp1 = cost_gp1, # face-to-face cons
      cost_gp2 = cost_gp2, # telephone cons after test 
      cost_nurse =cost_nurse,
      cost_ct = cost_ct
    ),
    
    accuracy_late = c(
      sens_ca125 = ca125_sample$byStage$sens_ca125_late[i],
      spec_ca125 = ca125_sample$byStage$spec_ca125_late[i],
      sens_ovt1 = ca125_sample$byStage$sens_ovt1_late[i],
      spec_ovt1 = ca125_sample$byStage$spec_ovt1_late[i],
      sens_ovt3 = ca125_sample$byStage$sens_ovt3_late[i],
      spec_ovt3 = ca125_sample$byStage$spec_ovt3_late[i],
      sens_ovt3_fl_us = ca125_sample$byStage$sens_ovt3_late[i],
      spec_ovt3_fl_us = ca125_sample$byStage$spec_ovt3_late[i],
      
      sens_us = uss_sample$uss_2step$overall$sens_us[i],
      spec_us = uss_sample$uss_2step$overall$spec_us[i],
      sens_us_fl_ca125 = uss_sample$uss_2step$overall$sens_us_ca125[i],
      spec_us_fl_ca125 = uss_sample$uss_2step$overall$spec_us_ca125[i],
      sens_us_fl_ovt = uss_sample$uss_2step$overall$sens_us_ca125[i],
      spec_us_fl_ovt = uss_sample$uss_2step$overall$spec_us_ca125[i],
      
      
      fp_uter_us = fp_uter_us[i], 
      fp_loGI_us = fp_loGI_us[i], 
      fp_uter_ovt3 = fp_uter_ovt3[i],       
      fp_loGI_ovt3 = fp_loGI_ovt3[i], 
      fp_lung_us=fp_lung_us,
      fp_panc_us=fp_panc_us,
      fp_lung_ovt3 = fp_lung_ovt3[i], 
      fp_panc_ovt3 = fp_panc_ovt3[i], 
      
      # no uncertainty for service item costs
      cost_us = cost_us,
      cost_us_cdc = cost_us_cdc,
      cost_us_op = cost_us_op,
      cost_ca125 = cost_ca125,
      cost_gp1 = cost_gp1, # face-to-face cons
      cost_gp2 = cost_gp2, # telephone cons after test 
      cost_nurse =cost_nurse,
      cost_ct = cost_ct
    ) 
  )
}
test_psa3(accuracy = accuracy)
saveRDS(accuracy, file.path(work_data, "psa", "accuracy_psa_ca125_byStage_2stepUS_overall.rds"))


##> CA125 varying threshold----
N <- 1001
accuracy <- list()
for (i in 1:N) {
  
  accuracy[[i]] <- list(
    
    # 49 or under
    u49 = list(
      accuracy_early = c(
        sens_ca125 = ca125_sample$byAge$sens_ca125_u49[i],
        spec_ca125 = ca125_sample$byAge$spec_ca125_u49[i],
        
        sens_us = uss_sample$uss_2step$byAge$sens_us_u49[i],
        spec_us = uss_sample$uss_2step$byAge$spec_us_u49[i],
        sens_us_fl_ca125 = uss_sample$uss_2step$byAge$sens_us_ca125_u49[i],
        spec_us_fl_ca125 = uss_sample$uss_2step$byAge$spec_us_ca125_u49[i],
        sens_us_fl_ovt = uss_sample$uss_2step$byAge$sens_us_ca125_u49[i],
        spec_us_fl_ovt = uss_sample$uss_2step$byAge$spec_us_ca125_u49[i],
        
        # in a new ca125 threshold modified pathway
        # used ca125 accuracy under modified threshold for age groups
        # use the pathway of decision tree 101  
        # so just replace sens and spec with ca125 modified accuracy corresponding to 1% and 3%
        sens_ovt1 = ca125_sample$byAge50$sens_ca125_u49_var[i],
        spec_ovt1 = ca125_sample$byAge50$spec_ca125_u49_var[i],
        
        sens_ovt3 = ca125_sample$byAge50$sens_ca125_u49_var_ovt3[i],
        spec_ovt3 = ca125_sample$byAge50$spec_ca125_u49_var_ovt3[i],
        sens_ovt3_fl_us = ca125_sample$byAge50$sens_ca125_u49_var_ovt3[i],
        spec_ovt3_fl_us = ca125_sample$byAge50$spec_ca125_u49_var_ovt3[i],
        
        fp_uter_us = fp_uter_us[i], 
        fp_loGI_us = fp_loGI_us[i], 
        fp_uter_ovt3 = fp_uter_ovt3[i],       
        fp_loGI_ovt3 = fp_loGI_ovt3[i], 
        fp_lung_us=fp_lung_us,
        fp_panc_us=fp_panc_us,
        fp_lung_ovt3 = fp_lung_ovt3[i], 
        fp_panc_ovt3 = fp_panc_ovt3[i], 
        
        # no uncertainty for service item costs
        cost_us = cost_us,
        cost_us_cdc = cost_us_cdc,
        cost_us_op = cost_us_op,
        cost_ca125 = cost_ca125,
        cost_gp1 = cost_gp1, # face-to-face cons
        cost_gp2 = cost_gp2, # telephone cons after test 
        cost_nurse =cost_nurse,
        cost_ct = cost_ct
      )
    ),
    
    # 51 or over
    o50 = list(
      accuracy_early = c(
        sens_ca125 = ca125_sample$byAge$sens_ca125_o50[i],
        spec_ca125 = ca125_sample$byAge$spec_ca125_o50[i],
        
        sens_us = uss_sample$uss_2step$byAge$sens_us_o50[i],
        spec_us = uss_sample$uss_2step$byAge$spec_us_o50[i],
        sens_us_fl_ca125 = uss_sample$uss_2step$byAge$sens_us_ca125_o50[i],
        spec_us_fl_ca125 = uss_sample$uss_2step$byAge$spec_us_ca125_o50[i],
        sens_us_fl_ovt = uss_sample$uss_2step$byAge$sens_us_ca125_o50[i],
        spec_us_fl_ovt = uss_sample$uss_2step$byAge$spec_us_ca125_o50[i],
        
        # in a new ca125 threshold modified pathway
        # used ca125 accuracy under modified threshold for age groups
        # use the pathway of decision tree 101  
        # so just replace sens and spec with ca125 modified accuracy corresponding to 1% and 3%
        sens_ovt1 = ca125_sample$byAge50$sens_ca125_o50_var[i],
        spec_ovt1 = ca125_sample$byAge50$spec_ca125_o50_var[i],
        
        sens_ovt3 = ca125_sample$byAge50$sens_ca125_o50_var_ovt3[i],
        spec_ovt3 = ca125_sample$byAge50$spec_ca125_o50_var_ovt3[i],
        sens_ovt3_fl_us = ca125_sample$byAge50$sens_ca125_o50_var_ovt3[i],
        spec_ovt3_fl_us = ca125_sample$byAge50$spec_ca125_o50_var_ovt3[i],
        
        fp_uter_us = fp_uter_us[i], 
        fp_loGI_us = fp_loGI_us[i], 
        fp_uter_ovt3 = fp_uter_ovt3[i],       
        fp_loGI_ovt3 = fp_loGI_ovt3[i], 
        fp_lung_us=fp_lung_us,
        fp_panc_us=fp_panc_us,
        fp_lung_ovt3 = fp_lung_ovt3[i], 
        fp_panc_ovt3 = fp_panc_ovt3[i], 
        
        # no uncertainty for service item costs
        cost_us = cost_us,
        cost_us_cdc = cost_us_cdc,
        cost_us_op = cost_us_op,
        cost_ca125 = cost_ca125,
        cost_gp1 = cost_gp1, # face-to-face cons
        cost_gp2 = cost_gp2, # telephone cons after test 
        cost_nurse =cost_nurse,
        cost_ct = cost_ct
      )
    )
  )
}

for(i in 1:N){
  
  accuracy[[i]]$u49$accuracy_late <- accuracy[[i]]$u49$accuracy_early
  accuracy[[i]]$o50$accuracy_late <- accuracy[[i]]$o50$accuracy_early
  
}

acc_u49 <- lapply(accuracy, function(x) x[[1]])
acc_o50 <- lapply(accuracy, function(x) x[[2]])
# 
test_psa3(accuracy = acc_u49)
test_psa3(accuracy = acc_o50)

saveRDS(accuracy, file.path(work_data, "psa", "accuracy_psa_CA125varyingThreshold_2stepUS.rds"))
