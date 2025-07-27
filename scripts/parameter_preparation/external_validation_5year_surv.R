rm(list = ls())

library(tidyverse)
source("//Qmcr/dfs/WIPH/Primary Care/OVATOOLS/r_project/scripts/file_path.R")
source(file.path(post, "fn_post_analysis.R"))

use_mort <- TRUE
ipd2 <- readRDS(file.path(work_data, "ipd2_20241126.rds"))
# ipd2_20241126.rds and ipd2_20240906.RDS only have different var "group"
# as the size and order are not changed, Markov model results are unchanged from ipd2_20240906.RDS are unchanged
# because Markov model results does not include variable "group"
id_dic3 <- readRDS(file.path(work_data, "id_dic3_20241126.rds"))
ipd <- readRDS(file.path(work_data, "ipd_new20241126.RDS")) 
ipd_new <- readRDS(file.path(work_data, "ipd_new20240906.RDS")) 

c_type <- "ova"
stage_late <- 1
output_name <- paste0("5years_", c_type, "_stage", stage_late,"_useMort_", as.character(use_mort))
ova_late <- read.csv(file.path(output, paste0(output_name,"_sum.csv")))

c_type <- "ova"
stage_late <- 0
output_name <- paste0("5years_", c_type, "_stage", stage_late,"_useMort_", as.character(use_mort))
ova_early <- read.csv(file.path(output, paste0(output_name,"_sum.csv")))

c_type <- "lung"
stage_late <- 1
output_name <- paste0("5years_", c_type, "_stage", stage_late,"_useMort_", as.character(use_mort))
lung_late <- read.csv(file.path(output, paste0(output_name,"_sum.csv")))

c_type <- "lung"
stage_late <- 0
output_name <- paste0("5years_", c_type, "_stage", stage_late,"_useMort_", as.character(use_mort))
lung_early <- read.csv(file.path(output, paste0(output_name,"_sum.csv")))

c_type <- "loGI"
stage_late <- 1
output_name <- paste0("5years_", c_type, "_stage", stage_late,"_useMort_", as.character(use_mort))
loGI_late <- read.csv(file.path(output, paste0(output_name,"_sum.csv")))

c_type <- "loGI"
stage_late <- 0
output_name <- paste0("5years_", c_type, "_stage", stage_late,"_useMort_", as.character(use_mort))
loGI_early <- read.csv(file.path(output, paste0(output_name,"_sum.csv")))

c_type <- "panc"
stage_late <- 1
output_name <- paste0("5years_", c_type, "_stage", stage_late,"_useMort_", as.character(use_mort))
panc_late <- read.csv(file.path(output, paste0(output_name,"_sum.csv")))

c_type <- "panc"
stage_late <- 0
output_name <- paste0("5years_", c_type, "_stage", stage_late,"_useMort_", as.character(use_mort))
panc_early <- read.csv(file.path(output, paste0(output_name,"_sum.csv")))

c_type <- "uter"
stage_late <- 1
output_name <- paste0("5years_", c_type, "_stage", stage_late,"_useMort_", as.character(use_mort))
uter_late <- read.csv(file.path(output, paste0(output_name,"_sum.csv")))

c_type <- "uter"
stage_late <- 0
output_name <- paste0("5years_", c_type, "_stage", stage_late,"_useMort_", as.character(use_mort))
uter_early <- read.csv(file.path(output, paste0(output_name,"_sum.csv")))


rst_list <- c("ova_early", "ova_late", "lung_early", "lung_late", 
              "panc_early", "panc_late", "uter_early", "uter_late", 
              "loGI_early", "loGI_late")

for (rst in rst_list) {
  to_update <- get(rst)
  updated <- map2pop(to_update, id_dic3, ipd2)
  updated <- updated %>% select(-group, -epatid)
  assign(rst, updated)
}

fn_5ys <- function(rst){
  
  rt <- 1-(mean(rst$cd)+mean(rst$ncd))
  
  return(rt)
}

fn_5ys(ova_late)*0.67+fn_5ys(ova_early)*0.33

fn_5ys(lung_late)*0.76+fn_5ys(lung_early)*0.24

fn_5ys(panc_late)*0.85+fn_5ys(panc_early)*0.15

fn_5ys(loGI_late)*0.63+fn_5ys(loGI_early)*0.37

fn_5ys(uter_late)*0.28+fn_5ys(uter_early)*0.72




