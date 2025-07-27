rm(list = ls())

source("//Qmcr/dfs/WIPH/Primary Care/OVATOOLS/r_project/scripts/file_path.R")
source(file.path(r_wd, "scripts/fn_MI.R"))

##> MI models, get coef and se----
mod_list <- readRDS(file.path(work_data, "mod_list_MI_stpm2_20240910.rds"))
ova <- print_coef(mod_list, "ova")
loGI <- print_coef(mod_list, "loGI")
lung <- print_coef(mod_list, "lung")
panc <- print_coef(mod_list, "panc")
uter <- print_coef(mod_list, "uter")
other <- print_coef(mod_list, "other")

##> risk model on stage taking average across imputations
mod_list <- readRDS(file.path(work_data, "mod_list_MI_takeMode_stpm2_20240910.rds"))
ova <- print_coef_stpm2(mod_list, "ova")
loGI <- print_coef_stpm2(mod_list, "loGI")
lung <- print_coef_stpm2(mod_list, "lung")
panc <- print_coef_stpm2(mod_list, "panc")
uter <- print_coef_stpm2(mod_list, "uter")
other <- print_coef_stpm2(mod_list, "other")

##> risk model on stage taking average across imputations
mod_list <- readRDS(file.path(output, "mod_cca_stpm2.rds"))
ova <- print_coef_stpm2(mod_list, "ova")
loGI <- print_coef_stpm2(mod_list, "loGI")
lung <- print_coef_stpm2(mod_list, "lung")
panc <- print_coef_stpm2(mod_list, "panc")
uter <- print_coef_stpm2(mod_list, "uter")
other <- print_coef_stpm2(mod_list, "other")

# sort by names----
names_all <- Reduce(union, list(names(ova), names(loGI), names(lung), names(panc),
                                names(uter), names(other)))

names_all[grepl("(Intercept)", names_all)] <- 
  paste0("n1@", names_all[grepl("(Intercept)", names_all)])

names_all[grepl("age_cent60$", names_all)] <- 
  paste0("n2@", names_all[grepl("age_cent60$", names_all)])

names_all[grepl("ethn", names_all)] <- 
  paste0("n3@", names_all[grepl("ethn", names_all)])

names_all[grepl("town", names_all)] <- 
  paste0("n4@", names_all[grepl("town", names_all)])

names_all[grepl("^stage_late$", names_all)] <- 
  paste0("n5@", names_all[grepl("^stage_late$", names_all)])

names_all[grepl("nsx", names_all)] <- 
  paste0("n6@", names_all[grepl("nsx", names_all)])

names_all <- sort(names_all)

names_all <- sapply(strsplit(names_all, "@"), "[", 2)

# assemble----

rt <- cbind(ova[names_all], loGI[names_all], lung[names_all], panc[names_all],
            uter[names_all], other[names_all])

rownames(rt) <- names_all
colnames(rt) <- c("ova", "loGI", "lung", "panc", "uter", "other") 

write.csv(rt, file.path(output, "printed_coef_20240910.csv") )
write.csv(rt, file.path(output, "printed_coef_takeMode20240910.csv") )
write.csv(rt, file.path(output, "printed_coef_cca20240910.csv") )
