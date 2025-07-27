library(tidyverse)
library(haven)
library(foreign)
library(parallel)
library(doSNOW)
library(dtplyr)
library(data.table)
library(survival)

# HRG cleaning ------------------------------------------------------------
# year: starting point
hrg_clean <- function(hes_hrg, year){
  
  # keep useful info
  dt <- hes_hrg %>% 
    select(e_patid, spno, epikey, suscorehrg, sushrg, hes_yr) %>% 
    arrange(e_patid, spno, epikey) %>% 
    mutate(hes_yr=as.integer(hes_yr)) %>% 
    filter(hes_yr>=year) # use at least 1 year before diagnosis for control
  # use 2011 to avoid some admission in 2012 use 2011 HRG
  
  # suscorehrg (core spell HRG) has 5226 ""
  # sushrg (episode HRG) has 5253 ""
  # suscorehrg == "" & sushrg !="" : 0
  # suscorehrg has 242397 "N/A", where all sushrg has meaningful values
  
  # for those non-missing suscorehrg
  # sort by patient id and spell id, remove duplicates
  dt <- dt %>% 
    lazy_dt() %>%
    group_by(e_patid, spno) %>% 
    
    # first, remove duplicated suscorehrg within one spell
    # in effect keep the first row of duplicates
    distinct(suscorehrg, .keep_all = TRUE) %>% 
    ungroup() %>% 
    as.data.frame()
  
  # second, for those suscorehrg has "" or "N/A"
  # if within the same spell, where there are at least one valid suscorehrg value,
  # ignore "" or N/A
  dt <- dt %>% 
    lazy_dt() %>%
    mutate(valid_n = if_else(suscorehrg!="" & suscorehrg!="N/A", 1, 0)) %>% 
    group_by(e_patid, spno) %>%   
    mutate(spell_valid_n = max(valid_n)) %>% 
    ungroup() %>% 
    filter(!((suscorehrg=="" | suscorehrg=="N/A") & spell_valid_n==1)) %>% 
    
    # third, for those remaining suscorehrg=="N/A", use sushrg to impute
    mutate(hrg=if_else(suscorehrg=="N/A", sushrg, suscorehrg)) %>% 
    rename(epatid=e_patid) %>% 
    select(epatid, spno, epikey, hrg, hes_yr) %>% 
    as.data.frame()
  
  # forth, for those remaining suscorehrg=="", all sushrg==""
  # code as NA
  dt$hrg[dt$hrg==""] <- NA
  
  return(dt)
}


# associate with costs ----------------------------------------------------
# input cleaned HRG code data
# study data, for narrow to the ids we are interested in

hrg_cost <- function(hes_hrg_clean){
  
  # read in NHS cost/tariff data
  library(readxl)
  
  # priority given to cost data 21/22 
  # for missing HRG, try through 20/21 to 15/16
  # if there are still missing HRG
  # use tariff data from 21/22 to 14/15
  
  # 1. cost data first
  hrg_cost2122 <- read_excel(file.path(r_wd, "working_data", "NHS_cost_2009_22.xlsx"), sheet = "21-22")
  # some code have double hrg price, keep those with the more activities. 
  hrg_cost2122$Unit_cost <- as.numeric(hrg_cost2122$Unit_cost)
  hrg_cost2122$valid_value <- NULL
  
  # two codes with negative values, i.e. LB60E & RN09Z, considered typos
  # because they are positive in previous files
  # use absolute values
  hrg_cost2122$Unit_cost[hrg_cost2122$hrg=="LB60E" | hrg_cost2122$hrg=="RN09Z"] <- 
    abs(hrg_cost2122$Unit_cost[hrg_cost2122$hrg=="LB60E" | hrg_cost2122$hrg=="RN09Z"])
  
  imp <- hes_hrg_clean %>% 
    # filter(epatid %in% he_ova$epatid) %>% 
    # do not limited in the study population for estimating costs
    left_join(hrg_cost2122, by="hrg") 
  
  imputed <- c()
  
  # check through 09-21
  for (i in 1:13) {
    
    keep <- imp %>% 
      filter(!is.na(Unit_cost))
    
    imputed <- rbind(imputed, keep)
    
    if (i==13) break
    
    hrg_cost <- read_excel(file.path(r_wd, "working_data", "NHS_cost_2009_22.xlsx"), 
                           sheet = paste0(21-i, "-", 22-i)
    )
    hrg_cost$Unit_cost <- as.numeric(hrg_cost$Unit_cost)
    hrg_cost$valid_value <- NULL
    
    miss <- imp %>% 
      filter(is.na(Unit_cost)) %>% 
      select(-Unit_cost, -cost_year) %>% 
      left_join(hrg_cost, by="hrg") 
    
    imp <- miss
    
    final_miss <- miss %>% filter(is.na(Unit_cost))
  }
  
  # 2. tariff data
  hrg_tariff2122 <- read_excel(file.path(r_wd, "working_data", "NHS_tariff_14_22.xlsx"), sheet = "21-22")
  # some code have double hrg price, keep those with the more activities. 
  hrg_tariff2122$Unit_cost <- as.numeric(hrg_tariff2122$Unit_cost)
  hrg_tariff2122$valid_value <- NULL
  
  
  imp <- final_miss %>% 
    select(-Unit_cost, -cost_year) %>% 
    left_join(hrg_tariff2122, by="hrg") 
  
  # check through 14-21
  for (i in 1:10) {
    
    keep <- imp %>% 
      filter(!is.na(Unit_cost))
    
    imputed <- rbind(imputed, keep)
    
    if (i==10) break
    
    if (i < 8) {
      hrg_tariff <- read_excel(file.path(r_wd, "working_data", "NHS_tariff_14_22.xlsx"), 
                               sheet = paste0(21-i, "-", 22-i))
    } else
      if (i==8) {
        hrg_tariff <- read_excel(file.path(r_wd, "working_data", "NHS_tariff_14_22.xlsx"), 
                                 sheet = "10-11")
      }else
        if (i==9)                           
          hrg_tariff <- read_excel(file.path(r_wd, "working_data", "NHS_tariff_14_22.xlsx"), 
                                   sheet = "10-11scots")                                                        
    
    hrg_tariff$Unit_cost <- as.numeric(hrg_tariff$Unit_cost)
    hrg_tariff$valid_value <- NULL
    
    miss <- imp %>% 
      filter(is.na(Unit_cost)) %>% 
      select(-Unit_cost, -cost_year) %>% 
      left_join(hrg_tariff, by="hrg") 
    
    imp <- miss
    
    final_miss2 <- miss %>% filter(is.na(Unit_cost))
  }
  
  # names(table(final_miss2$hrg) )
  # write.csv(table(final_miss2$hrg), file.path(work_data, "final_miss_hrg.csv"))
  
  # 3. manually found the remaining HRG price from different sources
  manual_hrg <- read_excel(file.path(r_wd, "working_data", "final_miss_hrg.xlsx"), sheet = 1)
  
  final_miss2 <- final_miss2 %>%
    select(-Unit_cost, -cost_year) %>%
    left_join(manual_hrg, by="hrg")
  
  # after extensively use national cost 09-22, no HRG codes with missing values after step 1-2  
  # however, after extend scope to 2011, some missingness happen, use 2010-2011 tariff
  # after that 6 HRG codes with missing values, NZ01E-H, NZ03D-E, imputed using the average across its family  
  
  # 4. remaining 1664 HRG codes themselves are missing
  # impute them with total cost/total activity in total HRG 21-22 cost
  imp_avg <- 53481020023/77014922
  
  final_miss2 <- final_miss2 %>% 
    mutate(hrg = if_else(is.na(hrg), "HRG_missing", hrg),
           Unit_cost = if_else(hrg=="HRG_missing", imp_avg, Unit_cost),
           cost_year = if_else(hrg=="HRG_missing", "21-22", cost_year)
    )
  
  rt <- rbind(imputed, final_miss2)
  
  return(rt)
}


# cancer HRG --------------------------------------------------------------
# hes_hrg_cost: cleaned HES HRG data with cost
# hes_admin_diag: HES admission dates
# he_ova: working dataset for HE analysis of ovarian and other cancers
# cancer_type: ova, lung, loGI, uter, panc, benign
# one_year_only: just use diagnosis within 1 year after index date?

cancer_filter <- function(he_ova, cancer_type, one_year_only=F){

  cancer_type_date <- if (cancer_type=="other") "cancer_other" else cancer_type
    
  if (one_year_only){
    
    # 1-year cancer after index date has excluded baseline cancer
    
    if (cancer_type!="benign"){
      pat <- he_ova[he_ova[[paste0("cancer_", cancer_type, "_1yr")]]==1, ]
      pat <- pat[, c("epatid", "censor_date", "death", paste0(cancer_type_date, "_1yr_date"), "index_date", "age")]
    } else {
      # do not restrict to those with a diagnosis of benign gynae
      # pat <- he_ova[he_ova[[paste0(cancer_type, "_1yr")]]==1, ]
      pat <- he_ova[he_ova[["cancer_1yr"]]==0, ]
      pat[["diag_date"]] <- pat[["index_date"]] # for those without cancer, use index date as diagnosis date, as they might not be admitted 
      pat <- pat[, c("epatid", "censor_date", "death", "diag_date", "index_date", "age")]
    }
  } else {
    
    # all cancer after index date may include those with baseline cancer, need to exclude them
    if (cancer_type!="benign") {
      pat <- he_ova[he_ova[[paste0("cancer_", cancer_type, "_icd")]]==1, ] # all incident cancer after index date
      pat <- pat[, c("epatid", "censor_date", "death", paste0(cancer_type_date, "_icd_date"), "index_date", "age")]
    } else {
      # do not restrict to those with a diagnosis of benign gynae
      # include all those not diagnosed with cancer
      # pat <- he_ova[he_ova[[paste0(cancer_type, "_gynae")]]==1, ] # do not exclude baseline benign gynae
      pat <- he_ova[he_ova[["cancer_icd"]]==0, ]
      pat[["diag_date"]] <- pat[["index_date"]] # for those without cancer, use index date as diagnosis date, as they might not be admitted 
      pat <- pat[, c("epatid", "censor_date", "death", "diag_date", "index_date", "age")]
    }
  }
  
  colnames(pat)[4] <- "diag_date"
  
  return(pat)
}

# create a by-year data set by cancer type from diagnosis
cancer_hrg <- function(hes_hrg_cost, hes_admin_diag, he_ova, cancer_type, one_year_only){
  
  hes_admin_diag1 <- hes_admin_diag %>% 
    filter(admidate>="2012-01-01") %>% # 1 year before the starting point
    select(epatid, spno, admidate, discharged)
  
  # incident cancer or benign
  pat <- cancer_filter(he_ova, cancer_type, one_year_only=one_year_only)
  
  rt <- hes_hrg_cost %>% 
    filter(epatid %in% pat$epatid) %>% 
    select(epatid, spno, Unit_cost_inf) %>% 
    group_by(epatid, spno) %>% 
    mutate(cost_n = paste0("cost_", 1:n())) %>% 
    spread(cost_n, Unit_cost_inf) 
  
  # dplyr system is too slow, so use basic R
  rt$cost <- rowSums(rt[grepl("cost_", names(rt))], na.rm = T)
  
  rt <- rt %>% select(-contains("cost_"))

  rt2 <- rt %>% 
    lazy_dt() %>% 
    # mutate(cost = sum(cost_1, cost_2, na.rm = T)) %>% 
    # select(-cost_1, -cost_2) %>% 
    left_join(hes_admin_diag1, by=c("epatid", "spno")) %>% 
    left_join(pat, by="epatid")  %>% 
    filter(admidate>=diag_date - 365.25) %>% # keep 1 year before diagnosis
    mutate(admi_year = case_when(admidate<diag_date ~ 0,
                                 admidate-diag_date<365.25 & admidate>=diag_date  ~ 1,
                                 admidate-diag_date<365.25*2 & admidate-diag_date>=365.25 ~ 2,
                                 admidate-diag_date<365.25*3 & admidate-diag_date>=365.25*2 ~ 3,
                                 admidate-diag_date<365.25*4 & admidate-diag_date>=365.25*3 ~ 4,
                                 admidate-diag_date<365.25*5 & admidate-diag_date>=365.25*4 ~ 5,
                                 admidate-diag_date<365.25*6 & admidate-diag_date>=365.25*5 ~ 6,
                                 admidate-diag_date<365.25*7 & admidate-diag_date>=365.25*6 ~ 7,
                                 admidate-diag_date<365.25*8 & admidate-diag_date>=365.25*7 ~ 8 
    )) %>% 
    mutate(discharge_year = case_when(discharged<diag_date ~ 0,
                                      discharged-diag_date<365.25 & discharged>=diag_date  ~ 1,
                                      discharged-diag_date<365.25*2 & discharged-diag_date>=365.25 ~ 2,
                                      discharged-diag_date<365.25*3 & discharged-diag_date>=365.25*2 ~ 3,
                                      discharged-diag_date<365.25*4 & discharged-diag_date>=365.25*3 ~ 4,
                                      discharged-diag_date<365.25*5 & discharged-diag_date>=365.25*4 ~ 5,
                                      discharged-diag_date<365.25*6 & discharged-diag_date>=365.25*5 ~ 6,
                                      discharged-diag_date<365.25*7 & discharged-diag_date>=365.25*6 ~ 7,
                                      discharged-diag_date<365.25*8 & discharged-diag_date>=365.25*7 ~ 8 
    )) %>% 
    # max(discharge_year - admi_year)=1
    # therefore if discharge and admission is not within the same year
    # only calculate the ratio in the first year
    mutate(spell_days = as.numeric(discharged - admidate) + 1, 
           admi_year_ratio = if_else(admi_year==discharge_year, 1, 
                                     # for admission before diagnosis and discharge after diagnosis, assume all cost in discharge year (year 1)
                                     if_else(admi_year!=discharge_year & admi_year==0, 0, 
                                             as.numeric((diag_date + admi_year*365.25) - admidate)/spell_days))
    ) %>% 
    select(epatid, cost, admi_year, admi_year_ratio) %>% 
    as.data.frame()
  
  # split costs across two years since diagnosis
  
  # costs assigned to the admission year
  rt21 <- rt2 %>% 
    mutate(cost_admi_year = cost * admi_year_ratio) %>% 
    filter(cost_admi_year !=0) %>% 
    select(-cost, -admi_year_ratio)
  
  # costs assigned to next year
  rt22 <- rt2 %>% 
    filter(admi_year_ratio<1) %>% 
    mutate(cost_admi_year = cost * (1-admi_year_ratio)) %>% 
    mutate(admi_year = admi_year + 1) %>% 
    select(-cost, -admi_year_ratio)
  
  # combine costs in the same year  
  rt3 <- rt21 %>% 
    bind_rows(rt22) %>% 
    group_by(epatid, admi_year) %>% 
    summarise(cost = sum(cost_admi_year, na.rm = T))
  
  # generate by year records
  
  # according to the HES document, the censoring date is 2021-03-31
  # the maximum hrg year is 2020
  # however, actually, the hrg code corresponding admission date can be as late as 2021-03-31
  # continue to use 2021-03-31 as the censoring date.
  pat <- pat %>% 
    mutate(fu_time = as.numeric(censor_date - diag_date)) %>% 
    mutate(fu_time = if_else(fu_time==0, 0.5, fu_time)) %>% 
    filter(fu_time > 0) # rid of negative fu time
  
  surv <- as.formula(Surv(fu_time, death) ~ .)
  
  fu_length <- ceiling(max(pat$fu_time)/365.25)
  
  dt <- survSplit(surv, pat, cut = 365.25*(1:fu_length), episode = "admi_year")
  
  # generate proportion of follow up in one year
  # insert admi_year = 0, 1 year before diagnosis
  admi0 <- cbind(pat[, c("epatid", "censor_date", "diag_date", "index_date", "age")], 
                 tstart=-365.25, fu_time=0, death=0, admi_year=0, fu_prop=1)
  
  dt <- dt %>% 
    mutate(fu_prop = (fu_time-tstart)/365.25) %>% 
    bind_rows(admi0) %>% 
    arrange(epatid, admi_year) %>% 
    mutate(fu_prop = 1 - fu_prop) %>%  # no lost fu = 0, as baseline
    # age at the beginning of the year
    mutate(curr_age = round(age + as.numeric(diag_date - index_date)/365.25 + admi_year - 1)) %>% 
    relocate(curr_age, .after=age)
  
  dt_cost <- dt %>% 
    full_join(rt3, join_by("epatid", "admi_year"))
  
  # two ids produced HRG long after death, remove them
  dt_cost <- dt_cost %>% 
    filter(!is.na(death)) %>% # because records from rt3 only do not contain death info
    mutate(cost = if_else(is.na(cost), 0 , cost)) # NA = no HRG at that year
  
  return(dt_cost)
}

# # for those without any HRG codes, generate cost by years until censor_date as 0
# # for benign only
# # the above code has included those without HRG codes and fill them as 0!!!
### no need the following function anymore

# fill_no_hrg <- function(hes_hrg_cost, he_ova, cancer_type, one_year_only){
#   
#   # incident cancer or benign
#   pat <- cancer_filter(he_ova, cancer_type, one_year_only=one_year_only)
#   
#   pat <- pat %>% 
#     filter(!epatid %in% hes_hrg_cost$epatid)  
#     
#   # according to the HES document, the censoring date is 2021-03-31
#   # the maximum hrg year is 2020
#   # however, actually, the hrg code corresponding admission date can be as late as 2021-03-31
#   # continue to use 2021-03-31 as the censoring date.
#   pat <- pat %>% 
#     mutate(fu_time = as.numeric(censor_date - diag_date)) %>% 
#     mutate(fu_time = if_else(fu_time==0, 0.5, fu_time)) %>% 
#     filter(fu_time > 0) # rid of negative fu time
#   
#   surv <- as.formula(Surv(fu_time, death) ~ .)
#   
#   fu_length <- ceiling(max(pat$fu_time)/365.25)
#   
#   dt <- survSplit(surv, pat, cut = 365.25*(1:fu_length), episode = "admi_year")
#     
#   # generate proportion of follow up in one year
#   # insert admi_year = 0, 1 year before diagnosis
#   admi0 <- cbind(pat[, c("epatid", "censor_date", "diag_date", "index_date", "age")], 
#                  tstart=-365.25, fu_time=0, death=0, admi_year=0, fu_prop=1)
#   
#   dt <- dt %>% 
#     mutate(fu_prop = (fu_time-tstart)/365.25) %>% 
#     bind_rows(admi0) %>% 
#     arrange(epatid, admi_year) %>% 
#     mutate(fu_prop = 1 - fu_prop) %>%  # no lost fu = 0, as baseline
#     # age at the beginning of the year
#     mutate(curr_age = round(age + as.numeric(diag_date - index_date)/365.25 + admi_year - 1)) %>% 
#     relocate(curr_age, .after=age) 
#     
#   dt$cost <- 0
#   
#   return(dt)
#  
# }
# 

# create a by-year data set for all from index date, including those w/o cancer
all_hrg <- function(hes_hrg_cost, hes_admin_diag, he_ova){
  
  hes_admin_diag1 <- hes_admin_diag %>% 
    filter(admidate>="2012-01-01") %>% # 1 year before the starting point
    select(epatid, spno, admidate, discharged)
  
  ids <- he_ova %>% 
    select(epatid, index_date)
  
  rt <- hes_hrg_cost %>%
    filter(epatid %in% ids$epatid) %>% 
    select(epatid, spno, Unit_cost_inf) %>% 
    group_by(epatid, spno) %>% 
    mutate(cost_n = paste0("cost_", 1:n())) %>% 
    spread(cost_n, Unit_cost_inf)
  
  rt2 <- rt %>% 
    lazy_dt() %>% 
    mutate(cost = sum(cost_1, cost_2, cost_3, na.rm = T)) %>% 
    select(-cost_1, -cost_2, -cost_3) %>% 
    left_join(hes_admin_diag1, by=c("epatid", "spno")) %>% 
    left_join(ids, by="epatid") %>%
    filter(admidate>=index_date - 365.25) %>% # keep 1 year before index date
    mutate(admi_year = case_when(admidate<index_date ~ 0,
                                 admidate-index_date<365.25 & admidate>=index_date  ~ 1,
                                 admidate-index_date<365.25*2 & admidate-index_date>=365.25 ~ 2,
                                 admidate-index_date<365.25*3 & admidate-index_date>=365.25*2 ~ 3,
                                 admidate-index_date<365.25*4 & admidate-index_date>=365.25*3 ~ 4,
                                 admidate-index_date<365.25*5 & admidate-index_date>=365.25*4 ~ 5,
                                 admidate-index_date<365.25*6 & admidate-index_date>=365.25*5 ~ 6,
                                 admidate-index_date<365.25*7 & admidate-index_date>=365.25*6 ~ 7,
                                 admidate-index_date<365.25*8 & admidate-index_date>=365.25*7 ~ 8 
    )) %>% 
    mutate(discharge_year = case_when(discharged<index_date ~ 0,
                                      discharged-index_date<365.25 & discharged>=index_date  ~ 1,
                                      discharged-index_date<365.25*2 & discharged-index_date>=365.25 ~ 2,
                                      discharged-index_date<365.25*3 & discharged-index_date>=365.25*2 ~ 3,
                                      discharged-index_date<365.25*4 & discharged-index_date>=365.25*3 ~ 4,
                                      discharged-index_date<365.25*5 & discharged-index_date>=365.25*4 ~ 5,
                                      discharged-index_date<365.25*6 & discharged-index_date>=365.25*5 ~ 6,
                                      discharged-index_date<365.25*7 & discharged-index_date>=365.25*6 ~ 7,
                                      discharged-index_date<365.25*8 & discharged-index_date>=365.25*7 ~ 8 
    )) %>% 
    mutate(spell_days = as.numeric(discharged - admidate) + 1, 
           # calculate the ratio in the first year 
           admi_year_ratio = if_else(admi_year==discharge_year, 1, 
                                     as.numeric((index_date + admi_year*365.25) - admidate)/spell_days)
    ) %>% 
    as.data.frame() 
  
  # split costs across years since index
  # max(discharge_year- admi_year)=4. but only a few >1
  
  # costs assigned to the admission year
  rt20 <- rt2 %>% 
    mutate(cost_at_year = cost * admi_year_ratio) %>% 
    filter(cost_at_year !=0) %>% 
    mutate(cost_year = admi_year) %>% 
    select(epatid, cost_at_year, cost_year)
  
  # costs assigned to next year
  rt21 <- rt2 %>% 
    filter(discharge_year - admi_year==1) %>% 
    mutate(cost_at_year = cost * (1-admi_year_ratio)) %>% 
    mutate(cost_year = admi_year + 1) %>% 
    select(epatid, cost_at_year, cost_year)
  
  rt211 <- rt2 %>% 
    filter(discharge_year - admi_year>1) %>% 
    # because this year as a whole in hospital if discharge_year- admi_year>1
    mutate(cost_at_year = cost * 365.25/spell_days) %>% 
    mutate(cost_year = admi_year + 1) %>% 
    select(epatid, cost_at_year, cost_year)
  
  # costs assigned to 3rd year 
  rt22 <- rt2 %>% 
    filter(discharge_year - admi_year==2) %>% 
    mutate(cost_at_year = cost * (1 - admi_year_ratio - 365.25/spell_days)) %>% 
    mutate(cost_year = admi_year + 2) %>% 
    select(epatid, cost_at_year, cost_year)
  
  rt221 <- rt2 %>% 
    filter(discharge_year - admi_year>2) %>% 
    # because this year as a whole in hospital if discharge_year- admi_year>2
    mutate(cost_at_year = cost * 365.25/spell_days) %>% 
    mutate(cost_year = admi_year + 2) %>% 
    select(epatid, cost_at_year, cost_year)
  
  # costs assinged to 4th year
  rt23 <- rt2 %>% 
    filter(discharge_year - admi_year==3) %>% 
    mutate(cost_at_year = cost * (1 - admi_year_ratio - 2*365.25/spell_days)) %>% 
    mutate(cost_year = admi_year + 3) %>% 
    select(epatid, cost_at_year, cost_year)
  
  rt231 <- rt2 %>% 
    filter(discharge_year - admi_year>3) %>% 
    # because this year as a whole in hospital if discharge_year- admi_year>2
    mutate(cost_at_year = cost * 365.25/spell_days) %>% 
    mutate(cost_year = admi_year + 3) %>% 
    select(epatid, cost_at_year, cost_year)
  
  # costs assinged to 5th year
  rt24 <- rt2 %>% 
    filter(discharge_year - admi_year==4) %>% 
    mutate(cost_at_year = cost * (1 - admi_year_ratio - 3*365.25/spell_days)) %>% 
    mutate(cost_year = admi_year + 4) %>% 
    select(epatid, cost_at_year, cost_year)
  
  # combine costs in the same year  
  rt3 <- rt20 %>% 
    bind_rows(rt21, rt211, rt22, rt221, rt23, rt231, rt24) %>% 
    group_by(epatid, cost_year) %>% 
    summarise(cost = sum(cost_at_year, na.rm = T))
  
  
  pat <- he_ova %>% 
    lazy_dt() %>% 
    select(epatid, index_date, age,  
           cancer_icd_date, cancer_icd, 
           ova_icd_date, cancer_ova_icd, 
           panc_icd_date, cancer_panc_icd, 
           lung_icd_date, cancer_lung_icd, 
           uter_icd_date, cancer_uter_icd,
           loGI_icd_date, cancer_loGI_icd, 
           benign_date, benign_gynae, # for benign incident, do not exclude benign history
           b_cancer_nonOva, b_benign_gynae, b_panc, b_lung, b_uter, b_loGI, 
           censor_date, death
    ) %>% 
    mutate(fu_time = as.numeric(censor_date - index_date)) %>% 
    mutate(fu_time = if_else(fu_time==0, 0.5, fu_time)) %>% 
    filter(fu_time > 0) %>%  # rid of censor date earlier than index date
    as.data.frame()
  
  # split by year from index date until censoring
  surv <- as.formula(Surv(fu_time, death) ~ .)
  
  fu_length <- ceiling(max(pat$fu_time)/365.25)
  
  dt <- survSplit(surv, pat, cut = 365.25*(1:fu_length), episode = "cost_year")
  
  # generate proportion of follow up in one year
  # insert cost_year = 0, 1 year before diagnosis
  admi0 <- cbind(pat, tstart=-365.25, cost_year=0, fu_prop=1)
  admi0$death <- 0
  admi0$fu_time <- 0
  
  dt <- dt %>% 
    mutate(fu_prop = (fu_time-tstart)/365.25) %>% 
    bind_rows(admi0) %>% 
    arrange(epatid, cost_year) %>% 
    mutate(fu_prop = 1 - fu_prop) # no lost fu = 0, as baseline
  
  dt <- dt %>% 
    # a few ids produced HRG long after death, left_join remove them
    left_join(rt3, join_by("epatid", "cost_year")) %>% 
    mutate(cost = if_else(is.na(cost), 0 , cost)) # NA = no HRG at that year
  
  # create cancer events history by year 
  dt$cir.start <- dt$index_date + dt$tstart
  
  events <- c("ova", "lung", "panc", "uter", "loGI")
  for (e in events) {
    
    dt[[paste0(e, "_0_1")]] <- dt[[paste0(e, "_1_2")]] <- dt[[paste0(e, "_2_3")]] <- 
      dt[[paste0(e, "_3_4")]] <- dt[[paste0(e, "_4_inf")]] <- 0
    
    c_ind <- paste0("cancer_", e, "_icd")
    c_date <- paste0(e, "_icd_date")
    
    dt[[paste0(e, "_0_1")]][dt[[c_ind]]==1 & 
                              dt[[c_date]] > dt$cir.start & 
                              dt[[c_date]] < dt$cir.start + 365.25] <- 1
    
    dt[[paste0(e, "_1_2")]][dt[[c_ind]]==1 & 
                              dt[[c_date]] + 365.25 > dt$cir.start & 
                              dt[[c_date]] + 365.25 <= dt$cir.start + 365.25] <- 1
    
    dt[[paste0(e, "_2_3")]][dt[[c_ind]]==1 & 
                              dt[[c_date]] + 365.25*2 > dt$cir.start & 
                              dt[[c_date]] + 365.25*2 <= dt$cir.start + 365.25] <- 1
    
    dt[[paste0(e, "_3_4")]][dt[[c_ind]]==1 & 
                              dt[[c_date]] + 365.25*3 > dt$cir.start & 
                              dt[[c_date]] + 365.25*3 <= dt$cir.start + 365.25] <- 1
    
    dt[[paste0(e, "_4_inf")]][dt[[c_ind]]==1 & 
                                dt[[c_date]] + 365.25*4 <= dt$cir.start + 365.25] <- 1
  }
  
  # for benign gynae, only code the same year event and other history
  dt$benign_0_1 <- dt$benign_1_inf <- 0
  dt$benign_0_1[dt$benign_gynae==1 & 
                  dt$benign_date > dt$cir.start & 
                  dt$benign_date < dt$cir.start + 365.25] <- 1
  
  dt$benign_1_inf[dt$benign_gynae==1 & 
                    dt$benign_date <= dt$cir.start] <- 1
  
  dt$curr_age <- (dt$age + dt$cost_year - 1)
  
  return(dt)
}  



