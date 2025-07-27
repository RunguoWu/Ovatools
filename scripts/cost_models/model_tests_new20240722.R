#################################################################
##                       cost model explore                    ##
#################################################################
rm(list = ls())

source("//Qmcr/dfs/WIPH/Primary Care/OVATOOLS/r_project/scripts/file_path.R")

source(file.path(r_wd, "scripts/cost/fn_model_tests.R"))

he_ova <- readRDS(file.path(work_data, "he_ova_step4_withSymptom_prescr_NCRASupdated.RDS"))

# Ova ---------------------------------------------------------------------
cancer_type <- "ova"

# prepare the analytical data for test
# model_1p_2p will prepare the data automatically
d4s <- model_prep(cancer_type, he_ova)

# variables in the cost model
var_x <- c("ethn3", "town5", "curr_age_cent", "stage_late", "admi_year_fct",
           "death_c", "death_nc", "fu_prop")

model_test(cancer_type, he_ova, var_x)

# loGI --------------------------------------------------------------------
cancer_type <- "loGI"

# prepare the analytical data for test
# model_1p_2p will prepare the data automatically
d4s <- model_prep(cancer_type, he_ova)

# variables in the cost model
var_x <- c("ethn3", "town5", "curr_age_cent", "stage_late", "admi_year_fct",
           "death_c", "death_nc", "fu_prop")

model_test(cancer_type, he_ova, var_x)

# panc --------------------------------------------------------------------
cancer_type <- "panc"

# prepare the analytical data for test
# model_1p_2p will prepare the data automatically
d4s <- model_prep(cancer_type, he_ova)

# variables in the cost model
var_x <- c("ethn3", "town5", "curr_age_cent", "stage_late", "admi_year_fct",
           "death_c", "death_nc", "fu_prop")

model_test(cancer_type, he_ova, var_x)

# uter --------------------------------------------------------------------
cancer_type <- "uter"

# prepare the analytical data for test
# model_1p_2p will prepare the data automatically
d4s <- model_prep(cancer_type, he_ova)

# variables in the cost model
var_x <- c("ethn3", "town5", "curr_age_cent", "stage_late", "admi_year_fct",
           "death_c", "death_nc", "fu_prop")

model_test(cancer_type, he_ova, var_x)

# lung --------------------------------------------------------------------
cancer_type <- "lung"

# prepare the analytical data for test
# model_1p_2p will prepare the data automatically
d4s <- model_prep(cancer_type, he_ova)

# variables in the cost model
var_x <- c("ethn3", "town5", "curr_age_cent", "stage_late", "admi_year_fct",
           "death_c", "death_nc", "fu_prop")

model_test(cancer_type, he_ova, var_x)

# other -------------------------------------------------------------------
cancer_type <- "other"

# prepare the analytical data for test
# model_1p_2p will prepare the data automatically
d4s <- model_prep(cancer_type, he_ova)

# variables in the cost model
var_x <- c("ethn3", "town5", "curr_age_cent", "stage_late", "admi_year_fct",
           "death_c", "death_nc", "fu_prop")

model_test(cancer_type, he_ova, var_x)


# benign ------------------------------------------------------------------
cancer_type <- "benign"

d4s <- model_prep(cancer_type, he_ova)

# variables in the cost model
var_x <- c("ethn3", "town5", "curr_age_cent", "admi_year_fct",
           "death_c", "death_nc", "fu_prop")

model_test(cancer_type, he_ova, var_x)

# the process is very slow; show the result below: gam_log
# [[1]]
# # A tibble: 7 × 7
# mod      mp_slope hl_p  pl_p     me   mae  rmse
# <chr>    <chr>    <chr> <chr> <dbl> <dbl> <dbl>
#   1 gau_id   1.86     <0.01 <0.01     0  2547  4527
# 2 gau_log  1.9      <0.01 0.97      0  2544  4525
# 3 poi_id   1.95     <0.01 <0.01     0  2546  4527
# 4 poi_log  1.93     <0.01 <0.01     0  2545  4525
# 5 gam_id   1.97     <0.01 <0.01     2  2545  4527
# 6 gam_log  1.95     <0.01 <0.01     1  2544  4525
# 7 gam_sqrt 1.96     <0.01 <0.01     2  2545  4526
# 
# [[2]]
# # A tibble: 8 × 4
# mod               me   mae  rmse
# <chr>          <dbl> <dbl> <dbl>
#   1 1p_gau_id          0  1181  2648
# 2 2p_p2_gau_id       0  1175  2644
# 3 2p_p2_gau_log      0  1174  2644
# 4 2p_p2_poi_id       0  1175  2644
# 5 2p_p2_poi_log      0  1175  2644
# 6 2p_p2_gam_id      -1  1175  2645
# 7 2p_p2_gam_log      0  1175  2644
# 8 2p_p2_gam_sqrt    -1  1175  2644









