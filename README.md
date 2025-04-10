# Ovatools
The health economic model codes for evaluating Ovatools

If you want to use the same directory structure as that in the Github, before running, change the general working directory address to your local address in file_path.R

The model code files are all saved in the scripts folder. 

# The decision tree model
The health economic model consists of a decision tree model and a Markov model. To use the decision tree model, first open the psa_decision_tree.R, which is used for probabilistic sensitivity analysis (PSA). When using it for deterministic analysis, set i = 1001. 

psa_decision_tree.R calls fn_master_decision_tree.R, which is the master function to run the decision tree. In it, 
I = 1001: means using the deterministic parameters
byAge = TRUE: use the by-age accuracy for the CA125 test. The default. 
accuracy_file: the name of the input accuracy file. Currently, we only provide the “byAgeOnly_niceUS” file. It has separate CA125 accuracy by age 50. You could make your accuracy file using the same structure of this file. 
adj_pop = TRUE: use the CA125 tested population data. It is the default, and currently we only provide this population data.
cost_us_cdc = FALSE: do not use the ultrasound cost in UK community diagnostic centre, which is much cheaper than the average ultrasound cost.
cost_us_op = FALSE: do not use the ultrasound cost in UK outpatient departments, which is slightly more expensive than the average.
cost_us_input = NA: You can input a certain amount of the ultrasound cost.
us_sens_adj = 0 and us_spec_adj = 0: the number used to modify the ultrasound sensitivity and specificity (ranging from 0 - 1). The modifier can be positive or negative. 

fn_master_decision_tree.R calls fn_decision_tree.R, where three types of primary care diagnostic pathways are modelled. The models run in probabilities, depending on cancer incidence and accuracy data. The main outcomes of the model are the probabilities in true positive (TP), false positive (FP), true negative (TN) and false negative (FN), for the detection result and ovarian cancer diagnosis. 
decision_tree100: the current practice. Patients are tested CA125 first. If CA125 ≥ 35 U/mL, an ultrasound scan is followed, and referral to hospital if ultrasound is abnormal.  
decision_tree101: the risk-based triage pathway. Patients are tested CA125 first, and ovarian cancer risk is calculated based CA125 values and age. If the moderate-risk threshold is reached (ovt1), an ultrasound scan is followed, and referral to hospital if ultrasound is abnormal. If the high-risk is reached (ovt3), referral to hospital. 

decision_tree110 and decision_tree112: patients are tested CA125 and ultrasound concurrently. Patients are referred to hospital, if CA125 or ultrasound result is abnormal. The threshold for CA125 is ≥ 35 U/mL in 110, and high risk in 112. A correlation coefficient rho is included to adjust the correlation between CA125 and ultrasound accuracy. The default is 0, which means no correlation. 

In alternative pathways, patients are potentially referred to hospital without an abnormal ultrasound result. The model estimates the abnormal ultrasound probability in a repeated ultrasound in hospital for patients referred to hospital without an abnormal ultrasound in primary care and without cancer diagnosed eventually. The probability could be used to estimate the surgery probability in the FP cases, assuming only patients with abnormal ultrasound would receive surgery. 

The model estimates proportions of other cancers, including uterine (uter), lower gastrointestinal (loGI), lung and pancreatic (panc) cancers, in FP and TN, assuming that no one are diagnosed with ovarian cancer and another cancer.
 
The model estimates the cost in the diagnostic pathway in primary care.

# The Markov model
