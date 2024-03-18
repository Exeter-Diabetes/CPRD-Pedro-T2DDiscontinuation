####################
## Description:   (this was run with bartMachine v1.3.4.1)
##  - In this file we:
##    - Fit a BART PS model to all variables.
##    - Perform variable selection on PS model and refit model.
##    - Refit BART PS model
#################### 


## increase memory usage to 100GB of RAM (needs to be run before library(bartMachine))
options(java.parameters = "-Xmx100g")

# load working directory
setwd("Samples/T2D_Discontinuation")

# load functions
source("code/00.set_up_data.R")

## Load libraries
library(tidyverse)
library(bartMachine)

bartMachine::set_bart_machine_num_cores(1)


###############################################################################
###############################################################################
############################### Read Data In ##################################
###############################################################################
###############################################################################

# load dataset
cprd_dataset <- set_up_data(
  raw_data = "20240308_t2d_1stinstance",
  diagnosis = FALSE,
  therapies = c("DPP4", "GLP1", "MFN", "SGLT2", "SU", "TZD"),
  dataset = "ps.dataset"
)


########
### Fit a propensity score model to all variables in dataset
###   Predicting whether they go on drug a vs others
########

### MFN

# Predicting  drug assignment

bart_ps_model_MFN <- bartMachine::bartMachine(X = cprd %>%
                                                select(
                                                  # Biomarkers
                                                  precreatinine_blood,
                                                  prealt, pretotalcholesterol, predbp, presbp, prehba1c,
                                                  preegfr, prebilirubin,
                                                  # # Commorbidities
                                                  preckdstage, 
                                                  # predrug_frailty_mild,
                                                  predrug_primary_hhf, predrug_af, predrug_angina,
                                                  # predrug_asthma, predrug_bronchiectasis, predrug_cld, predrug_copd,
                                                  # predrug_cysticfibrosis,
                                                  # predrug_dementia,
                                                  predrug_diabeticnephropathy,
                                                  # predrug_fh_premature_cvd,
                                                  # predrug_haem_cancer,
                                                  predrug_heartfailure,
                                                  predrug_hypertension, predrug_ihd, predrug_myocardialinfarction,
                                                  predrug_neuropathy, predrug_otherneuroconditions, predrug_pad,
                                                  # predrug_pulmonaryfibrosis,
                                                  # predrug_pulmonaryhypertension,
                                                  predrug_retinopathy, predrug_revasc,
                                                  # predrug_rheumatoidarthritis,
                                                  # predrug_solid_cancer,
                                                  # predrug_solidorgantransplant,
                                                  predrug_stroke,
                                                  predrug_tia,
                                                  # predrug_anxiety_disorders, predrug_medspecific_gi,
                                                  # predrug_benignprostatehyperplasia, predrug_micturition_control,
                                                  # predrug_volume_depletion, predrug_urinary_frequency,
                                                  predrug_falls,
                                                  # predrug_lowerlimbfracture, predrug_incident_mi, predrug_incident_stroke,
                                                  predrug_dka, predrug_osteoporosis, predrug_unstableangina,
                                                  # predrug_amputation,
                                                  hosp_admission_prev_year,
                                                  # hosp_admission_prev_year_count,
                                                  # # Extra info
                                                  gender, 
                                                  # prac_region, 
                                                  ethnicity_5cat, imd2015_10, dm_diag_age,
                                                  ins_in_1_year, prebmi, smoking_cat, drugline, stopdrug_3m_3mFU_MFN_hist,
                                                  # alcohol_cat, 
                                                  dstartdate_age, dstartdate_dm_dur, 
                                                  # dstartmonth,
                                                  CCI_index
                                                ),
                                              y = cprd %>%
                                                mutate(drugclass = ifelse(drugclass == "MFN", "MFN", "Other"),
                                                       drugclass = factor(drugclass)) %>%
                                                select(drugclass) %>%
                                                unlist(),
                                              num_trees = 50,
                                              use_missing_data = TRUE,
                                              num_burn_in = 2500,
                                              num_iterations_after_burn_in = 2000,
                                              serialize = TRUE)


saveRDS(bart_ps_model_MFN, "PS_model/bart_ps_model_MFN.rds")


# Variable selection

pdf(file = "figures/Prop Score Var Selection/bart_ps_model_var_select_MFN.pdf", width = 18, height = 11)
# error with cv
vs_bart_ps_model_MFN <- var_selection_by_permute(bart_ps_model_MFN)
dev.off()

saveRDS(vs_bart_ps_model_MFN, "PS_model/vs_bart_ps_model_MFN.rds")


### DPP4

# Predicting  drug assignment

bart_ps_model_DPP4 <- bartMachine::bartMachine(X = cprd %>%
                                                 select(
                                                   # Biomarkers
                                                   precreatinine_blood,
                                                   prealt, pretotalcholesterol, predbp, presbp, prehba1c,
                                                   preegfr, prebilirubin,
                                                   # # Commorbidities
                                                   preckdstage, 
                                                   # predrug_frailty_mild,
                                                   predrug_primary_hhf, predrug_af, predrug_angina,
                                                   # predrug_asthma, predrug_bronchiectasis, predrug_cld, predrug_copd,
                                                   # predrug_cysticfibrosis,
                                                   # predrug_dementia,
                                                   predrug_diabeticnephropathy,
                                                   # predrug_fh_premature_cvd,
                                                   # predrug_haem_cancer,
                                                   predrug_heartfailure,
                                                   predrug_hypertension, predrug_ihd, predrug_myocardialinfarction,
                                                   predrug_neuropathy, predrug_otherneuroconditions, predrug_pad,
                                                   # predrug_pulmonaryfibrosis,
                                                   # predrug_pulmonaryhypertension,
                                                   predrug_retinopathy, predrug_revasc,
                                                   # predrug_rheumatoidarthritis,
                                                   # predrug_solid_cancer,
                                                   # predrug_solidorgantransplant,
                                                   predrug_stroke,
                                                   predrug_tia,
                                                   # predrug_anxiety_disorders, predrug_medspecific_gi,
                                                   # predrug_benignprostatehyperplasia, predrug_micturition_control,
                                                   # predrug_volume_depletion, predrug_urinary_frequency,
                                                   predrug_falls,
                                                   # predrug_lowerlimbfracture, predrug_incident_mi, predrug_incident_stroke,
                                                   predrug_dka, predrug_osteoporosis, predrug_unstableangina,
                                                   # predrug_amputation,
                                                   hosp_admission_prev_year,
                                                   # hosp_admission_prev_year_count,
                                                   # # Extra info
                                                   gender, 
                                                   # prac_region, 
                                                   ethnicity_5cat, imd2015_10, dm_diag_age,
                                                   ins_in_1_year, prebmi, smoking_cat, drugline, stopdrug_3m_3mFU_MFN_hist,
                                                   # alcohol_cat, 
                                                   dstartdate_age, dstartdate_dm_dur, 
                                                   # dstartmonth,
                                                   CCI_index
                                                 ),
                                               y = cprd %>%
                                                 mutate(drugclass = ifelse(drugclass == "DPP4", "DPP4", "Other"),
                                                        drugclass = factor(drugclass)) %>%
                                                 select(drugclass) %>%
                                                 unlist(),
                                               num_trees = 50,
                                               use_missing_data = TRUE,
                                               num_burn_in = 2500,
                                               num_iterations_after_burn_in = 2000,
                                               serialize = TRUE)



saveRDS(bart_ps_model_DPP4, "PS_model/bart_ps_model_DPP4.rds")


# Variable selection

pdf(file = "figures/Prop Score Var Selection/bart_ps_model_var_select_DPP4.pdf", width = 18, height = 11)
# error with cv
vs_bart_ps_model_DPP4 <- var_selection_by_permute(bart_ps_model_DPP4)
dev.off()

saveRDS(vs_bart_ps_model_DPP4, "PS_model/vs_bart_ps_model_DPP4.rds")



### GLP1

# Predicting  drug assignment

bart_ps_model_GLP1 <- bartMachine::bartMachine(X = cprd %>%
                                                 select(
                                                   # Biomarkers
                                                   precreatinine_blood,
                                                   prealt, pretotalcholesterol, predbp, presbp, prehba1c,
                                                   preegfr, prebilirubin,
                                                   # # Commorbidities
                                                   preckdstage, 
                                                   # predrug_frailty_mild,
                                                   predrug_primary_hhf, predrug_af, predrug_angina,
                                                   # predrug_asthma, predrug_bronchiectasis, predrug_cld, predrug_copd,
                                                   # predrug_cysticfibrosis,
                                                   # predrug_dementia,
                                                   predrug_diabeticnephropathy,
                                                   # predrug_fh_premature_cvd,
                                                   # predrug_haem_cancer,
                                                   predrug_heartfailure,
                                                   predrug_hypertension, predrug_ihd, predrug_myocardialinfarction,
                                                   predrug_neuropathy, predrug_otherneuroconditions, predrug_pad,
                                                   # predrug_pulmonaryfibrosis,
                                                   # predrug_pulmonaryhypertension,
                                                   predrug_retinopathy, predrug_revasc,
                                                   # predrug_rheumatoidarthritis,
                                                   # predrug_solid_cancer,
                                                   # predrug_solidorgantransplant,
                                                   predrug_stroke,
                                                   predrug_tia,
                                                   # predrug_anxiety_disorders, predrug_medspecific_gi,
                                                   # predrug_benignprostatehyperplasia, predrug_micturition_control,
                                                   # predrug_volume_depletion, predrug_urinary_frequency,
                                                   predrug_falls,
                                                   # predrug_lowerlimbfracture, predrug_incident_mi, predrug_incident_stroke,
                                                   predrug_dka, predrug_osteoporosis, predrug_unstableangina,
                                                   # predrug_amputation,
                                                   hosp_admission_prev_year,
                                                   # hosp_admission_prev_year_count,
                                                   # # Extra info
                                                   gender, 
                                                   # prac_region, 
                                                   ethnicity_5cat, imd2015_10, dm_diag_age,
                                                   ins_in_1_year, prebmi, smoking_cat, drugline, stopdrug_3m_3mFU_MFN_hist,
                                                   # alcohol_cat, 
                                                   dstartdate_age, dstartdate_dm_dur, 
                                                   # dstartmonth,
                                                   CCI_index
                                                 ),
                                               y = cprd %>%
                                                 mutate(drugclass = ifelse(drugclass == "GLP1", "GLP1", "Other"),
                                                        drugclass = factor(drugclass)) %>%
                                                 select(drugclass) %>%
                                                 unlist(),
                                               num_trees = 50,
                                               use_missing_data = TRUE,
                                               num_burn_in = 2500,
                                               num_iterations_after_burn_in = 2000,
                                               serialize = TRUE)



saveRDS(bart_ps_model_GLP1, "PS_model/bart_ps_model_GLP1.rds")


# Variable selection

pdf(file = "figures/Prop Score Var Selection/bart_ps_model_var_select_GLP1.pdf", width = 18, height = 11)
# error with cv
vs_bart_ps_model_GLP1 <- var_selection_by_permute(bart_ps_model_GLP1)
dev.off()

saveRDS(vs_bart_ps_model_GLP1, "PS_model/vs_bart_ps_model_GLP1.rds")



### SGLT2

# Predicting  drug assignment

bart_ps_model_SGLT2 <- bartMachine::bartMachine(X = cprd %>%
                                                  select(
                                                    # Biomarkers
                                                    precreatinine_blood,
                                                    prealt, pretotalcholesterol, predbp, presbp, prehba1c,
                                                    preegfr, prebilirubin,
                                                    # # Commorbidities
                                                    preckdstage, 
                                                    # predrug_frailty_mild,
                                                    predrug_primary_hhf, predrug_af, predrug_angina,
                                                    # predrug_asthma, predrug_bronchiectasis, predrug_cld, predrug_copd,
                                                    # predrug_cysticfibrosis,
                                                    # predrug_dementia,
                                                    predrug_diabeticnephropathy,
                                                    # predrug_fh_premature_cvd,
                                                    # predrug_haem_cancer,
                                                    predrug_heartfailure,
                                                    predrug_hypertension, predrug_ihd, predrug_myocardialinfarction,
                                                    predrug_neuropathy, predrug_otherneuroconditions, predrug_pad,
                                                    # predrug_pulmonaryfibrosis,
                                                    # predrug_pulmonaryhypertension,
                                                    predrug_retinopathy, predrug_revasc,
                                                    # predrug_rheumatoidarthritis,
                                                    # predrug_solid_cancer,
                                                    # predrug_solidorgantransplant,
                                                    predrug_stroke,
                                                    predrug_tia,
                                                    # predrug_anxiety_disorders, predrug_medspecific_gi,
                                                    # predrug_benignprostatehyperplasia, predrug_micturition_control,
                                                    # predrug_volume_depletion, predrug_urinary_frequency,
                                                    predrug_falls,
                                                    # predrug_lowerlimbfracture, predrug_incident_mi, predrug_incident_stroke,
                                                    predrug_dka, predrug_osteoporosis, predrug_unstableangina,
                                                    # predrug_amputation,
                                                    hosp_admission_prev_year,
                                                    # hosp_admission_prev_year_count,
                                                    # # Extra info
                                                    gender, 
                                                    # prac_region, 
                                                    ethnicity_5cat, imd2015_10, dm_diag_age,
                                                    ins_in_1_year, prebmi, smoking_cat, drugline, stopdrug_3m_3mFU_MFN_hist,
                                                    # alcohol_cat, 
                                                    dstartdate_age, dstartdate_dm_dur, 
                                                    # dstartmonth,
                                                    CCI_index
                                                  ),
                                                y = cprd %>%
                                                  mutate(drugclass = ifelse(drugclass == "SGLT2", "SGLT2", "Other"),
                                                         drugclass = factor(drugclass)) %>%
                                                  select(drugclass) %>%
                                                  unlist(),
                                                num_trees = 50,
                                                use_missing_data = TRUE,
                                                num_burn_in = 2500,
                                                num_iterations_after_burn_in = 2000,
                                                serialize = TRUE)



saveRDS(bart_ps_model_SGLT2, "PS_model/bart_ps_model_SGLT2.rds")


# Variable selection

pdf(file = "figures/Prop Score Var Selection/bart_ps_model_var_select_SGLT2.pdf", width = 18, height = 11)
# error with cv
vs_bart_ps_model_SGLT2 <- var_selection_by_permute(bart_ps_model_SGLT2)
dev.off()

saveRDS(vs_bart_ps_model_SGLT2, "PS_model/vs_bart_ps_model_SGLT2.rds")



### SU

# Predicting  drug assignment

bart_ps_model_SU <- bartMachine::bartMachine(X = cprd %>%
                                               select(
                                                 # Biomarkers
                                                 precreatinine_blood,
                                                 prealt, pretotalcholesterol, predbp, presbp, prehba1c,
                                                 preegfr, prebilirubin,
                                                 # # Commorbidities
                                                 preckdstage, 
                                                 # predrug_frailty_mild,
                                                 predrug_primary_hhf, predrug_af, predrug_angina,
                                                 # predrug_asthma, predrug_bronchiectasis, predrug_cld, predrug_copd,
                                                 # predrug_cysticfibrosis,
                                                 # predrug_dementia,
                                                 predrug_diabeticnephropathy,
                                                 # predrug_fh_premature_cvd,
                                                 # predrug_haem_cancer,
                                                 predrug_heartfailure,
                                                 predrug_hypertension, predrug_ihd, predrug_myocardialinfarction,
                                                 predrug_neuropathy, predrug_otherneuroconditions, predrug_pad,
                                                 # predrug_pulmonaryfibrosis,
                                                 # predrug_pulmonaryhypertension,
                                                 predrug_retinopathy, predrug_revasc,
                                                 # predrug_rheumatoidarthritis,
                                                 # predrug_solid_cancer,
                                                 # predrug_solidorgantransplant,
                                                 predrug_stroke,
                                                 predrug_tia,
                                                 # predrug_anxiety_disorders, predrug_medspecific_gi,
                                                 # predrug_benignprostatehyperplasia, predrug_micturition_control,
                                                 # predrug_volume_depletion, predrug_urinary_frequency,
                                                 predrug_falls,
                                                 # predrug_lowerlimbfracture, predrug_incident_mi, predrug_incident_stroke,
                                                 predrug_dka, predrug_osteoporosis, predrug_unstableangina,
                                                 # predrug_amputation,
                                                 hosp_admission_prev_year,
                                                 # hosp_admission_prev_year_count,
                                                 # # Extra info
                                                 gender, 
                                                 # prac_region, 
                                                 ethnicity_5cat, imd2015_10, dm_diag_age,
                                                 ins_in_1_year, prebmi, smoking_cat, drugline, stopdrug_3m_3mFU_MFN_hist,
                                                 # alcohol_cat, 
                                                 dstartdate_age, dstartdate_dm_dur, 
                                                 # dstartmonth,
                                                 CCI_index
                                               ),
                                             y = cprd %>%
                                               mutate(drugclass = ifelse(drugclass == "SU", "SU", "Other"),
                                                      drugclass = factor(drugclass)) %>%
                                               select(drugclass) %>%
                                               unlist(),
                                             num_trees = 50,
                                             use_missing_data = TRUE,
                                             num_burn_in = 2500,
                                             num_iterations_after_burn_in = 2000,
                                             serialize = TRUE)


saveRDS(bart_ps_model_SU, "PS_model/bart_ps_model_SU.rds")


# Variable selection

pdf(file = "figures/Prop Score Var Selection/bart_ps_model_var_select_SU.pdf", width = 18, height = 11)
# error with cv
vs_bart_ps_model_SU <- var_selection_by_permute(bart_ps_model_SU)
dev.off()

saveRDS(vs_bart_ps_model_SU, "PS_model/vs_bart_ps_model_SU.rds")




### TZD

# Predicting  drug assignment

bart_ps_model_TZD <- bartMachine::bartMachine(X = cprd %>%
                                                select(
                                                  # Biomarkers
                                                  precreatinine_blood,
                                                  prealt, pretotalcholesterol, predbp, presbp, prehba1c,
                                                  preegfr, prebilirubin,
                                                  # # Commorbidities
                                                  preckdstage, 
                                                  # predrug_frailty_mild,
                                                  predrug_primary_hhf, predrug_af, predrug_angina,
                                                  # predrug_asthma, predrug_bronchiectasis, predrug_cld, predrug_copd,
                                                  # predrug_cysticfibrosis,
                                                  # predrug_dementia,
                                                  predrug_diabeticnephropathy,
                                                  # predrug_fh_premature_cvd,
                                                  # predrug_haem_cancer,
                                                  predrug_heartfailure,
                                                  predrug_hypertension, predrug_ihd, predrug_myocardialinfarction,
                                                  predrug_neuropathy, predrug_otherneuroconditions, predrug_pad,
                                                  # predrug_pulmonaryfibrosis,
                                                  # predrug_pulmonaryhypertension,
                                                  predrug_retinopathy, predrug_revasc,
                                                  # predrug_rheumatoidarthritis,
                                                  # predrug_solid_cancer,
                                                  # predrug_solidorgantransplant,
                                                  predrug_stroke,
                                                  predrug_tia,
                                                  # predrug_anxiety_disorders, predrug_medspecific_gi,
                                                  # predrug_benignprostatehyperplasia, predrug_micturition_control,
                                                  # predrug_volume_depletion, predrug_urinary_frequency,
                                                  predrug_falls,
                                                  # predrug_lowerlimbfracture, predrug_incident_mi, predrug_incident_stroke,
                                                  predrug_dka, predrug_osteoporosis, predrug_unstableangina,
                                                  # predrug_amputation,
                                                  hosp_admission_prev_year,
                                                  # hosp_admission_prev_year_count,
                                                  # # Extra info
                                                  gender, 
                                                  # prac_region, 
                                                  ethnicity_5cat, imd2015_10, dm_diag_age,
                                                  ins_in_1_year, prebmi, smoking_cat, drugline, stopdrug_3m_3mFU_MFN_hist,
                                                  # alcohol_cat, 
                                                  dstartdate_age, dstartdate_dm_dur, 
                                                  # dstartmonth,
                                                  CCI_index
                                                ),
                                              y = cprd %>%
                                                mutate(drugclass = ifelse(drugclass == "TZD", "TZD", "Other"),
                                                       drugclass = factor(drugclass)) %>%
                                                select(drugclass) %>%
                                                unlist(),
                                              num_trees = 50,
                                              use_missing_data = TRUE,
                                              num_burn_in = 2500,
                                              num_iterations_after_burn_in = 2000,
                                              serialize = TRUE)



saveRDS(bart_ps_model_TZD, "PS_model/bart_ps_model_TZD.rds")


# Variable selection

pdf(file = "figures/Prop Score Var Selection/bart_ps_model_var_select_TZD.pdf", width = 18, height = 11)
# error with cv
vs_bart_ps_model_TZD <- var_selection_by_permute(bart_ps_model_TZD)
dev.off()

saveRDS(vs_bart_ps_model_TZD, "PS_model/vs_bart_ps_model_TZD.rds")


