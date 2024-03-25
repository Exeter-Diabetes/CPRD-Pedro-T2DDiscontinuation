####################
## Description: 
##  - In this file we:
##    - Fit the already deployed SGLT2 vs GLP1 model
#################### 



# load working directory
setwd("Samples/T2D_Discontinuation")

# load functions
source("code/00.set_up_data.R")

# load libraries
library(tidyverse)
library(bcf)

###############################################################################
###############################################################################
############################### Read Data In ##################################
###############################################################################
###############################################################################

# load dataset
cprd_dataset <- set_up_data(
  raw_data = "20240308_t2d_1stinstance",
  diagnosis = FALSE,
  therapies = c("GLP1", "SGLT2"),
  dataset = "full.dataset"
) %>%
  drop_na(stopdrug_3m_6mFU)



cprd_dataset_model <- cprd_dataset %>%
  select(
    patid,
    dstartdate,
    drugclass,
    gender,
    dstartdate_age,
    prebmi,
    prehba1c,
    preegfr,
    dstartdate_dm_dur,
    prealt,
    predrug_pad,
    predrug_heartfailure,
    predrug_ihd,
    predrug_neuropathy,
    predrug_retinopathy,
    drugline,
    numdrugs
  ) %>%
  mutate(
    drugclass = factor(drugclass, levels = c("GLP1", "SGLT2")),
    gender = factor(gender, levels = c(2, 1), labels = c("Female", "Male")),
    predrug_pad = factor(predrug_pad, levels = c(0, 1), labels = c("No", "Yes")),
    predrug_heartfailure = factor(predrug_heartfailure, levels = c(0, 1), labels = c("No", "Yes")),
    predrug_ihd = factor(predrug_ihd, levels = c(0, 1), labels = c("No", "Yes")),
    predrug_neuropathy = factor(predrug_neuropathy, levels = c(0, 1), labels = c("No", "Yes")),
    predrug_retinopathy = factor(predrug_retinopathy, levels = c(0, 1), labels = c("No", "Yes")),
    hba1cmonth = as.numeric(12)
  ) %>%
  rename(
    "sex" = "gender",
    "agetx" = "dstartdate_age",
    "t2dmduration" = "dstartdate_dm_dur",
    "prepad" = "predrug_pad",
    "preheartfailure" = "predrug_heartfailure",
    "preihd" = "predrug_ihd",
    "preneuropathy" = "predrug_neuropathy",
    "preretinopathy" = "predrug_retinopathy",
    "ncurrtx" = "numdrugs"
  ) %>%
  as.data.frame() %>%
  drop_na()



bcf_model <- readRDS("code/Models/SGLT2vsGLP1/bcf_model.rds")


prediction_dataset <- cprd_dataset_model %>%
  mutate(drugclass = "SGLT2") %>%
  rbind(
    cprd_dataset_model %>%
      mutate(drugclass = "GLP1")
  )



predictions <-  predict(object = bcf_model,
                        x_predict_control = prediction_dataset %>%
                          select(
                            c("agetx", "t2dmduration", "drugline", "ncurrtx", "hba1cmonth", "prehba1c", "preegfr", "prealt", "prepad")
                          ) %>%
                          mutate_all(list(~as.numeric(.))) %>%
                          as.matrix(),
                        x_predict_moderate = prediction_dataset %>%
                          select(
                            c("agetx", "sex", "ncurrtx", "prehba1c", "prebmi", "preegfr", "preheartfailure" ,"preihd", "preneuropathy", "prepad" ,"preretinopathy")
                          ) %>%
                          mutate_all(list(~as.numeric(.))) %>%
                          as.matrix(),
                        pi_pred = rep(0.5, nrow(prediction_dataset)),
                        z_pred = prediction_dataset %>%
                          select(drugclass) %>%
                          mutate(drugclass = ifelse(drugclass == "GLP1", 0, 1)) %>%
                          unlist(),
                        save_tree_directory = "code/Models/SGLT2vsGLP1/", 
                        verbose = FALSE,
                        log_file = NULL)



model_predictions_dataset <- data.frame(
  patid = cprd_dataset_model$patid,
  dstartdate = cprd_dataset_model$dstartdate,
  pred.SGLT2 = colMeans(predictions$yhat)[1:nrow(cprd_dataset_model)],
  pred.GLP1 = colMeans(predictions$yhat)[1:nrow(cprd_dataset_model) + nrow(cprd_dataset_model)],
  pred.tau = colMeans(predictions$yhat)[1:nrow(cprd_dataset_model)] - colMeans(predictions$yhat)[1:nrow(cprd_dataset_model)+nrow(cprd_dataset_model)]
) %>% 
  as.data.frame()


saveRDS(model_predictions_dataset, "results/Models/SGLT2vsGLP1/model_predictions_dataset.rds")



##################################################

## Tables comparing above 5 mmol/mol benefit for those who discontinue or not

library(tableone)

# load dataset
cprd_dataset <- set_up_data(
  raw_data = "20240308_t2d_1stinstance",
  diagnosis = FALSE,
  therapies = c("DPP4", "GLP1", "MFN", "SGLT2", "SU", "TZD"),
  dataset = "full.dataset"
) %>%
  left_join(
    readRDS("results/Models/SGLT2vsGLP1/model_predictions_dataset.rds") %>%
      select(
        patid, dstartdate, pred.tau
      ), by = c("patid", "dstartdate")
  ) %>%
  mutate(
    benefit = ifelse(pred.tau > 5, "GLP-1RA >5", ifelse(pred.tau < 5 & pred.tau > 0, "GLP-1RA <5", ifelse(pred.tau < 0 & pred.tau > -5, "SGLT2i <5", ifelse(pred.tau < -5, "SGLT2i >5", "Test"))))
    # benefit = ifelse(pred.tau > 5, "GLP-1RA", ifelse(pred.tau < -5, "SGLT2i", NA))
  ) %>%
  drop_na(benefit)




output_path <- "results/tables"


cprd_tables <- cprd_dataset %>%
  mutate(
    prefastingglucose_na = ifelse(!is.na(prefastingglucose), "No", NA),
    prehdl_na = ifelse(!is.na(prehdl), "No", NA),
    pretriglyceride_na = ifelse(!is.na(pretriglyceride), "No", NA),
    precreatinine_blood_na = ifelse(!is.na(precreatinine_blood), "No", NA),
    preldl_na = ifelse(!is.na(preldl), "No", NA),
    prealt_na = ifelse(!is.na(prealt), "No", NA),
    preast_na = ifelse(!is.na(preast), "No", NA),
    pretotalcholesterol_na = ifelse(!is.na(pretotalcholesterol), "No", NA),
    predbp_na = ifelse(!is.na(predbp), "No", NA),
    presbp_na = ifelse(!is.na(presbp), "No", NA),
    preacr_na = ifelse(!is.na(preacr), "No", NA),
    prehba1c_na = ifelse(!is.na(prehba1c), "No", NA),
    preegfr_na = ifelse(!is.na(preegfr), "No", NA),
    prealbumin_blood_na = ifelse(!is.na(prealbumin_blood), "No", NA),
    prebilirubin_na = ifelse(!is.na(prebilirubin), "No", NA),
    prehaematocrit_na = ifelse(!is.na(prehaematocrit), "No", NA),
    prehaemoglobin_na = ifelse(!is.na(prehaemoglobin), "No", NA),
    prepcr_na = ifelse(!is.na(prepcr), "No", NA),
    dm_diag_age_na = ifelse(!is.na(dm_diag_age), "No", NA),
    prebmi_na = ifelse(!is.na(prebmi), "No", NA),
    dstartdate_age_na = ifelse(!is.na(dstartdate_age), "No", NA),
    dstartdate_dm_dur_na = ifelse(!is.na(dstartdate_dm_dur), "No", NA),
    dstartmonth_na = ifelse(!is.na(dstartmonth), "No", NA)
  )



vars <- c(
  # Outcome
  "stopdrug_3m_6mFU", "stopdrug_6m_6mFU",
  "pred.tau", "benefit",
  # Drug taken
  "drugline", "numdrugs", " drugclass",
  # Biomarkers
  "prefastingglucose", "prefastingglucose_na", "prehdl", "prehdl_na", 
  "pretriglyceride", "pretriglyceride_na", "precreatinine_blood", "precreatinine_blood_na", 
  "preldl", "preldl_na", "prealt", "prealt_na", "preast", "preast_na",
  "pretotalcholesterol", "pretotalcholesterol_na", "predbp", "predbp_na", 
  "presbp", "presbp_na", "preacr", "preacr_na", "prehba1c", "prehba1c_na", 
  "preegfr", "preegfr_na", "prealbumin_blood", "prealbumin_blood_na", 
  "prebilirubin", "prebilirubin_na", "prehaematocrit", "prehaematocrit_na",
  "prehaemoglobin", "prehaemoglobin_na", "prepcr", "prepcr_na", 
  # Commorbidities
  "preckdstage", "predrug_frailty_mild", "predrug_frailty_moderate", 
  "predrug_frailty_severe", "predrug_primary_hhf", "predrug_af", "predrug_angina", 
  "predrug_asthma", "predrug_bronchiectasis", "predrug_cld", "predrug_copd", 
  "predrug_cysticfibrosis", "predrug_dementia", "predrug_diabeticnephropathy", 
  "predrug_fh_premature_cvd", "predrug_haem_cancer", "predrug_heartfailure", 
  "predrug_hypertension", "predrug_ihd", "predrug_myocardialinfarction", 
  "predrug_neuropathy", "predrug_otherneuroconditions", "predrug_pad", 
  "predrug_pulmonaryfibrosis", "predrug_pulmonaryhypertension", 
  "predrug_retinopathy", "predrug_revasc", "predrug_rheumatoidarthritis", 
  "predrug_solid_cancer", "predrug_solidorgantransplant", "predrug_stroke", 
  "predrug_tia", "predrug_anxiety_disorders", "predrug_medspecific_gi",
  "predrug_benignprostatehyperplasia", "predrug_micturition_control",
  "predrug_volume_depletion", "predrug_urinary_frequency", "predrug_falls",
  "predrug_lowerlimbfracture", "predrug_incident_mi", "predrug_incident_stroke",
  "predrug_dka", "predrug_osteoporosis", "predrug_unstableangina", 
  "predrug_amputation",
  "hosp_admission_prev_year", "hosp_admission_prev_year_count",
  # Extra info
  "gender", "prac_region", "ethnicity_5cat", "imd2015_10", "dm_diag_age", "dm_diag_age_na",
  "ins_in_1_year", "prebmi", "prebmi_na", "smoking_cat", "stopdrug_3m_3mFU_MFN_hist",
  "alcohol_cat", "fh_diabetes", "dstartdate_age", "dstartdate_age_na", 
  "dstartdate_dm_dur", "dstartdate_dm_dur_na", "dstartmonth", "dstartmonth_na", "CCI_index"
)



vars_cat <- c(
  # Outcome
  "stopdrug_3m_6mFU", "stopdrug_6m_6mFU",
  "benefit",
  # Drug taken
  "drugline", "numdrugs", "drugclass",
  # Biomarkers
  "prefastingglucose_na", "prehdl_na", 
  "pretriglyceride_na", "precreatinine_blood_na", 
  "preldl_na", "prealt_na", "preast_na",
  "pretotalcholesterol_na", "predbp_na", 
  "presbp_na", "preacr_na", "prehba1c_na", 
  "preegfr_na", "prealbumin_blood_na", 
  "prebilirubin_na", "prehaematocrit_na",
  "prehaemoglobin_na", "prepcr_na",
  # Commorbidities
  "preckdstage", "predrug_frailty_mild", "predrug_frailty_moderate", 
  "predrug_frailty_severe", "predrug_primary_hhf", "predrug_af", "predrug_angina", 
  "predrug_asthma", "predrug_bronchiectasis", "predrug_cld", "predrug_copd", 
  "predrug_cysticfibrosis", "predrug_dementia", "predrug_diabeticnephropathy", 
  "predrug_fh_premature_cvd", "predrug_haem_cancer", "predrug_heartfailure", 
  "predrug_hypertension", "predrug_ihd", "predrug_myocardialinfarction", 
  "predrug_neuropathy", "predrug_otherneuroconditions", "predrug_pad", 
  "predrug_pulmonaryfibrosis", "predrug_pulmonaryhypertension", 
  "predrug_retinopathy", "predrug_revasc", "predrug_rheumatoidarthritis", 
  "predrug_solid_cancer", "predrug_solidorgantransplant", "predrug_stroke", 
  "predrug_tia", "predrug_anxiety_disorders", "predrug_medspecific_gi",
  "predrug_benignprostatehyperplasia", "predrug_micturition_control",
  "predrug_volume_depletion", "predrug_urinary_frequency", "predrug_falls",
  "predrug_lowerlimbfracture", "predrug_incident_mi", "predrug_incident_stroke",
  "predrug_dka", "predrug_osteoporosis", "predrug_unstableangina", 
  "predrug_amputation",
  "hosp_admission_prev_year",
  # Extra info
  "gender", "prac_region", "ethnicity_5cat", "imd2015_10", "dm_diag_age_na",
  "ins_in_1_year", "prebmi_na", "smoking_cat", "stopdrug_3m_3mFU_MFN_hist",
  "alcohol_cat", "fh_diabetes", "dstartdate_age_na", 
  "dstartdate_dm_dur_na", "dstartmonth_na"
)


### DPP4

table_characteristics_benefit <- CreateTableOne(
  vars = vars,
  factorVars = vars_cat,
  includeNA = TRUE,
  strata = c("benefit", "drugclass"),
  data = cprd_tables,
  test = TRUE,
  smd = TRUE
)

table_characteristics_benefit_print <- print(table_characteristics_benefit, exact = "stage", quote = FALSE, noSpaces = TRUE, printToggle = FALSE)

write.csv(table_characteristics_benefit_print, file = paste0(output_path, "/table_characteristics_sglt2_vs_glp1_5m_benefit.csv"))







