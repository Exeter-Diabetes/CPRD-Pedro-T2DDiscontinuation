####################
## Description: 
##  - In this file we:
##    - Fit a linear PS model to all variables.
#################### 



# load working directory
setwd("Samples/T2D_Discontinuation")

# load functions
source("code/00.set_up_data.R")

# load libraries
library(MASS)
library(tidyverse)
library(PSweight)

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
### Features for a propensity score
########


drug.pscores <- SumStat(
  ps.formula = as.formula(paste0("drugclass ~ ", paste0(c(
    # Biomarkers
    "precreatinine_blood",
    "prealt", "pretotalcholesterol", "predbp", "presbp", "prehba1c",
    "preegfr", "prebilirubin",
    # # Commorbidities
    "preckdstage", 
    # "predrug_frailty_mild",
    "predrug_primary_hhf", "predrug_af", "predrug_angina",
    # "predrug_asthma", "predrug_bronchiectasis", "predrug_cld", "predrug_copd",
    # "predrug_cysticfibrosis",
    # "predrug_dementia",
    "predrug_diabeticnephropathy",
    # "predrug_fh_premature_cvd",
    # "predrug_haem_cancer",
    "predrug_heartfailure",
    "predrug_hypertension", "predrug_ihd", "predrug_myocardialinfarction",
    "predrug_neuropathy", "predrug_otherneuroconditions", "predrug_pad",
    # "predrug_pulmonaryfibrosis",
    # "predrug_pulmonaryhypertension",
    "predrug_retinopathy", "predrug_revasc",
    # "predrug_rheumatoidarthritis",
    # "predrug_solid_cancer",
    # "predrug_solidorgantransplant",
    "predrug_stroke",
    "predrug_tia",
    # "predrug_anxiety_disorders", "predrug_medspecific_gi",
    # "predrug_benignprostatehyperplasia", "predrug_micturition_control",
    # "predrug_volume_depletion", "predrug_urinary_frequency",
    "predrug_falls",
    # "predrug_lowerlimbfracture", "predrug_incident_mi", "predrug_incident_stroke",
    "predrug_dka", "predrug_osteoporosis", "predrug_unstableangina",
    # "predrug_amputation",
    "hosp_admission_prev_year",
    # "hosp_admission_prev_year_count",
    # # Extra info
    "gender", 
    # "prac_region", 
    "ethnicity_5cat", "imd2015_10", "dm_diag_age",
    "ins_in_1_year", "prebmi", "smoking_cat", "drugline", "stopdrug_3m_3mFU_MFN_hist",
    # "alcohol_cat", 
    "dstartdate_age", "dstartdate_dm_dur", 
    # "dstartmonth",
    "CCI_index"
  ), collapse = "+"))),
  data = cprd_dataset,
  weight = c("overlap")
)

saveRDS(drug.pscores, "results/PS_model/drug.pscores_model.rds")


##########################################
#  Save Propensity scores
#

ps.only_dataset <- data.frame(
  patid = cprd_dataset$patid,
  dstartdate = cprd_dataset$dstartdate,
  prop.score.MFN = drug.pscores$propensity[,1],
  prop.score.GLP1 = drug.pscores$propensity[,2],
  prop.score.DPP4 = drug.pscores$propensity[,3],
  prop.score.SGLT2 = drug.pscores$propensity[,4],
  prop.score.SU = drug.pscores$propensity[,5],
  prop.score.TZD = drug.pscores$propensity[,6],
  prop.score = drug.pscores$`ps.weights` %>% as.data.frame() %>% select(overlap) %>% unlist()
)

saveRDS(ps.only_dataset, "results/PS_model/ps.dataset_lm_all.rds")



##########################################
#  Check balance
#

# drug.pscores <- readRDS("results/PS_model/drug.pscores_model.rds")

pdf("results/figures/covariate_balance.pdf", width = 12, height = 15)
plot(drug.pscores)
dev.off()









# lm.hba1c <- lm(formula = as.formula(paste0("stopdrug_3m_6mFU ~ ", paste0(c(
#   # Biomarkers
#   "precreatinine_blood",
#   "prealt", "pretotalcholesterol", "predbp", "presbp", "prehba1c",
#   "preegfr", "prebilirubin",
#   # Commorbidities
#   "preckdstage", "predrug_frailty_mild", "predrug_frailty_moderate",
#   "predrug_frailty_severe", "predrug_primary_hhf", "predrug_af", "predrug_angina",
#   "predrug_asthma", "predrug_bronchiectasis", "predrug_cld", "predrug_copd",
#   "predrug_cysticfibrosis", "predrug_dementia", "predrug_diabeticnephropathy",
#   "predrug_fh_premature_cvd",
#   "predrug_haem_cancer", "predrug_heartfailure",
#   "predrug_hypertension", "predrug_ihd", "predrug_myocardialinfarction",
#   "predrug_neuropathy", "predrug_otherneuroconditions", "predrug_pad",
#   "predrug_pulmonaryfibrosis", "predrug_pulmonaryhypertension",
#   "predrug_retinopathy", "predrug_revasc", "predrug_rheumatoidarthritis",
#   "predrug_solid_cancer", "predrug_solidorgantransplant", "predrug_stroke",
#   "predrug_tia", "predrug_anxiety_disorders", "predrug_medspecific_gi",
#   "predrug_benignprostatehyperplasia", "predrug_micturition_control",
#   "predrug_volume_depletion", "predrug_urinary_frequency", "predrug_falls",
#   "predrug_lowerlimbfracture", "predrug_incident_mi", "predrug_incident_stroke",
#   "predrug_dka", "predrug_osteoporosis", "predrug_unstableangina",
#   "predrug_amputation",
#   "hosp_admission_prev_year", "hosp_admission_prev_year_count",
#   # Extra info
#   "gender", "prac_region", "ethnicity_5cat", "imd2015_10", "dm_diag_age",
#   "ins_in_1_year", "prebmi", "smoking_cat", "drugline",
#   "alcohol_cat", "dstartdate_age", "dstartdate_dm_dur", "dstartmonth",
#   "CCI_index"
# ), collapse = "+"))),
#                data = ps.dataset,
#                weights = drug.pscores$pw.weights$overlap)


