####################
## Description: 
##  - In this file we:
##    - Fit a linear model all the variables into one model.
#################### 


# load working directory
setwd("Samples/T2D_Discontinuation")

# load functions
source("code/00.set_up_data.R")

# load libraries
library(tidyverse)
library(pROC)
library(predtools)

###############################################################################
###############################################################################
############################### Read Data In ##################################
###############################################################################
###############################################################################

# load dataset
cprd_dataset.dev <- set_up_data(
  raw_data = "20240308_t2d_1stinstance",
  diagnosis = FALSE,
  therapies = c("DPP4", "GLP1", "MFN", "SGLT2", "SU", "TZD"),
  dataset = "3m.disc.dataset.dev"
) %>%
  drop_na()


# modify some variables for fitting
cprd_dataset.dev <- cprd_dataset.dev %>%
  mutate(
    ethnicity_5cat = factor(ethnicity_5cat, levels = c(0, 1, 2, 3, 4), labels = c("White", "Other", "Other" , "Other" , "Other")),
    smoking_cat = factor(smoking_cat, levels = c("Active smoker", "Ex-smoker", "Non-smoker"), labels = c("Active smoker", "Other", "Other")),
    alcohol_cat = factor(alcohol_cat, levels = c("Excess", "Harmful", "None", "Within limits"), labels = c("Excess/Harmful", "Excess/Harmful", "Within limits/None", "Within limits/None")),
    preckdstage = factor(preckdstage, levels = c("stage_0", "stage_1", "stage_2", "stage_3a", "stage_3b", "stage_4", "stage_5"), labels = c("stage 0/1/2", "stage 0/1/2", "stage 0/1/2", "stage 3/4/5", "stage 3/4/5", "stage 3/4/5", "stage 3/4/5"))
  )


# load dataset
cprd_dataset.val <- set_up_data(
  raw_data = "20240308_t2d_1stinstance",
  diagnosis = FALSE,
  therapies = c("DPP4", "GLP1", "MFN", "SGLT2", "SU", "TZD"),
  dataset = "3m.disc.dataset.val"
) %>%
  drop_na()


# modify some variables for fitting
cprd_dataset.val <- cprd_dataset.val %>%
  mutate(
    ethnicity_5cat = factor(ethnicity_5cat, levels = c(0, 1, 2, 3, 4), labels = c("White", "Other", "Other" , "Other" , "Other")),
    smoking_cat = factor(smoking_cat, levels = c("Active smoker", "Ex-smoker", "Non-smoker"), labels = c("Active smoker", "Other", "Other")),
    alcohol_cat = factor(alcohol_cat, levels = c("Excess", "Harmful", "None", "Within limits"), labels = c("Excess/Harmful", "Excess/Harmful", "Within limits/None", "Within limits/None")),
    preckdstage = factor(preckdstage, levels = c("stage_0", "stage_1", "stage_2", "stage_3a", "stage_3b", "stage_4", "stage_5"), labels = c("stage 0/1/2", "stage 0/1/2", "stage 0/1/2", "stage 3/4/5", "stage 3/4/5", "stage 3/4/5", "stage 3/4/5"))
  )


###############################################################################
###############################################################################
########################## Multivariate analysis ##############################
###############################################################################
###############################################################################



###############################################################################
###############################################################################
########################### UBER model analysis ###############################
###############################################################################
###############################################################################

##:--------- All drugs in the model

# roc(cprd_dataset.dev %>% select(stopdrug_3m_6mFU) %>% unlist(), predict(lm_uber_model, type = "response"), ci = TRUE)

formula_glm <- "stopdrug_3m_6mFU ~ 
                            drugclass +
                            dstartdate_dm_dur*drugclass + 
                            dstartdate_age*drugclass + 
                            prehba1c*drugclass + 
                            preegfr*drugclass + 
                            prebmi*drugclass + 
                            prealt*drugclass + 
                            pretotalcholesterol*drugclass + 
                            gender*drugclass + 
                            stopdrug_3m_3mFU_MFN_hist +
                            ethnicity_5cat +
                            smoking_cat +
                            alcohol_cat +
                            predrug_statins +
                            predrug_bloodmed + 
                            predrug_cardio_event + 
                            predrug_heart_event + 
                            predrug_micro_event + 
                            preckdstage + 
                            predrug_cld +
                            predrug_frailty_proxy + 
                            drugline + 
                            numdrugs + 
                            imd2015_10"



lm_uber_model.no_weight_adjust <- glm(as.formula(formula_glm), data = cprd_dataset.dev, family = binomial())


roc(cprd_dataset.dev %>% select(stopdrug_3m_6mFU) %>% unlist(), predict(lm_uber_model.no_weight_adjust, type = "response"), ci = TRUE)



probabilities.only_dataset <- data.frame(
  
  ## Development dataset
  patid = cprd_dataset.dev$patid,
  dstartdate = cprd_dataset.dev$dstartdate,

  pred.no_weight = predict(lm_uber_model.no_weight_adjust, newdata = cprd_dataset.dev, type = "response"),
  pred.no_weight.MFN = predict(lm_uber_model.no_weight_adjust, newdata = cprd_dataset.dev %>% mutate(drugclass = factor("MFN", levels = levels(cprd_dataset.dev$drugclass))) %>% as.data.frame(), type = "response"),
  pred.no_weight.GLP1 = predict(lm_uber_model.no_weight_adjust, newdata = cprd_dataset.dev %>% mutate(drugclass = factor("GLP1", levels = levels(cprd_dataset.dev$drugclass))) %>% as.data.frame(), type = "response"),
  pred.no_weight.DPP4 = predict(lm_uber_model.no_weight_adjust, newdata = cprd_dataset.dev %>% mutate(drugclass = factor("DPP4", levels = levels(cprd_dataset.dev$drugclass))) %>% as.data.frame(), type = "response"),
  pred.no_weight.SGLT2 = predict(lm_uber_model.no_weight_adjust, newdata = cprd_dataset.dev %>% mutate(drugclass = factor("SGLT2", levels = levels(cprd_dataset.dev$drugclass))) %>% as.data.frame(), type = "response"),
  pred.no_weight.SU = predict(lm_uber_model.no_weight_adjust, newdata = cprd_dataset.dev %>% mutate(drugclass = factor("SU", levels = levels(cprd_dataset.dev$drugclass))) %>% as.data.frame(), type = "response"),
  pred.no_weight.TZD = predict(lm_uber_model.no_weight_adjust, newdata = cprd_dataset.dev %>% mutate(drugclass = factor("TZD", levels = levels(cprd_dataset.dev$drugclass))) %>% as.data.frame(), type = "response")
  
) %>%
  
  ## Validation dataset
  rbind(
    data.frame(
      patid = cprd_dataset.val$patid,
      dstartdate = cprd_dataset.val$dstartdate,
      
      pred.no_weight = predict(lm_uber_model.no_weight_adjust, newdata = cprd_dataset.val, type = "response"),
      pred.no_weight.MFN = predict(lm_uber_model.no_weight_adjust, newdata = cprd_dataset.val %>% mutate(drugclass = factor("MFN", levels = levels(cprd_dataset.val$drugclass))) %>% as.data.frame(), type = "response"),
      pred.no_weight.GLP1 = predict(lm_uber_model.no_weight_adjust, newdata = cprd_dataset.val %>% mutate(drugclass = factor("GLP1", levels = levels(cprd_dataset.val$drugclass))) %>% as.data.frame(), type = "response"),
      pred.no_weight.DPP4 = predict(lm_uber_model.no_weight_adjust, newdata = cprd_dataset.val %>% mutate(drugclass = factor("DPP4", levels = levels(cprd_dataset.val$drugclass))) %>% as.data.frame(), type = "response"),
      pred.no_weight.SGLT2 = predict(lm_uber_model.no_weight_adjust, newdata = cprd_dataset.val %>% mutate(drugclass = factor("SGLT2", levels = levels(cprd_dataset.val$drugclass))) %>% as.data.frame(), type = "response"),
      pred.no_weight.SU = predict(lm_uber_model.no_weight_adjust, newdata = cprd_dataset.val %>% mutate(drugclass = factor("SU", levels = levels(cprd_dataset.val$drugclass))) %>% as.data.frame(), type = "response"),
      pred.no_weight.TZD = predict(lm_uber_model.no_weight_adjust, newdata = cprd_dataset.val %>% mutate(drugclass = factor("TZD", levels = levels(cprd_dataset.val$drugclass))) %>% as.data.frame(), type = "response")
      
    )
  )

saveRDS(probabilities.only_dataset, "results/Models/Predictions/model_predictions_3m_all_drugs.rds")


cprd_dataset.dev <- cprd_dataset.dev %>%
  mutate(
    pred.no_weight = predict(lm_uber_model.no_weight_adjust, type = "response")
  )

cprd_dataset.val <- cprd_dataset.val %>%
  mutate(
    pred.no_weight = predict(lm_uber_model.no_weight_adjust, newdata = cprd_dataset.val, type = "response")
  )

pdf("results/figures/uber_model_calibration_3m_all_drugs.pdf", width = 4, height = 4)

calibration_plot(
  data = cprd_dataset.dev %>% mutate(stopdrug_3m_6mFU = as.numeric(stopdrug_3m_6mFU)-1), 
  obs = "stopdrug_3m_6mFU", 
  pred = "pred.no_weight", 
  title = "Dev: Model variables adjustment", 
  y_lim = c(0, 0.3), 
  x_lim=c(0, 0.3), 
  nTiles = 10)

calibration_plot(
  data = cprd_dataset.val %>% mutate(stopdrug_3m_6mFU = as.numeric(stopdrug_3m_6mFU)-1), 
  obs = "stopdrug_3m_6mFU", 
  pred = "pred.no_weight", 
  title = "Val: Model variables adjustment", 
  y_lim = c(0, 0.3), 
  x_lim=c(0, 0.3), 
  nTiles = 10)

dev.off()




##:----------- Only 2nd line therapies

# load dataset
cprd_dataset.dev <- set_up_data(
  raw_data = "20240308_t2d_1stinstance",
  diagnosis = FALSE,
  therapies = c("DPP4", "GLP1", "SGLT2", "SU", "TZD"),
  dataset = "3m.disc.dataset.dev"
) %>%
  drop_na()


# modify some variables for fitting
cprd_dataset.dev <- cprd_dataset.dev %>%
  mutate(
    ethnicity_5cat = factor(ethnicity_5cat, levels = c(0, 1, 2, 3, 4), labels = c("White", "Other", "Other" , "Other" , "Other")),
    smoking_cat = factor(smoking_cat, levels = c("Active smoker", "Ex-smoker", "Non-smoker"), labels = c("Active smoker", "Other", "Other")),
    alcohol_cat = factor(alcohol_cat, levels = c("Excess", "Harmful", "None", "Within limits"), labels = c("Excess/Harmful", "Excess/Harmful", "Within limits/None", "Within limits/None")),
    preckdstage = factor(preckdstage, levels = c("stage_0", "stage_1", "stage_2", "stage_3a", "stage_3b", "stage_4", "stage_5"), labels = c("stage 0/1/2", "stage 0/1/2", "stage 0/1/2", "stage 3/4/5", "stage 3/4/5", "stage 3/4/5", "stage 3/4/5"))
  )


# load dataset
cprd_dataset.val <- set_up_data(
  raw_data = "20240308_t2d_1stinstance",
  diagnosis = FALSE,
  therapies = c("DPP4", "GLP1", "SGLT2", "SU", "TZD"),
  dataset = "3m.disc.dataset.val"
) %>%
  drop_na()


# modify some variables for fitting
cprd_dataset.val <- cprd_dataset.val %>%
  mutate(
    ethnicity_5cat = factor(ethnicity_5cat, levels = c(0, 1, 2, 3, 4), labels = c("White", "Other", "Other" , "Other" , "Other")),
    smoking_cat = factor(smoking_cat, levels = c("Active smoker", "Ex-smoker", "Non-smoker"), labels = c("Active smoker", "Other", "Other")),
    alcohol_cat = factor(alcohol_cat, levels = c("Excess", "Harmful", "None", "Within limits"), labels = c("Excess/Harmful", "Excess/Harmful", "Within limits/None", "Within limits/None")),
    preckdstage = factor(preckdstage, levels = c("stage_0", "stage_1", "stage_2", "stage_3a", "stage_3b", "stage_4", "stage_5"), labels = c("stage 0/1/2", "stage 0/1/2", "stage 0/1/2", "stage 3/4/5", "stage 3/4/5", "stage 3/4/5", "stage 3/4/5"))
  )


lm_uber_model.no_weight_adjust <- glm(as.formula(formula_glm), data = cprd_dataset.dev, family = binomial())


roc(cprd_dataset.dev %>% select(stopdrug_3m_6mFU) %>% unlist(), predict(lm_uber_model.no_weight_adjust, type = "response"), ci = TRUE)



probabilities.only_dataset <- data.frame(
  
  ## Development dataset
  patid = cprd_dataset.dev$patid,
  dstartdate = cprd_dataset.dev$dstartdate,
  
  pred.no_weight = predict(lm_uber_model.no_weight_adjust, newdata = cprd_dataset.dev, type = "response"),
  pred.no_weight.GLP1 = predict(lm_uber_model.no_weight_adjust, newdata = cprd_dataset.dev %>% mutate(drugclass = factor("GLP1", levels = levels(cprd_dataset.dev$drugclass))) %>% as.data.frame(), type = "response"),
  pred.no_weight.DPP4 = predict(lm_uber_model.no_weight_adjust, newdata = cprd_dataset.dev %>% mutate(drugclass = factor("DPP4", levels = levels(cprd_dataset.dev$drugclass))) %>% as.data.frame(), type = "response"),
  pred.no_weight.SGLT2 = predict(lm_uber_model.no_weight_adjust, newdata = cprd_dataset.dev %>% mutate(drugclass = factor("SGLT2", levels = levels(cprd_dataset.dev$drugclass))) %>% as.data.frame(), type = "response"),
  pred.no_weight.SU = predict(lm_uber_model.no_weight_adjust, newdata = cprd_dataset.dev %>% mutate(drugclass = factor("SU", levels = levels(cprd_dataset.dev$drugclass))) %>% as.data.frame(), type = "response"),
  pred.no_weight.TZD = predict(lm_uber_model.no_weight_adjust, newdata = cprd_dataset.dev %>% mutate(drugclass = factor("TZD", levels = levels(cprd_dataset.dev$drugclass))) %>% as.data.frame(), type = "response")
  
) %>%
  
  ## Validation dataset
  rbind(
    data.frame(
      patid = cprd_dataset.val$patid,
      dstartdate = cprd_dataset.val$dstartdate,
      
      pred.no_weight = predict(lm_uber_model.no_weight_adjust, newdata = cprd_dataset.val, type = "response"),
      pred.no_weight.GLP1 = predict(lm_uber_model.no_weight_adjust, newdata = cprd_dataset.val %>% mutate(drugclass = factor("GLP1", levels = levels(cprd_dataset.val$drugclass))) %>% as.data.frame(), type = "response"),
      pred.no_weight.DPP4 = predict(lm_uber_model.no_weight_adjust, newdata = cprd_dataset.val %>% mutate(drugclass = factor("DPP4", levels = levels(cprd_dataset.val$drugclass))) %>% as.data.frame(), type = "response"),
      pred.no_weight.SGLT2 = predict(lm_uber_model.no_weight_adjust, newdata = cprd_dataset.val %>% mutate(drugclass = factor("SGLT2", levels = levels(cprd_dataset.val$drugclass))) %>% as.data.frame(), type = "response"),
      pred.no_weight.SU = predict(lm_uber_model.no_weight_adjust, newdata = cprd_dataset.val %>% mutate(drugclass = factor("SU", levels = levels(cprd_dataset.val$drugclass))) %>% as.data.frame(), type = "response"),
      pred.no_weight.TZD = predict(lm_uber_model.no_weight_adjust, newdata = cprd_dataset.val %>% mutate(drugclass = factor("TZD", levels = levels(cprd_dataset.val$drugclass))) %>% as.data.frame(), type = "response")
      
    )
  )

saveRDS(probabilities.only_dataset, "results/Models/Predictions/model_predictions_3m_2nd_line_drugs.rds")


cprd_dataset.dev <- cprd_dataset.dev %>%
  mutate(
    pred.no_weight = predict(lm_uber_model.no_weight_adjust, type = "response")
  )

cprd_dataset.val <- cprd_dataset.val %>%
  mutate(
    pred.no_weight = predict(lm_uber_model.no_weight_adjust, newdata = cprd_dataset.val, type = "response")
  )

pdf("results/figures/uber_model_calibration_3m_2nd_line_drugs.pdf", width = 4, height = 4)

calibration_plot(
  data = cprd_dataset.dev %>% mutate(stopdrug_3m_6mFU = as.numeric(stopdrug_3m_6mFU)-1), 
  obs = "stopdrug_3m_6mFU", 
  pred = "pred.no_weight", 
  title = "Dev: Model variables adjustment", 
  y_lim = c(0, 0.3), 
  x_lim=c(0, 0.3), 
  nTiles = 10)

calibration_plot(
  data = cprd_dataset.val %>% mutate(stopdrug_3m_6mFU = as.numeric(stopdrug_3m_6mFU)-1), 
  obs = "stopdrug_3m_6mFU", 
  pred = "pred.no_weight", 
  title = "Val: Model variables adjustment", 
  y_lim = c(0, 0.3), 
  x_lim=c(0, 0.3), 
  nTiles = 10)

dev.off()




