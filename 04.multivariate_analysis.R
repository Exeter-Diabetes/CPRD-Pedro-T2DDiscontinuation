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
cprd_dataset <- set_up_data(
  raw_data = "20240308_t2d_1stinstance",
  diagnosis = FALSE,
  therapies = c("DPP4", "GLP1", "MFN", "SGLT2", "SU", "TZD"),
  dataset = "full.dataset"
) %>%
  drop_na()


# load propensity scores
ps.only_dataset <- readRDS("results/PS_model/ps.dataset_lm_all.rds")


# join propensity scores into dataset
cprd_dataset <- cprd_dataset %>%
  left_join(
    ps.only_dataset %>%
      select(patid, dstartdate, prop.score.overlap, prop.score.IPW),
    by = c("patid", "dstartdate")
  )


###############################################################################
###############################################################################
########################## Multivariate analysis ##############################
###############################################################################
###############################################################################

# stopdrug_3m_3mFU_MFN_hist
# All therapies pooled
lm_multi_all <- glm(stopdrug_3m_6mFU ~ dstartdate_dm_dur + dstartdate_age + prehba1c + preegfr + prebmi + prealt + pretotalcholesterol + gender + predrug_bloodmed + predrug_cardio_event + predrug_heart_event + predrug_micro_event + preckdstage + predrug_cld, data = cprd_dataset, weights = prop.score.overlap, family=quasibinomial)

# Each therapy
lm_multi_DPP4 <- glm(stopdrug_3m_6mFU ~ dstartdate_dm_dur + dstartdate_age + prehba1c + preegfr + prebmi + prealt + pretotalcholesterol + gender + predrug_bloodmed + predrug_cardio_event + predrug_heart_event + predrug_micro_event + preckdstage + predrug_cld, data = cprd_dataset %>% filter(drugclass == "DPP4"), weights = prop.score.overlap, family=quasibinomial)
lm_multi_GLP1 <- glm(stopdrug_3m_6mFU ~ dstartdate_dm_dur + dstartdate_age + prehba1c + preegfr + prebmi + prealt + pretotalcholesterol + gender + predrug_bloodmed + predrug_cardio_event + predrug_heart_event + predrug_micro_event + preckdstage + predrug_cld, data = cprd_dataset %>% filter(drugclass == "GLP1"), weights = prop.score.overlap, family=quasibinomial)
lm_multi_MFN <- glm(stopdrug_3m_6mFU ~ dstartdate_dm_dur + dstartdate_age + prehba1c + preegfr + prebmi + prealt + pretotalcholesterol + gender + predrug_bloodmed + predrug_cardio_event + predrug_heart_event + predrug_micro_event + preckdstage + predrug_cld, data = cprd_dataset %>% filter(drugclass == "MFN"), weights = prop.score.overlap, family=quasibinomial)
lm_multi_SGLT2 <- glm(stopdrug_3m_6mFU ~ dstartdate_dm_dur + dstartdate_age + prehba1c + preegfr + prebmi + prealt + pretotalcholesterol + gender + predrug_bloodmed + predrug_cardio_event + predrug_heart_event + predrug_micro_event + preckdstage + predrug_cld, data = cprd_dataset %>% filter(drugclass == "SGLT2"), weights = prop.score.overlap, family=quasibinomial)
lm_multi_SU <- glm(stopdrug_3m_6mFU ~ dstartdate_dm_dur + dstartdate_age + prehba1c + preegfr + prebmi + prealt + pretotalcholesterol + gender + predrug_bloodmed + predrug_cardio_event + predrug_heart_event + predrug_micro_event + preckdstage + predrug_cld, data = cprd_dataset %>% filter(drugclass == "SU"), weights = prop.score.overlap, family=quasibinomial)
lm_multi_TZD <- glm(stopdrug_3m_6mFU ~ dstartdate_dm_dur + dstartdate_age + prehba1c + preegfr + prebmi + prealt + pretotalcholesterol + gender + predrug_bloodmed + predrug_cardio_event + predrug_heart_event + predrug_micro_event + preckdstage + predrug_cld, data = cprd_dataset %>% filter(drugclass == "TZD"), weights = prop.score.overlap, family=quasibinomial)



# ROC
plot_roc <- data.frame(
  Model = "Pooled",
  AUC = as.numeric(roc(cprd_dataset %>% select(stopdrug_3m_6mFU) %>% unlist(), predict(lm_multi_all, type = "response"), ci = TRUE)$ci)[2],
  LCI = as.numeric(roc(cprd_dataset %>% select(stopdrug_3m_6mFU) %>% unlist(), predict(lm_multi_all, type = "response"), ci = TRUE)$ci)[1],
  UCI = as.numeric(roc(cprd_dataset %>% select(stopdrug_3m_6mFU) %>% unlist(), predict(lm_multi_all, type = "response"), ci = TRUE)$ci)[3]
) %>%
  rbind(
    data.frame(
      Model = "DPP4",
      AUC = as.numeric(roc(cprd_dataset %>% filter(drugclass == "DPP4") %>% select(stopdrug_3m_6mFU) %>% unlist(), predict(lm_multi_DPP4, type = "response"), ci = TRUE)$ci)[2],
      LCI = as.numeric(roc(cprd_dataset %>% filter(drugclass == "DPP4") %>% select(stopdrug_3m_6mFU) %>% unlist(), predict(lm_multi_DPP4, type = "response"), ci = TRUE)$ci)[1],
      UCI = as.numeric(roc(cprd_dataset %>% filter(drugclass == "DPP4") %>% select(stopdrug_3m_6mFU) %>% unlist(), predict(lm_multi_DPP4, type = "response"), ci = TRUE)$ci)[3]
    ),
    data.frame(
      Model = "GLP1",
      AUC = as.numeric(roc(cprd_dataset %>% filter(drugclass == "GLP1") %>% select(stopdrug_3m_6mFU) %>% unlist(), predict(lm_multi_GLP1, type = "response"), ci = TRUE)$ci)[2],
      LCI = as.numeric(roc(cprd_dataset %>% filter(drugclass == "GLP1") %>% select(stopdrug_3m_6mFU) %>% unlist(), predict(lm_multi_GLP1, type = "response"), ci = TRUE)$ci)[1],
      UCI = as.numeric(roc(cprd_dataset %>% filter(drugclass == "GLP1") %>% select(stopdrug_3m_6mFU) %>% unlist(), predict(lm_multi_GLP1, type = "response"), ci = TRUE)$ci)[3]
    ),
    data.frame(
      Model = "MFN",
      AUC = as.numeric(roc(cprd_dataset %>% filter(drugclass == "MFN") %>% select(stopdrug_3m_6mFU) %>% unlist(), predict(lm_multi_MFN, type = "response"), ci = TRUE)$ci)[2],
      LCI = as.numeric(roc(cprd_dataset %>% filter(drugclass == "MFN") %>% select(stopdrug_3m_6mFU) %>% unlist(), predict(lm_multi_MFN, type = "response"), ci = TRUE)$ci)[1],
      UCI = as.numeric(roc(cprd_dataset %>% filter(drugclass == "MFN") %>% select(stopdrug_3m_6mFU) %>% unlist(), predict(lm_multi_MFN, type = "response"), ci = TRUE)$ci)[3]
    ),
    data.frame(
      Model = "SGLT2",
      AUC = as.numeric(roc(cprd_dataset %>% filter(drugclass == "SGLT2") %>% select(stopdrug_3m_6mFU) %>% unlist(), predict(lm_multi_SGLT2, type = "response"), ci = TRUE)$ci)[2],
      LCI = as.numeric(roc(cprd_dataset %>% filter(drugclass == "SGLT2") %>% select(stopdrug_3m_6mFU) %>% unlist(), predict(lm_multi_SGLT2, type = "response"), ci = TRUE)$ci)[1],
      UCI = as.numeric(roc(cprd_dataset %>% filter(drugclass == "SGLT2") %>% select(stopdrug_3m_6mFU) %>% unlist(), predict(lm_multi_SGLT2, type = "response"), ci = TRUE)$ci)[3]
    ),
    data.frame(
      Model = "SU",
      AUC = as.numeric(roc(cprd_dataset %>% filter(drugclass == "SU") %>% select(stopdrug_3m_6mFU) %>% unlist(), predict(lm_multi_SU, type = "response"), ci = TRUE)$ci)[2],
      LCI = as.numeric(roc(cprd_dataset %>% filter(drugclass == "SU") %>% select(stopdrug_3m_6mFU) %>% unlist(), predict(lm_multi_SU, type = "response"), ci = TRUE)$ci)[1],
      UCI = as.numeric(roc(cprd_dataset %>% filter(drugclass == "SU") %>% select(stopdrug_3m_6mFU) %>% unlist(), predict(lm_multi_SU, type = "response"), ci = TRUE)$ci)[3]
    ),
    data.frame(
      Model = "TZD",
      AUC = as.numeric(roc(cprd_dataset %>% filter(drugclass == "TZD") %>% select(stopdrug_3m_6mFU) %>% unlist(), predict(lm_multi_TZD, type = "response"), ci = TRUE)$ci)[2],
      LCI = as.numeric(roc(cprd_dataset %>% filter(drugclass == "TZD") %>% select(stopdrug_3m_6mFU) %>% unlist(), predict(lm_multi_TZD, type = "response"), ci = TRUE)$ci)[1],
      UCI = as.numeric(roc(cprd_dataset %>% filter(drugclass == "TZD") %>% select(stopdrug_3m_6mFU) %>% unlist(), predict(lm_multi_TZD, type = "response"), ci = TRUE)$ci)[3]
    )
  ) %>%
  mutate(
    Model = factor(Model, levels = rev(c("Pooled", "DPP4", "GLP1", "MFN", "SGLT2", "SU", "TZD"))),
    AUC = as.numeric(AUC),
    LCI = as.numeric(LCI),
    UCI = as.numeric(UCI)
  )  %>%
  ggplot() + 
  geom_errorbar(aes(x = AUC, xmin = LCI, xmax = UCI, y = Model, colour = Model)) +
  geom_pointrange(aes(x = AUC, xmin = LCI, xmax = UCI, y = Model, colour = Model)) +
  scale_colour_manual(values = c("Pooled" = "black", "SGLT2" = "#E69F00", "GLP1" = "#56B4E9", "SU" = "#CC79A7", "DPP4" = "#0072B2", "TZD" = "#D55E00")) +
  theme_bw() +
  theme(
    legend.position = "none"
  )




###############################################################################
###############################################################################
########################### UBER model analysis ###############################
###############################################################################
###############################################################################

# roc(cprd_dataset %>% select(stopdrug_3m_6mFU) %>% unlist(), predict(lm_uber_model, type = "response"), ci = TRUE)

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
                            stopdrug_3m_3mFU_MFN_hist*drugclass +
                            predrug_statins +
                            predrug_bloodmed + 
                            predrug_cardio_event + 
                            predrug_heart_event + 
                            predrug_micro_event + 
                            preckdstage + 
                            predrug_cld"


# Overlap weights
lm_uber_model.overlap_adjust <- glm(as.formula(paste(formula_glm, "+ drugline + numdrugs + smoking_cat + imd2015_10")), 
                     data = cprd_dataset, weights = prop.score.overlap, family=quasibinomial)

lm_uber_model.overlap <- glm(as.formula(formula_glm), 
                                    data = cprd_dataset, weights = prop.score.overlap, family=quasibinomial)

roc(cprd_dataset %>% select(stopdrug_3m_6mFU) %>% unlist(), predict(lm_uber_model.overlap, type = "response"), ci = TRUE)


# IPW weights
lm_uber_model.IPW_adjust <- glm(as.formula(paste(formula_glm, "+ drugline + numdrugs + smoking_cat + imd2015_10")), 
                             data = cprd_dataset, weights = prop.score.IPW, family=quasibinomial)

lm_uber_model.IPW <- glm(as.formula(formula_glm), 
                                data = cprd_dataset, weights = prop.score.IPW, family=quasibinomial)

roc(cprd_dataset %>% select(stopdrug_3m_6mFU) %>% unlist(), predict(lm_uber_model.IPW, type = "response"), ci = TRUE)

# No weights
lm_uber_model.no_weight_adjust <- glm(as.formula(paste(formula_glm, "+ drugline + numdrugs + smoking_cat + imd2015_10")), 
                         data = cprd_dataset, family=quasibinomial)


roc(cprd_dataset %>% select(stopdrug_3m_6mFU) %>% unlist(), predict(lm_uber_model.no_weight_adjust, type = "response"), ci = TRUE)


probabilities.only_dataset <- data.frame(
  patid = cprd_dataset$patid,
  dstartdate = cprd_dataset$dstartdate,
  
  pred.overlap.MFN = predict(lm_uber_model.overlap, newdata = cprd_dataset %>% mutate(drugclass = factor("MFN", levels = levels(cprd_dataset$drugclass))) %>% as.data.frame(), type = "response"),
  pred.overlap.GLP1 = predict(lm_uber_model.overlap, newdata = cprd_dataset %>% mutate(drugclass = factor("GLP1", levels = levels(cprd_dataset$drugclass))) %>% as.data.frame(), type = "response"),
  pred.overlap.DPP4 = predict(lm_uber_model.overlap, newdata = cprd_dataset %>% mutate(drugclass = factor("DPP4", levels = levels(cprd_dataset$drugclass))) %>% as.data.frame(), type = "response"),
  pred.overlap.SGLT2 = predict(lm_uber_model.overlap, newdata = cprd_dataset %>% mutate(drugclass = factor("SGLT2", levels = levels(cprd_dataset$drugclass))) %>% as.data.frame(), type = "response"),
  pred.overlap.SU = predict(lm_uber_model.overlap, newdata = cprd_dataset %>% mutate(drugclass = factor("SU", levels = levels(cprd_dataset$drugclass))) %>% as.data.frame(), type = "response"),
  pred.overlap.TZD = predict(lm_uber_model.overlap, newdata = cprd_dataset %>% mutate(drugclass = factor("TZD", levels = levels(cprd_dataset$drugclass))) %>% as.data.frame(), type = "response"),
  
  pred.overlap_adjust.MFN = predict(lm_uber_model.overlap_adjust, newdata = cprd_dataset %>% mutate(drugclass = factor("MFN", levels = levels(cprd_dataset$drugclass))) %>% as.data.frame(), type = "response"),
  pred.overlap_adjust.GLP1 = predict(lm_uber_model.overlap_adjust, newdata = cprd_dataset %>% mutate(drugclass = factor("GLP1", levels = levels(cprd_dataset$drugclass))) %>% as.data.frame(), type = "response"),
  pred.overlap_adjust.DPP4 = predict(lm_uber_model.overlap_adjust, newdata = cprd_dataset %>% mutate(drugclass = factor("DPP4", levels = levels(cprd_dataset$drugclass))) %>% as.data.frame(), type = "response"),
  pred.overlap_adjust.SGLT2 = predict(lm_uber_model.overlap_adjust, newdata = cprd_dataset %>% mutate(drugclass = factor("SGLT2", levels = levels(cprd_dataset$drugclass))) %>% as.data.frame(), type = "response"),
  pred.overlap_adjust.SU = predict(lm_uber_model.overlap_adjust, newdata = cprd_dataset %>% mutate(drugclass = factor("SU", levels = levels(cprd_dataset$drugclass))) %>% as.data.frame(), type = "response"),
  pred.overlap_adjust.TZD = predict(lm_uber_model.overlap_adjust, newdata = cprd_dataset %>% mutate(drugclass = factor("TZD", levels = levels(cprd_dataset$drugclass))) %>% as.data.frame(), type = "response"),
  
  pred.IPW.MFN = predict(lm_uber_model.IPW, newdata = cprd_dataset %>% mutate(drugclass = factor("MFN", levels = levels(cprd_dataset$drugclass))) %>% as.data.frame(), type = "response"),
  pred.IPW.GLP1 = predict(lm_uber_model.IPW, newdata = cprd_dataset %>% mutate(drugclass = factor("GLP1", levels = levels(cprd_dataset$drugclass))) %>% as.data.frame(), type = "response"),
  pred.IPW.DPP4 = predict(lm_uber_model.IPW, newdata = cprd_dataset %>% mutate(drugclass = factor("DPP4", levels = levels(cprd_dataset$drugclass))) %>% as.data.frame(), type = "response"),
  pred.IPW.SGLT2 = predict(lm_uber_model.IPW, newdata = cprd_dataset %>% mutate(drugclass = factor("SGLT2", levels = levels(cprd_dataset$drugclass))) %>% as.data.frame(), type = "response"),
  pred.IPW.SU = predict(lm_uber_model.IPW, newdata = cprd_dataset %>% mutate(drugclass = factor("SU", levels = levels(cprd_dataset$drugclass))) %>% as.data.frame(), type = "response"),
  pred.IPW.TZD = predict(lm_uber_model.IPW, newdata = cprd_dataset %>% mutate(drugclass = factor("TZD", levels = levels(cprd_dataset$drugclass))) %>% as.data.frame(), type = "response"),
  
  pred.IPW_adjust.MFN = predict(lm_uber_model.IPW_adjust, newdata = cprd_dataset %>% mutate(drugclass = factor("MFN", levels = levels(cprd_dataset$drugclass))) %>% as.data.frame(), type = "response"),
  pred.IPW_adjust.GLP1 = predict(lm_uber_model.IPW_adjust, newdata = cprd_dataset %>% mutate(drugclass = factor("GLP1", levels = levels(cprd_dataset$drugclass))) %>% as.data.frame(), type = "response"),
  pred.IPW_adjust.DPP4 = predict(lm_uber_model.IPW_adjust, newdata = cprd_dataset %>% mutate(drugclass = factor("DPP4", levels = levels(cprd_dataset$drugclass))) %>% as.data.frame(), type = "response"),
  pred.IPW_adjust.SGLT2 = predict(lm_uber_model.IPW_adjust, newdata = cprd_dataset %>% mutate(drugclass = factor("SGLT2", levels = levels(cprd_dataset$drugclass))) %>% as.data.frame(), type = "response"),
  pred.IPW_adjust.SU = predict(lm_uber_model.IPW_adjust, newdata = cprd_dataset %>% mutate(drugclass = factor("SU", levels = levels(cprd_dataset$drugclass))) %>% as.data.frame(), type = "response"),
  pred.IPW_adjust.TZD = predict(lm_uber_model.IPW_adjust, newdata = cprd_dataset %>% mutate(drugclass = factor("TZD", levels = levels(cprd_dataset$drugclass))) %>% as.data.frame(), type = "response"),
  
  pred.no_weight.MFN = predict(lm_uber_model.no_weight_adjust, newdata = cprd_dataset %>% mutate(drugclass = factor("MFN", levels = levels(cprd_dataset$drugclass))) %>% as.data.frame(), type = "response"),
  pred.no_weight.GLP1 = predict(lm_uber_model.no_weight_adjust, newdata = cprd_dataset %>% mutate(drugclass = factor("GLP1", levels = levels(cprd_dataset$drugclass))) %>% as.data.frame(), type = "response"),
  pred.no_weight.DPP4 = predict(lm_uber_model.no_weight_adjust, newdata = cprd_dataset %>% mutate(drugclass = factor("DPP4", levels = levels(cprd_dataset$drugclass))) %>% as.data.frame(), type = "response"),
  pred.no_weight.SGLT2 = predict(lm_uber_model.no_weight_adjust, newdata = cprd_dataset %>% mutate(drugclass = factor("SGLT2", levels = levels(cprd_dataset$drugclass))) %>% as.data.frame(), type = "response"),
  pred.no_weight.SU = predict(lm_uber_model.no_weight_adjust, newdata = cprd_dataset %>% mutate(drugclass = factor("SU", levels = levels(cprd_dataset$drugclass))) %>% as.data.frame(), type = "response"),
  pred.no_weight.TZD = predict(lm_uber_model.no_weight_adjust, newdata = cprd_dataset %>% mutate(drugclass = factor("TZD", levels = levels(cprd_dataset$drugclass))) %>% as.data.frame(), type = "response")

)

saveRDS(probabilities.only_dataset, "results/Models/Predictions/model_predictions.rds")




cprd_dataset <- cprd_dataset %>%
  mutate(
    pred.overlap = predict(lm_uber_model.overlap, type = "response"),
    pred.overlap_adjust = predict(lm_uber_model.overlap_adjust, type = "response"),
    pred.IPW = predict(lm_uber_model.IPW, type = "response"),
    pred.IPW_adjust = predict(lm_uber_model.IPW_adjust, type = "response"),
    pred.no_weight = predict(lm_uber_model.no_weight_adjust, type = "response")
    )


pdf("results/figures/uber_model_calibration.pdf", width = 4, height = 4)

calibration_plot(
  data = cprd_dataset %>% mutate(stopdrug_3m_6mFU = as.numeric(stopdrug_3m_6mFU)-1), 
  obs = "stopdrug_3m_6mFU", 
  pred = "pred.overlap", 
  title = "Overlap weighting", 
  y_lim = c(0, 0.3), 
  x_lim=c(0, 0.3), 
  nTiles = 10)
calibration_plot(
  data = cprd_dataset %>% mutate(stopdrug_3m_6mFU = as.numeric(stopdrug_3m_6mFU)-1), 
  obs = "stopdrug_3m_6mFU", 
  pred = "pred.overlap_adjust", 
  title = "Overlap weighting + estimate adjustment", 
  y_lim = c(0, 0.3), 
  x_lim=c(0, 0.3), 
  nTiles = 10)
calibration_plot(
  data = cprd_dataset %>% mutate(stopdrug_3m_6mFU = as.numeric(stopdrug_3m_6mFU)-1), 
  obs = "stopdrug_3m_6mFU", 
  pred = "pred.IPW", 
  title = "IPW weighting", 
  y_lim = c(0, 0.3), 
  x_lim=c(0, 0.3), 
  nTiles = 10)
calibration_plot(
  data = cprd_dataset %>% mutate(stopdrug_3m_6mFU = as.numeric(stopdrug_3m_6mFU)-1), 
  obs = "stopdrug_3m_6mFU", 
  pred = "pred.IPW_adjust", 
  title = "IPW weighting + estimate adjustment", 
  y_lim = c(0, 0.3), 
  x_lim=c(0, 0.3), 
  nTiles = 10)
calibration_plot(
  data = cprd_dataset %>% mutate(stopdrug_3m_6mFU = as.numeric(stopdrug_3m_6mFU)-1), 
  obs = "stopdrug_3m_6mFU", 
  pred = "pred.no_weight", 
  title = "Model variables adjustment", 
  y_lim = c(0, 0.3), 
  x_lim=c(0, 0.3), 
  nTiles = 10)

dev.off()


