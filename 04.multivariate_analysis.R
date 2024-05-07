####################
## Description: 
##  - In this file we:
##    - Fit a linear model to each of the variables for pooled therapies and by therapy.
#################### 


# load working directory
setwd("Samples/T2D_Discontinuation")

# load functions
source("code/00.set_up_data.R")

# load libraries
library(tidyverse)
library(pROC)

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
      select(patid, dstartdate, prop.score),
    by = c("patid", "dstartdate")
  )


###############################################################################
###############################################################################
########################## Multivariate analysis ##############################
###############################################################################
###############################################################################

# stopdrug_3m_3mFU_MFN_hist
# All therapies pooled
lm_multi_all <- glm(stopdrug_3m_6mFU ~ dstartdate_dm_dur + dstartdate_age + prehba1c + preegfr + prebmi + prealt + pretotalcholesterol + gender + predrug_bloodmed + predrug_cardio_event + predrug_heart_event + predrug_micro_event + preckdstage + predrug_cld, data = cprd_dataset, weights = prop.score, family=quasibinomial)

# Each therapy
lm_multi_DPP4 <- glm(stopdrug_3m_6mFU ~ dstartdate_dm_dur + dstartdate_age + prehba1c + preegfr + prebmi + prealt + pretotalcholesterol + gender + predrug_bloodmed + predrug_cardio_event + predrug_heart_event + predrug_micro_event + preckdstage + predrug_cld, data = cprd_dataset %>% filter(drugclass == "DPP4"), weights = prop.score, family=quasibinomial)
lm_multi_GLP1 <- glm(stopdrug_3m_6mFU ~ dstartdate_dm_dur + dstartdate_age + prehba1c + preegfr + prebmi + prealt + pretotalcholesterol + gender + predrug_bloodmed + predrug_cardio_event + predrug_heart_event + predrug_micro_event + preckdstage + predrug_cld, data = cprd_dataset %>% filter(drugclass == "GLP1"), weights = prop.score, family=quasibinomial)
lm_multi_MFN <- glm(stopdrug_3m_6mFU ~ dstartdate_dm_dur + dstartdate_age + prehba1c + preegfr + prebmi + prealt + pretotalcholesterol + gender + predrug_bloodmed + predrug_cardio_event + predrug_heart_event + predrug_micro_event + preckdstage + predrug_cld, data = cprd_dataset %>% filter(drugclass == "MFN"), weights = prop.score, family=quasibinomial)
lm_multi_SGLT2 <- glm(stopdrug_3m_6mFU ~ dstartdate_dm_dur + dstartdate_age + prehba1c + preegfr + prebmi + prealt + pretotalcholesterol + gender + predrug_bloodmed + predrug_cardio_event + predrug_heart_event + predrug_micro_event + preckdstage + predrug_cld, data = cprd_dataset %>% filter(drugclass == "SGLT2"), weights = prop.score, family=quasibinomial)
lm_multi_SU <- glm(stopdrug_3m_6mFU ~ dstartdate_dm_dur + dstartdate_age + prehba1c + preegfr + prebmi + prealt + pretotalcholesterol + gender + predrug_bloodmed + predrug_cardio_event + predrug_heart_event + predrug_micro_event + preckdstage + predrug_cld, data = cprd_dataset %>% filter(drugclass == "SU"), weights = prop.score, family=quasibinomial)
lm_multi_TZD <- glm(stopdrug_3m_6mFU ~ dstartdate_dm_dur + dstartdate_age + prehba1c + preegfr + prebmi + prealt + pretotalcholesterol + gender + predrug_bloodmed + predrug_cardio_event + predrug_heart_event + predrug_micro_event + preckdstage + predrug_cld, data = cprd_dataset %>% filter(drugclass == "TZD"), weights = prop.score, family=quasibinomial)



# ROC
data.frame(
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












