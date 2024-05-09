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
library(patchwork)

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
########################### Univariate analysis ###############################
###############################################################################
###############################################################################

cprd_dataset <- cprd_dataset %>%
  mutate(dstartdate_dm_dur = dstartdate_dm_dur/10)


# dstartdate_dm_dur
lm_uni_all.dstartdate_dm_dur <- glm(stopdrug_3m_6mFU ~ dstartdate_dm_dur, data = cprd_dataset, weights = prop.score, family=quasibinomial)

lm_uni_DPP4.dstartdate_dm_dur <- glm(stopdrug_3m_6mFU ~ dstartdate_dm_dur, data = cprd_dataset %>% filter(drugclass == "DPP4"), weights = prop.score, family=quasibinomial)
lm_uni_GLP1.dstartdate_dm_dur <- glm(stopdrug_3m_6mFU ~ dstartdate_dm_dur, data = cprd_dataset %>% filter(drugclass == "GLP1"), weights = prop.score, family=quasibinomial)
lm_uni_MFN.dstartdate_dm_dur <- glm(stopdrug_3m_6mFU ~ dstartdate_dm_dur, data = cprd_dataset %>% filter(drugclass == "MFN"), weights = prop.score, family=quasibinomial)
lm_uni_SGLT2.dstartdate_dm_dur <- glm(stopdrug_3m_6mFU ~ dstartdate_dm_dur, data = cprd_dataset %>% filter(drugclass == "SGLT2"), weights = prop.score, family=quasibinomial)
lm_uni_SU.dstartdate_dm_dur <- glm(stopdrug_3m_6mFU ~ dstartdate_dm_dur, data = cprd_dataset %>% filter(drugclass == "SU"), weights = prop.score, family=quasibinomial)
lm_uni_TZD.dstartdate_dm_dur <- glm(stopdrug_3m_6mFU ~ dstartdate_dm_dur, data = cprd_dataset %>% filter(drugclass == "TZD"), weights = prop.score, family=quasibinomial)


cprd_dataset <- cprd_dataset %>%
  mutate(dstartdate_age = dstartdate_age/10)

# dstartdate_age
lm_uni_all.dstartdate_age <- glm(stopdrug_3m_6mFU ~ dstartdate_age, data = cprd_dataset, weights = prop.score, family=quasibinomial)

lm_uni_DPP4.dstartdate_age <- glm(stopdrug_3m_6mFU ~ dstartdate_age, data = cprd_dataset %>% filter(drugclass == "DPP4"), weights = prop.score, family=quasibinomial)
lm_uni_GLP1.dstartdate_age <- glm(stopdrug_3m_6mFU ~ dstartdate_age, data = cprd_dataset %>% filter(drugclass == "GLP1"), weights = prop.score, family=quasibinomial)
lm_uni_MFN.dstartdate_age <- glm(stopdrug_3m_6mFU ~ dstartdate_age, data = cprd_dataset %>% filter(drugclass == "MFN"), weights = prop.score, family=quasibinomial)
lm_uni_SGLT2.dstartdate_age <- glm(stopdrug_3m_6mFU ~ dstartdate_age, data = cprd_dataset %>% filter(drugclass == "SGLT2"), weights = prop.score, family=quasibinomial)
lm_uni_SU.dstartdate_age <- glm(stopdrug_3m_6mFU ~ dstartdate_age, data = cprd_dataset %>% filter(drugclass == "SU"), weights = prop.score, family=quasibinomial)
lm_uni_TZD.dstartdate_age <- glm(stopdrug_3m_6mFU ~ dstartdate_age, data = cprd_dataset %>% filter(drugclass == "TZD"), weights = prop.score, family=quasibinomial)


cprd_dataset <- cprd_dataset %>%
  mutate(prehba1c = prehba1c/10)

# prehba1c
lm_uni_all.prehba1c <- glm(stopdrug_3m_6mFU ~ prehba1c, data = cprd_dataset, weights = prop.score, family=quasibinomial)

lm_uni_DPP4.prehba1c <- glm(stopdrug_3m_6mFU ~ prehba1c, data = cprd_dataset %>% filter(drugclass == "DPP4"), weights = prop.score, family=quasibinomial)
lm_uni_GLP1.prehba1c <- glm(stopdrug_3m_6mFU ~ prehba1c, data = cprd_dataset %>% filter(drugclass == "GLP1"), weights = prop.score, family=quasibinomial)
lm_uni_MFN.prehba1c <- glm(stopdrug_3m_6mFU ~ prehba1c, data = cprd_dataset %>% filter(drugclass == "MFN"), weights = prop.score, family=quasibinomial)
lm_uni_SGLT2.prehba1c <- glm(stopdrug_3m_6mFU ~ prehba1c, data = cprd_dataset %>% filter(drugclass == "SGLT2"), weights = prop.score, family=quasibinomial)
lm_uni_SU.prehba1c <- glm(stopdrug_3m_6mFU ~ prehba1c, data = cprd_dataset %>% filter(drugclass == "SU"), weights = prop.score, family=quasibinomial)
lm_uni_TZD.prehba1c <- glm(stopdrug_3m_6mFU ~ prehba1c, data = cprd_dataset %>% filter(drugclass == "TZD"), weights = prop.score, family=quasibinomial)


cprd <- cprd_dataset %>%
  mutate_at(c("preegfr"), ~(scale(.) %>% as.vector))


# preegfr
lm_uni_all.preegfr <- glm(stopdrug_3m_6mFU ~ preegfr, data = cprd_dataset, weights = prop.score, family=quasibinomial)

lm_uni_DPP4.preegfr <- glm(stopdrug_3m_6mFU ~ preegfr, data = cprd_dataset %>% filter(drugclass == "DPP4"), weights = prop.score, family=quasibinomial)
lm_uni_GLP1.preegfr <- glm(stopdrug_3m_6mFU ~ preegfr, data = cprd_dataset %>% filter(drugclass == "GLP1"), weights = prop.score, family=quasibinomial)
lm_uni_MFN.preegfr <- glm(stopdrug_3m_6mFU ~ preegfr, data = cprd_dataset %>% filter(drugclass == "MFN"), weights = prop.score, family=quasibinomial)
lm_uni_SGLT2.preegfr <- glm(stopdrug_3m_6mFU ~ preegfr, data = cprd_dataset %>% filter(drugclass == "SGLT2"), weights = prop.score, family=quasibinomial)
lm_uni_SU.preegfr <- glm(stopdrug_3m_6mFU ~ preegfr, data = cprd_dataset %>% filter(drugclass == "SU"), weights = prop.score, family=quasibinomial)
lm_uni_TZD.preegfr <- glm(stopdrug_3m_6mFU ~ preegfr, data = cprd_dataset %>% filter(drugclass == "TZD"), weights = prop.score, family=quasibinomial)


cprd_dataset <- cprd_dataset %>%
  mutate(prebmi = prebmi/5)

# prebmi
lm_uni_all.prebmi <- glm(stopdrug_3m_6mFU ~ prebmi, data = cprd_dataset, weights = prop.score, family=quasibinomial)

lm_uni_DPP4.prebmi <- glm(stopdrug_3m_6mFU ~ prebmi, data = cprd_dataset %>% filter(drugclass == "DPP4"), weights = prop.score, family=quasibinomial)
lm_uni_GLP1.prebmi <- glm(stopdrug_3m_6mFU ~ prebmi, data = cprd_dataset %>% filter(drugclass == "GLP1"), weights = prop.score, family=quasibinomial)
lm_uni_MFN.prebmi <- glm(stopdrug_3m_6mFU ~ prebmi, data = cprd_dataset %>% filter(drugclass == "MFN"), weights = prop.score, family=quasibinomial)
lm_uni_SGLT2.prebmi <- glm(stopdrug_3m_6mFU ~ prebmi, data = cprd_dataset %>% filter(drugclass == "SGLT2"), weights = prop.score, family=quasibinomial)
lm_uni_SU.prebmi <- glm(stopdrug_3m_6mFU ~ prebmi, data = cprd_dataset %>% filter(drugclass == "SU"), weights = prop.score, family=quasibinomial)
lm_uni_TZD.prebmi <- glm(stopdrug_3m_6mFU ~ prebmi, data = cprd_dataset %>% filter(drugclass == "TZD"), weights = prop.score, family=quasibinomial)


cprd <- cprd_dataset %>%
  mutate_at(c("prealt"), ~(scale(.) %>% as.vector))

# prealt
lm_uni_all.prealt <- glm(stopdrug_3m_6mFU ~ prealt, data = cprd_dataset, weights = prop.score, family=quasibinomial)

lm_uni_DPP4.prealt <- glm(stopdrug_3m_6mFU ~ prealt, data = cprd_dataset %>% filter(drugclass == "DPP4"), weights = prop.score, family=quasibinomial)
lm_uni_GLP1.prealt <- glm(stopdrug_3m_6mFU ~ prealt, data = cprd_dataset %>% filter(drugclass == "GLP1"), weights = prop.score, family=quasibinomial)
lm_uni_MFN.prealt <- glm(stopdrug_3m_6mFU ~ prealt, data = cprd_dataset %>% filter(drugclass == "MFN"), weights = prop.score, family=quasibinomial)
lm_uni_SGLT2.prealt <- glm(stopdrug_3m_6mFU ~ prealt, data = cprd_dataset %>% filter(drugclass == "SGLT2"), weights = prop.score, family=quasibinomial)
lm_uni_SU.prealt <- glm(stopdrug_3m_6mFU ~ prealt, data = cprd_dataset %>% filter(drugclass == "SU"), weights = prop.score, family=quasibinomial)
lm_uni_TZD.prealt <- glm(stopdrug_3m_6mFU ~ prealt, data = cprd_dataset %>% filter(drugclass == "TZD"), weights = prop.score, family=quasibinomial)



cprd <- cprd_dataset %>%
  mutate_at(c("pretotalcholesterol"), ~(scale(.) %>% as.vector))


# pretotalcholesterol
lm_uni_all.pretotalcholesterol <- glm(stopdrug_3m_6mFU ~ pretotalcholesterol, data = cprd_dataset, weights = prop.score, family=quasibinomial)

lm_uni_DPP4.pretotalcholesterol <- glm(stopdrug_3m_6mFU ~ pretotalcholesterol, data = cprd_dataset %>% filter(drugclass == "DPP4"), weights = prop.score, family=quasibinomial)
lm_uni_GLP1.pretotalcholesterol <- glm(stopdrug_3m_6mFU ~ pretotalcholesterol, data = cprd_dataset %>% filter(drugclass == "GLP1"), weights = prop.score, family=quasibinomial)
lm_uni_MFN.pretotalcholesterol <- glm(stopdrug_3m_6mFU ~ pretotalcholesterol, data = cprd_dataset %>% filter(drugclass == "MFN"), weights = prop.score, family=quasibinomial)
lm_uni_SGLT2.pretotalcholesterol <- glm(stopdrug_3m_6mFU ~ pretotalcholesterol, data = cprd_dataset %>% filter(drugclass == "SGLT2"), weights = prop.score, family=quasibinomial)
lm_uni_SU.pretotalcholesterol <- glm(stopdrug_3m_6mFU ~ pretotalcholesterol, data = cprd_dataset %>% filter(drugclass == "SU"), weights = prop.score, family=quasibinomial)
lm_uni_TZD.pretotalcholesterol <- glm(stopdrug_3m_6mFU ~ pretotalcholesterol, data = cprd_dataset %>% filter(drugclass == "TZD"), weights = prop.score, family=quasibinomial)


# gender
lm_uni_all.gender <- glm(stopdrug_3m_6mFU ~ gender, data = cprd_dataset, weights = prop.score, family=quasibinomial)

lm_uni_DPP4.gender <- glm(stopdrug_3m_6mFU ~ gender, data = cprd_dataset %>% filter(drugclass == "DPP4"), weights = prop.score, family=quasibinomial)
lm_uni_GLP1.gender <- glm(stopdrug_3m_6mFU ~ gender, data = cprd_dataset %>% filter(drugclass == "GLP1"), weights = prop.score, family=quasibinomial)
lm_uni_MFN.gender <- glm(stopdrug_3m_6mFU ~ gender, data = cprd_dataset %>% filter(drugclass == "MFN"), weights = prop.score, family=quasibinomial)
lm_uni_SGLT2.gender <- glm(stopdrug_3m_6mFU ~ gender, data = cprd_dataset %>% filter(drugclass == "SGLT2"), weights = prop.score, family=quasibinomial)
lm_uni_SU.gender <- glm(stopdrug_3m_6mFU ~ gender, data = cprd_dataset %>% filter(drugclass == "SU"), weights = prop.score, family=quasibinomial)
lm_uni_TZD.gender <- glm(stopdrug_3m_6mFU ~ gender, data = cprd_dataset %>% filter(drugclass == "TZD"), weights = prop.score, family=quasibinomial)


# predrug_bloodmed
lm_uni_all.predrug_bloodmed <- glm(stopdrug_3m_6mFU ~ predrug_bloodmed, data = cprd_dataset, weights = prop.score, family=quasibinomial)

lm_uni_DPP4.predrug_bloodmed <- glm(stopdrug_3m_6mFU ~ predrug_bloodmed, data = cprd_dataset %>% filter(drugclass == "DPP4"), weights = prop.score, family=quasibinomial)
lm_uni_GLP1.predrug_bloodmed <- glm(stopdrug_3m_6mFU ~ predrug_bloodmed, data = cprd_dataset %>% filter(drugclass == "GLP1"), weights = prop.score, family=quasibinomial)
lm_uni_MFN.predrug_bloodmed <- glm(stopdrug_3m_6mFU ~ predrug_bloodmed, data = cprd_dataset %>% filter(drugclass == "MFN"), weights = prop.score, family=quasibinomial)
lm_uni_SGLT2.predrug_bloodmed <- glm(stopdrug_3m_6mFU ~ predrug_bloodmed, data = cprd_dataset %>% filter(drugclass == "SGLT2"), weights = prop.score, family=quasibinomial)
lm_uni_SU.predrug_bloodmed <- glm(stopdrug_3m_6mFU ~ predrug_bloodmed, data = cprd_dataset %>% filter(drugclass == "SU"), weights = prop.score, family=quasibinomial)
lm_uni_TZD.predrug_bloodmed <- glm(stopdrug_3m_6mFU ~ predrug_bloodmed, data = cprd_dataset %>% filter(drugclass == "TZD"), weights = prop.score, family=quasibinomial)


# predrug_statins
lm_uni_all.predrug_statins <- glm(stopdrug_3m_6mFU ~ predrug_statins, data = cprd_dataset, weights = prop.score, family=quasibinomial)

lm_uni_DPP4.predrug_statins <- glm(stopdrug_3m_6mFU ~ predrug_statins, data = cprd_dataset %>% filter(drugclass == "DPP4"), weights = prop.score, family=quasibinomial)
lm_uni_GLP1.predrug_statins <- glm(stopdrug_3m_6mFU ~ predrug_statins, data = cprd_dataset %>% filter(drugclass == "GLP1"), weights = prop.score, family=quasibinomial)
lm_uni_MFN.predrug_statins <- glm(stopdrug_3m_6mFU ~ predrug_statins, data = cprd_dataset %>% filter(drugclass == "MFN"), weights = prop.score, family=quasibinomial)
lm_uni_SGLT2.predrug_statins <- glm(stopdrug_3m_6mFU ~ predrug_statins, data = cprd_dataset %>% filter(drugclass == "SGLT2"), weights = prop.score, family=quasibinomial)
lm_uni_SU.predrug_statins <- glm(stopdrug_3m_6mFU ~ predrug_statins, data = cprd_dataset %>% filter(drugclass == "SU"), weights = prop.score, family=quasibinomial)
lm_uni_TZD.predrug_statins <- glm(stopdrug_3m_6mFU ~ predrug_statins, data = cprd_dataset %>% filter(drugclass == "TZD"), weights = prop.score, family=quasibinomial)


# stopdrug_3m_3mFU_MFN_hist
lm_uni_all.stopdrug_3m_3mFU_MFN_hist <- glm(stopdrug_3m_6mFU ~ stopdrug_3m_3mFU_MFN_hist, data = cprd_dataset %>% filter(drugclass != "MFN"), weights = prop.score, family=quasibinomial)

lm_uni_DPP4.stopdrug_3m_3mFU_MFN_hist <- glm(stopdrug_3m_6mFU ~ stopdrug_3m_3mFU_MFN_hist, data = cprd_dataset %>% filter(drugclass == "DPP4"), weights = prop.score, family=quasibinomial)
lm_uni_GLP1.stopdrug_3m_3mFU_MFN_hist <- glm(stopdrug_3m_6mFU ~ stopdrug_3m_3mFU_MFN_hist, data = cprd_dataset %>% filter(drugclass == "GLP1"), weights = prop.score, family=quasibinomial)
# lm_uni_MFN.stopdrug_3m_3mFU_MFN_hist <- glm(stopdrug_3m_6mFU ~ stopdrug_3m_3mFU_MFN_hist, data = cprd_dataset %>% filter(drugclass == "MFN"), weights = prop.score, family=quasibinomial)
lm_uni_SGLT2.stopdrug_3m_3mFU_MFN_hist <- glm(stopdrug_3m_6mFU ~ stopdrug_3m_3mFU_MFN_hist, data = cprd_dataset %>% filter(drugclass == "SGLT2"), weights = prop.score, family=quasibinomial)
lm_uni_SU.stopdrug_3m_3mFU_MFN_hist <- glm(stopdrug_3m_6mFU ~ stopdrug_3m_3mFU_MFN_hist, data = cprd_dataset %>% filter(drugclass == "SU"), weights = prop.score, family=quasibinomial)
lm_uni_TZD.stopdrug_3m_3mFU_MFN_hist <- glm(stopdrug_3m_6mFU ~ stopdrug_3m_3mFU_MFN_hist, data = cprd_dataset %>% filter(drugclass == "TZD"), weights = prop.score, family=quasibinomial)


# predrug_cardio_event
lm_uni_all.predrug_cardio_event <- glm(stopdrug_3m_6mFU ~ predrug_cardio_event, data = cprd_dataset, weights = prop.score, family=quasibinomial)

lm_uni_DPP4.predrug_cardio_event <- glm(stopdrug_3m_6mFU ~ predrug_cardio_event, data = cprd_dataset %>% filter(drugclass == "DPP4"), weights = prop.score, family=quasibinomial)
lm_uni_GLP1.predrug_cardio_event <- glm(stopdrug_3m_6mFU ~ predrug_cardio_event, data = cprd_dataset %>% filter(drugclass == "GLP1"), weights = prop.score, family=quasibinomial)
lm_uni_MFN.predrug_cardio_event <- glm(stopdrug_3m_6mFU ~ predrug_cardio_event, data = cprd_dataset %>% filter(drugclass == "MFN"), weights = prop.score, family=quasibinomial)
lm_uni_SGLT2.predrug_cardio_event <- glm(stopdrug_3m_6mFU ~ predrug_cardio_event, data = cprd_dataset %>% filter(drugclass == "SGLT2"), weights = prop.score, family=quasibinomial)
lm_uni_SU.predrug_cardio_event <- glm(stopdrug_3m_6mFU ~ predrug_cardio_event, data = cprd_dataset %>% filter(drugclass == "SU"), weights = prop.score, family=quasibinomial)
lm_uni_TZD.predrug_cardio_event <- glm(stopdrug_3m_6mFU ~ predrug_cardio_event, data = cprd_dataset %>% filter(drugclass == "TZD"), weights = prop.score, family=quasibinomial)


# predrug_heart_event
lm_uni_all.predrug_heart_event <- glm(stopdrug_3m_6mFU ~ predrug_heart_event, data = cprd_dataset, weights = prop.score, family=quasibinomial)

lm_uni_DPP4.predrug_heart_event <- glm(stopdrug_3m_6mFU ~ predrug_heart_event, data = cprd_dataset %>% filter(drugclass == "DPP4"), weights = prop.score, family=quasibinomial)
lm_uni_GLP1.predrug_heart_event <- glm(stopdrug_3m_6mFU ~ predrug_heart_event, data = cprd_dataset %>% filter(drugclass == "GLP1"), weights = prop.score, family=quasibinomial)
lm_uni_MFN.predrug_heart_event <- glm(stopdrug_3m_6mFU ~ predrug_heart_event, data = cprd_dataset %>% filter(drugclass == "MFN"), weights = prop.score, family=quasibinomial)
lm_uni_SGLT2.predrug_heart_event <- glm(stopdrug_3m_6mFU ~ predrug_heart_event, data = cprd_dataset %>% filter(drugclass == "SGLT2"), weights = prop.score, family=quasibinomial)
lm_uni_SU.predrug_heart_event <- glm(stopdrug_3m_6mFU ~ predrug_heart_event, data = cprd_dataset %>% filter(drugclass == "SU"), weights = prop.score, family=quasibinomial)
lm_uni_TZD.predrug_heart_event <- glm(stopdrug_3m_6mFU ~ predrug_heart_event, data = cprd_dataset %>% filter(drugclass == "TZD"), weights = prop.score, family=quasibinomial)


# predrug_micro_event
lm_uni_all.predrug_micro_event <- glm(stopdrug_3m_6mFU ~ predrug_micro_event, data = cprd_dataset, weights = prop.score, family=quasibinomial)

lm_uni_DPP4.predrug_micro_event <- glm(stopdrug_3m_6mFU ~ predrug_micro_event, data = cprd_dataset %>% filter(drugclass == "DPP4"), weights = prop.score, family=quasibinomial)
lm_uni_GLP1.predrug_micro_event <- glm(stopdrug_3m_6mFU ~ predrug_micro_event, data = cprd_dataset %>% filter(drugclass == "GLP1"), weights = prop.score, family=quasibinomial)
lm_uni_MFN.predrug_micro_event <- glm(stopdrug_3m_6mFU ~ predrug_micro_event, data = cprd_dataset %>% filter(drugclass == "MFN"), weights = prop.score, family=quasibinomial)
lm_uni_SGLT2.predrug_micro_event <- glm(stopdrug_3m_6mFU ~ predrug_micro_event, data = cprd_dataset %>% filter(drugclass == "SGLT2"), weights = prop.score, family=quasibinomial)
lm_uni_SU.predrug_micro_event <- glm(stopdrug_3m_6mFU ~ predrug_micro_event, data = cprd_dataset %>% filter(drugclass == "SU"), weights = prop.score, family=quasibinomial)
lm_uni_TZD.predrug_micro_event <- glm(stopdrug_3m_6mFU ~ predrug_micro_event, data = cprd_dataset %>% filter(drugclass == "TZD"), weights = prop.score, family=quasibinomial)


# preckdstage
lm_uni_all.preckdstage <- glm(stopdrug_3m_6mFU ~ preckdstage, data = cprd_dataset, weights = prop.score, family=quasibinomial)

lm_uni_DPP4.preckdstage <- glm(stopdrug_3m_6mFU ~ preckdstage, data = cprd_dataset %>% filter(drugclass == "DPP4"), weights = prop.score, family=quasibinomial)
lm_uni_GLP1.preckdstage <- glm(stopdrug_3m_6mFU ~ preckdstage, data = cprd_dataset %>% filter(drugclass == "GLP1"), weights = prop.score, family=quasibinomial)
lm_uni_MFN.preckdstage <- glm(stopdrug_3m_6mFU ~ preckdstage, data = cprd_dataset %>% filter(drugclass == "MFN"), weights = prop.score, family=quasibinomial)
lm_uni_SGLT2.preckdstage <- glm(stopdrug_3m_6mFU ~ preckdstage, data = cprd_dataset %>% filter(drugclass == "SGLT2"), weights = prop.score, family=quasibinomial)
lm_uni_SU.preckdstage <- glm(stopdrug_3m_6mFU ~ preckdstage, data = cprd_dataset %>% filter(drugclass == "SU"), weights = prop.score, family=quasibinomial)
lm_uni_TZD.preckdstage <- glm(stopdrug_3m_6mFU ~ preckdstage, data = cprd_dataset %>% filter(drugclass == "TZD"), weights = prop.score, family=quasibinomial)


# predrug_cld
lm_uni_all.predrug_cld <- glm(stopdrug_3m_6mFU ~ predrug_cld, data = cprd_dataset, weights = prop.score, family=quasibinomial)

lm_uni_DPP4.predrug_cld <- glm(stopdrug_3m_6mFU ~ predrug_cld, data = cprd_dataset %>% filter(drugclass == "DPP4"), weights = prop.score, family=quasibinomial)
lm_uni_GLP1.predrug_cld <- glm(stopdrug_3m_6mFU ~ predrug_cld, data = cprd_dataset %>% filter(drugclass == "GLP1"), weights = prop.score, family=quasibinomial)
lm_uni_MFN.predrug_cld <- glm(stopdrug_3m_6mFU ~ predrug_cld, data = cprd_dataset %>% filter(drugclass == "MFN"), weights = prop.score, family=quasibinomial)
lm_uni_SGLT2.predrug_cld <- glm(stopdrug_3m_6mFU ~ predrug_cld, data = cprd_dataset %>% filter(drugclass == "SGLT2"), weights = prop.score, family=quasibinomial)
lm_uni_SU.predrug_cld <- glm(stopdrug_3m_6mFU ~ predrug_cld, data = cprd_dataset %>% filter(drugclass == "SU"), weights = prop.score, family=quasibinomial)
lm_uni_TZD.predrug_cld <- glm(stopdrug_3m_6mFU ~ predrug_cld, data = cprd_dataset %>% filter(drugclass == "TZD"), family=quasibinomial)




## Plot all results into one table
drugs = c("Pooled", "DPP4", "GLP1", "MFN", "SGLT2", "SU", "TZD")
coefficients <- data.frame(
  cbind(
    variable = rep("Diabetes duration (per 10 years)", 7),
    type = rep("Biological", 7),
    drug = drugs,
    coef = c(lm_uni_all.dstartdate_dm_dur$coef[2], lm_uni_DPP4.dstartdate_dm_dur$coef[2], lm_uni_GLP1.dstartdate_dm_dur$coef[2] , lm_uni_MFN.dstartdate_dm_dur$coef[2] , lm_uni_SGLT2.dstartdate_dm_dur$coef[2] , lm_uni_SU.dstartdate_dm_dur$coef[2] , lm_uni_TZD.dstartdate_dm_dur$coef[2]),
    LCI = c(confint(lm_uni_all.dstartdate_dm_dur)[2,1], confint(lm_uni_DPP4.dstartdate_dm_dur)[2,1], confint(lm_uni_GLP1.dstartdate_dm_dur)[2,1] , confint(lm_uni_MFN.dstartdate_dm_dur)[2,1] , confint(lm_uni_SGLT2.dstartdate_dm_dur)[2,1] , confint(lm_uni_SU.dstartdate_dm_dur)[2,1] , confint(lm_uni_TZD.dstartdate_dm_dur)[2,1]),
    UCI = c(confint(lm_uni_all.dstartdate_dm_dur)[2,2], confint(lm_uni_DPP4.dstartdate_dm_dur)[2,2], confint(lm_uni_GLP1.dstartdate_dm_dur)[2,2] , confint(lm_uni_MFN.dstartdate_dm_dur)[2,2] , confint(lm_uni_SGLT2.dstartdate_dm_dur)[2,2] , confint(lm_uni_SU.dstartdate_dm_dur)[2,2] , confint(lm_uni_TZD.dstartdate_dm_dur)[2,2])
  )
) %>%
  rbind(
    data.frame(
      cbind(
        variable = rep("Age (per 10 years)", 7),
        type = rep("Biological", 7),
        drug = drugs,
        coef = c(lm_uni_all.dstartdate_age$coef[2], lm_uni_DPP4.dstartdate_age$coef[2], lm_uni_GLP1.dstartdate_age$coef[2] , lm_uni_MFN.dstartdate_age$coef[2] , lm_uni_SGLT2.dstartdate_age$coef[2] , lm_uni_SU.dstartdate_age$coef[2] , lm_uni_TZD.dstartdate_age$coef[2]),
        LCI = c(confint(lm_uni_all.dstartdate_age)[2,1], confint(lm_uni_DPP4.dstartdate_age)[2,1], confint(lm_uni_GLP1.dstartdate_age)[2,1] , confint(lm_uni_MFN.dstartdate_age)[2,1] , confint(lm_uni_SGLT2.dstartdate_age)[2,1] , confint(lm_uni_SU.dstartdate_age)[2,1] , confint(lm_uni_TZD.dstartdate_age)[2,1]),
        UCI = c(confint(lm_uni_all.dstartdate_age)[2,2], confint(lm_uni_DPP4.dstartdate_age)[2,2], confint(lm_uni_GLP1.dstartdate_age)[2,2] , confint(lm_uni_MFN.dstartdate_age)[2,2] , confint(lm_uni_SGLT2.dstartdate_age)[2,2] , confint(lm_uni_SU.dstartdate_age)[2,2] , confint(lm_uni_TZD.dstartdate_age)[2,2])
      )
    ),
    data.frame(
      cbind(
        variable = rep("HbA1c (per 10mmol/mol)", 7),
        type = rep("Biological", 7),
        drug = drugs,
        coef = c(lm_uni_all.prehba1c$coef[2], lm_uni_DPP4.prehba1c$coef[2], lm_uni_GLP1.prehba1c$coef[2] , lm_uni_MFN.prehba1c$coef[2] , lm_uni_SGLT2.prehba1c$coef[2] , lm_uni_SU.prehba1c$coef[2] , lm_uni_TZD.prehba1c$coef[2]),
        LCI = c(confint(lm_uni_all.prehba1c)[2,1], confint(lm_uni_DPP4.prehba1c)[2,1], confint(lm_uni_GLP1.prehba1c)[2,1] , confint(lm_uni_MFN.prehba1c)[2,1] , confint(lm_uni_SGLT2.prehba1c)[2,1] , confint(lm_uni_SU.prehba1c)[2,1] , confint(lm_uni_TZD.prehba1c)[2,1]),
        UCI = c(confint(lm_uni_all.prehba1c)[2,2], confint(lm_uni_DPP4.prehba1c)[2,2], confint(lm_uni_GLP1.prehba1c)[2,2] , confint(lm_uni_MFN.prehba1c)[2,2] , confint(lm_uni_SGLT2.prehba1c)[2,2] , confint(lm_uni_SU.prehba1c)[2,2] , confint(lm_uni_TZD.prehba1c)[2,2])
      )
    ),
    data.frame(
      cbind(
        variable = rep("eGFR (per SD)", 7),
        type = rep("Biological", 7),
        drug = drugs,
        coef = c(lm_uni_all.preegfr$coef[2], lm_uni_DPP4.preegfr$coef[2], lm_uni_GLP1.preegfr$coef[2] , lm_uni_MFN.preegfr$coef[2] , lm_uni_SGLT2.preegfr$coef[2] , lm_uni_SU.preegfr$coef[2] , lm_uni_TZD.preegfr$coef[2]),
        LCI = c(confint(lm_uni_all.preegfr)[2,1], confint(lm_uni_DPP4.preegfr)[2,1], confint(lm_uni_GLP1.preegfr)[2,1] , confint(lm_uni_MFN.preegfr)[2,1] , confint(lm_uni_SGLT2.preegfr)[2,1] , confint(lm_uni_SU.preegfr)[2,1] , confint(lm_uni_TZD.preegfr)[2,1]),
        UCI = c(confint(lm_uni_all.preegfr)[2,2], confint(lm_uni_DPP4.preegfr)[2,2], confint(lm_uni_GLP1.preegfr)[2,2] , confint(lm_uni_MFN.preegfr)[2,2] , confint(lm_uni_SGLT2.preegfr)[2,2] , confint(lm_uni_SU.preegfr)[2,2] , confint(lm_uni_TZD.preegfr)[2,2])
      )
    ),
    data.frame(
      cbind(
        variable = rep("BMI (per 5kg/m^2)", 7),
        type = rep("Biological", 7),
        drug = drugs,
        coef = c(lm_uni_all.prebmi$coef[2], lm_uni_DPP4.prebmi$coef[2], lm_uni_GLP1.prebmi$coef[2] , lm_uni_MFN.prebmi$coef[2] , lm_uni_SGLT2.prebmi$coef[2] , lm_uni_SU.prebmi$coef[2] , lm_uni_TZD.prebmi$coef[2]),
        LCI = c(confint(lm_uni_all.prebmi)[2,1], confint(lm_uni_DPP4.prebmi)[2,1], confint(lm_uni_GLP1.prebmi)[2,1] , confint(lm_uni_MFN.prebmi)[2,1] , confint(lm_uni_SGLT2.prebmi)[2,1] , confint(lm_uni_SU.prebmi)[2,1] , confint(lm_uni_TZD.prebmi)[2,1]),
        UCI = c(confint(lm_uni_all.prebmi)[2,2], confint(lm_uni_DPP4.prebmi)[2,2], confint(lm_uni_GLP1.prebmi)[2,2] , confint(lm_uni_MFN.prebmi)[2,2] , confint(lm_uni_SGLT2.prebmi)[2,2] , confint(lm_uni_SU.prebmi)[2,2] , confint(lm_uni_TZD.prebmi)[2,2])
      )
    ),
    data.frame(
      cbind(
        variable = rep("ALT (per SD)", 7),
        type = rep("Biological", 7),
        drug = drugs,
        coef = c(lm_uni_all.prealt$coef[2], lm_uni_DPP4.prealt$coef[2], lm_uni_GLP1.prealt$coef[2] , lm_uni_MFN.prealt$coef[2] , lm_uni_SGLT2.prealt$coef[2] , lm_uni_SU.prealt$coef[2] , lm_uni_TZD.prealt$coef[2]),
        LCI = c(confint(lm_uni_all.prealt)[2,1], confint(lm_uni_DPP4.prealt)[2,1], confint(lm_uni_GLP1.prealt)[2,1] , confint(lm_uni_MFN.prealt)[2,1] , confint(lm_uni_SGLT2.prealt)[2,1] , confint(lm_uni_SU.prealt)[2,1] , confint(lm_uni_TZD.prealt)[2,1]),
        UCI = c(confint(lm_uni_all.prealt)[2,2], confint(lm_uni_DPP4.prealt)[2,2], confint(lm_uni_GLP1.prealt)[2,2] , confint(lm_uni_MFN.prealt)[2,2] , confint(lm_uni_SGLT2.prealt)[2,2] , confint(lm_uni_SU.prealt)[2,2] , confint(lm_uni_TZD.prealt)[2,2])
      )
    ),
    data.frame(
      cbind(
        variable = rep("Total Cholesterol (per SD)", 7),
        type = rep("Biological", 7),
        drug = drugs,
        coef = c(lm_uni_all.pretotalcholesterol$coef[2], lm_uni_DPP4.pretotalcholesterol$coef[2], lm_uni_GLP1.pretotalcholesterol$coef[2] , lm_uni_MFN.pretotalcholesterol$coef[2] , lm_uni_SGLT2.pretotalcholesterol$coef[2] , lm_uni_SU.pretotalcholesterol$coef[2] , lm_uni_TZD.pretotalcholesterol$coef[2]),
        LCI = c(confint(lm_uni_all.pretotalcholesterol)[2,1], confint(lm_uni_DPP4.pretotalcholesterol)[2,1], confint(lm_uni_GLP1.pretotalcholesterol)[2,1] , confint(lm_uni_MFN.pretotalcholesterol)[2,1] , confint(lm_uni_SGLT2.pretotalcholesterol)[2,1] , confint(lm_uni_SU.pretotalcholesterol)[2,1] , confint(lm_uni_TZD.pretotalcholesterol)[2,1]),
        UCI = c(confint(lm_uni_all.pretotalcholesterol)[2,2], confint(lm_uni_DPP4.pretotalcholesterol)[2,2], confint(lm_uni_GLP1.pretotalcholesterol)[2,2] , confint(lm_uni_MFN.pretotalcholesterol)[2,2] , confint(lm_uni_SGLT2.pretotalcholesterol)[2,2] , confint(lm_uni_SU.pretotalcholesterol)[2,2] , confint(lm_uni_TZD.pretotalcholesterol)[2,2])
      )
    ),
    data.frame(
      cbind(
        variable = rep("Sex (ref Male)", 7),
        type = rep("Biological", 7),
        drug = drugs,
        coef = c(lm_uni_all.gender$coef[2], lm_uni_DPP4.gender$coef[2], lm_uni_GLP1.gender$coef[2] , lm_uni_MFN.gender$coef[2] , lm_uni_SGLT2.gender$coef[2] , lm_uni_SU.gender$coef[2] , lm_uni_TZD.gender$coef[2]),
        LCI = c(confint(lm_uni_all.gender)[2,1], confint(lm_uni_DPP4.gender)[2,1], confint(lm_uni_GLP1.gender)[2,1] , confint(lm_uni_MFN.gender)[2,1] , confint(lm_uni_SGLT2.gender)[2,1] , confint(lm_uni_SU.gender)[2,1] , confint(lm_uni_TZD.gender)[2,1]),
        UCI = c(confint(lm_uni_all.gender)[2,2], confint(lm_uni_DPP4.gender)[2,2], confint(lm_uni_GLP1.gender)[2,2] , confint(lm_uni_MFN.gender)[2,2] , confint(lm_uni_SGLT2.gender)[2,2] , confint(lm_uni_SU.gender)[2,2] , confint(lm_uni_TZD.gender)[2,2])
      )
    ),
    data.frame(
      cbind(
        variable = rep("Blood medication", 7),
        type = rep("Behavioural", 7),
        drug = drugs,
        coef = c(lm_uni_all.predrug_bloodmed$coef[2], lm_uni_DPP4.predrug_bloodmed$coef[2], lm_uni_GLP1.predrug_bloodmed$coef[2] , lm_uni_MFN.predrug_bloodmed$coef[2] , lm_uni_SGLT2.predrug_bloodmed$coef[2] , lm_uni_SU.predrug_bloodmed$coef[2] , lm_uni_TZD.predrug_bloodmed$coef[2]),
        LCI = c(confint(lm_uni_all.predrug_bloodmed)[2,1], confint(lm_uni_DPP4.predrug_bloodmed)[2,1], confint(lm_uni_GLP1.predrug_bloodmed)[2,1] , confint(lm_uni_MFN.predrug_bloodmed)[2,1] , confint(lm_uni_SGLT2.predrug_bloodmed)[2,1] , confint(lm_uni_SU.predrug_bloodmed)[2,1] , confint(lm_uni_TZD.predrug_bloodmed)[2,1]),
        UCI = c(confint(lm_uni_all.predrug_bloodmed)[2,2], confint(lm_uni_DPP4.predrug_bloodmed)[2,2], confint(lm_uni_GLP1.predrug_bloodmed)[2,2] , confint(lm_uni_MFN.predrug_bloodmed)[2,2] , confint(lm_uni_SGLT2.predrug_bloodmed)[2,2] , confint(lm_uni_SU.predrug_bloodmed)[2,2] , confint(lm_uni_TZD.predrug_bloodmed)[2,2])
      )
    ),
    data.frame(
      cbind(
        variable = rep("Statins", 7),
        type = rep("Behavioural", 7),
        drug = drugs,
        coef = c(lm_uni_all.predrug_statins$coef[2], lm_uni_DPP4.predrug_statins$coef[2], lm_uni_GLP1.predrug_statins$coef[2] , lm_uni_MFN.predrug_statins$coef[2] , lm_uni_SGLT2.predrug_statins$coef[2] , lm_uni_SU.predrug_statins$coef[2] , lm_uni_TZD.predrug_statins$coef[2]),
        LCI = c(confint(lm_uni_all.predrug_statins)[2,1], confint(lm_uni_DPP4.predrug_statins)[2,1], confint(lm_uni_GLP1.predrug_statins)[2,1] , confint(lm_uni_MFN.predrug_statins)[2,1] , confint(lm_uni_SGLT2.predrug_statins)[2,1] , confint(lm_uni_SU.predrug_statins)[2,1] , confint(lm_uni_TZD.predrug_statins)[2,1]),
        UCI = c(confint(lm_uni_all.predrug_statins)[2,2], confint(lm_uni_DPP4.predrug_statins)[2,2], confint(lm_uni_GLP1.predrug_statins)[2,2] , confint(lm_uni_MFN.predrug_statins)[2,2] , confint(lm_uni_SGLT2.predrug_statins)[2,2] , confint(lm_uni_SU.predrug_statins)[2,2] , confint(lm_uni_TZD.predrug_statins)[2,2])
      )
    ),
    data.frame(
      cbind(
        variable = rep("MFN discontinuation", 6),
        type = rep("Behavioural", 6),
        drug = c("Pooled", "DPP4", "GLP1", "SGLT2", "SU", "TZD"),
        coef = c(lm_uni_all.stopdrug_3m_3mFU_MFN_hist$coef[2], lm_uni_DPP4.stopdrug_3m_3mFU_MFN_hist$coef[2], lm_uni_GLP1.stopdrug_3m_3mFU_MFN_hist$coef[2] , lm_uni_SGLT2.stopdrug_3m_3mFU_MFN_hist$coef[2] , lm_uni_SU.stopdrug_3m_3mFU_MFN_hist$coef[2] , lm_uni_TZD.stopdrug_3m_3mFU_MFN_hist$coef[2]),
        LCI = c(confint(lm_uni_all.stopdrug_3m_3mFU_MFN_hist)[2,1], confint(lm_uni_DPP4.stopdrug_3m_3mFU_MFN_hist)[2,1], confint(lm_uni_GLP1.stopdrug_3m_3mFU_MFN_hist)[2,1] , confint(lm_uni_SGLT2.stopdrug_3m_3mFU_MFN_hist)[2,1] , confint(lm_uni_SU.stopdrug_3m_3mFU_MFN_hist)[2,1] , confint(lm_uni_TZD.stopdrug_3m_3mFU_MFN_hist)[2,1]),
        UCI = c(confint(lm_uni_all.stopdrug_3m_3mFU_MFN_hist)[2,2], confint(lm_uni_DPP4.stopdrug_3m_3mFU_MFN_hist)[2,2], confint(lm_uni_GLP1.stopdrug_3m_3mFU_MFN_hist)[2,2] ,  confint(lm_uni_SGLT2.stopdrug_3m_3mFU_MFN_hist)[2,2] , confint(lm_uni_SU.stopdrug_3m_3mFU_MFN_hist)[2,2] , confint(lm_uni_TZD.stopdrug_3m_3mFU_MFN_hist)[2,2])
      )
    ),
    data.frame(
      cbind(
        variable = rep("Cardiovascular event", 7),
        type = rep("Biological", 7),
        drug = drugs,
        coef = c(lm_uni_all.predrug_cardio_event$coef[2], lm_uni_DPP4.predrug_cardio_event$coef[2], lm_uni_GLP1.predrug_cardio_event$coef[2] , lm_uni_MFN.predrug_cardio_event$coef[2] , lm_uni_SGLT2.predrug_cardio_event$coef[2] , lm_uni_SU.predrug_cardio_event$coef[2] , lm_uni_TZD.predrug_cardio_event$coef[2]),
        LCI = c(confint(lm_uni_all.predrug_cardio_event)[2,1], confint(lm_uni_DPP4.predrug_cardio_event)[2,1], confint(lm_uni_GLP1.predrug_cardio_event)[2,1] , confint(lm_uni_MFN.predrug_cardio_event)[2,1] , confint(lm_uni_SGLT2.predrug_cardio_event)[2,1] , confint(lm_uni_SU.predrug_cardio_event)[2,1] , confint(lm_uni_TZD.predrug_cardio_event)[2,1]),
        UCI = c(confint(lm_uni_all.predrug_cardio_event)[2,2], confint(lm_uni_DPP4.predrug_cardio_event)[2,2], confint(lm_uni_GLP1.predrug_cardio_event)[2,2] , confint(lm_uni_MFN.predrug_cardio_event)[2,2] , confint(lm_uni_SGLT2.predrug_cardio_event)[2,2] , confint(lm_uni_SU.predrug_cardio_event)[2,2] , confint(lm_uni_TZD.predrug_cardio_event)[2,2])
      )
    ),
    data.frame(
      cbind(
        variable = rep("Heart event", 7),
        type = rep("Biological", 7),
        drug = drugs,
        coef = c(lm_uni_all.predrug_heart_event$coef[2], lm_uni_DPP4.predrug_heart_event$coef[2], lm_uni_GLP1.predrug_heart_event$coef[2] , lm_uni_MFN.predrug_heart_event$coef[2] , lm_uni_SGLT2.predrug_heart_event$coef[2] , lm_uni_SU.predrug_heart_event$coef[2] , lm_uni_TZD.predrug_heart_event$coef[2]),
        LCI = c(confint(lm_uni_all.predrug_heart_event)[2,1], confint(lm_uni_DPP4.predrug_heart_event)[2,1], confint(lm_uni_GLP1.predrug_heart_event)[2,1] , confint(lm_uni_MFN.predrug_heart_event)[2,1] , confint(lm_uni_SGLT2.predrug_heart_event)[2,1] , confint(lm_uni_SU.predrug_heart_event)[2,1] , confint(lm_uni_TZD.predrug_heart_event)[2,1]),
        UCI = c(confint(lm_uni_all.predrug_heart_event)[2,2], confint(lm_uni_DPP4.predrug_heart_event)[2,2], confint(lm_uni_GLP1.predrug_heart_event)[2,2] , confint(lm_uni_MFN.predrug_heart_event)[2,2] , confint(lm_uni_SGLT2.predrug_heart_event)[2,2] , confint(lm_uni_SU.predrug_heart_event)[2,2] , confint(lm_uni_TZD.predrug_heart_event)[2,2])
      )
    ),
    data.frame(
      cbind(
        variable = rep("Microvascular event", 7),
        type = rep("Biological", 7),
        drug = drugs,
        coef = c(lm_uni_all.predrug_micro_event$coef[2], lm_uni_DPP4.predrug_micro_event$coef[2], lm_uni_GLP1.predrug_micro_event$coef[2] , lm_uni_MFN.predrug_micro_event$coef[2] , lm_uni_SGLT2.predrug_micro_event$coef[2] , lm_uni_SU.predrug_micro_event$coef[2] , lm_uni_TZD.predrug_micro_event$coef[2]),
        LCI = c(confint(lm_uni_all.predrug_micro_event)[2,1], confint(lm_uni_DPP4.predrug_micro_event)[2,1], confint(lm_uni_GLP1.predrug_micro_event)[2,1] , confint(lm_uni_MFN.predrug_micro_event)[2,1] , confint(lm_uni_SGLT2.predrug_micro_event)[2,1] , confint(lm_uni_SU.predrug_micro_event)[2,1] , confint(lm_uni_TZD.predrug_micro_event)[2,1]),
        UCI = c(confint(lm_uni_all.predrug_micro_event)[2,2], confint(lm_uni_DPP4.predrug_micro_event)[2,2], confint(lm_uni_GLP1.predrug_micro_event)[2,2] , confint(lm_uni_MFN.predrug_micro_event)[2,2] , confint(lm_uni_SGLT2.predrug_micro_event)[2,2] , confint(lm_uni_SU.predrug_micro_event)[2,2] , confint(lm_uni_TZD.predrug_micro_event)[2,2])
      )
    ),
    data.frame(
      cbind(
        variable = rep("CKD Stage 1 (ref Stage 0)", 7),
        type = rep("Biological", 7),
        drug = drugs,
        coef = c(lm_uni_all.preckdstage$coef[2], lm_uni_DPP4.preckdstage$coef[2], lm_uni_GLP1.preckdstage$coef[2] , lm_uni_MFN.preckdstage$coef[2] , lm_uni_SGLT2.preckdstage$coef[2] , lm_uni_SU.preckdstage$coef[2] , lm_uni_TZD.preckdstage$coef[2]),
        LCI = c(confint(lm_uni_all.preckdstage)[2,1], confint(lm_uni_DPP4.preckdstage)[2,1], confint(lm_uni_GLP1.preckdstage)[2,1] , confint(lm_uni_MFN.preckdstage)[2,1] , confint(lm_uni_SGLT2.preckdstage)[2,1] , confint(lm_uni_SU.preckdstage)[2,1] , confint(lm_uni_TZD.preckdstage)[2,1]),
        UCI = c(confint(lm_uni_all.preckdstage)[2,2], confint(lm_uni_DPP4.preckdstage)[2,2], confint(lm_uni_GLP1.preckdstage)[2,2] , confint(lm_uni_MFN.preckdstage)[2,2] , confint(lm_uni_SGLT2.preckdstage)[2,2] , confint(lm_uni_SU.preckdstage)[2,2] , confint(lm_uni_TZD.preckdstage)[2,2])
      )
    ),
    data.frame(
      cbind(
        variable = rep("CKD Stage 2 (ref Stage 0)", 7),
        type = rep("Biological", 7),
        drug = drugs,
        coef = c(lm_uni_all.preckdstage$coef[3], lm_uni_DPP4.preckdstage$coef[3], lm_uni_GLP1.preckdstage$coef[3] , lm_uni_MFN.preckdstage$coef[3] , lm_uni_SGLT2.preckdstage$coef[3] , lm_uni_SU.preckdstage$coef[3] , lm_uni_TZD.preckdstage$coef[3]),
        LCI = c(confint(lm_uni_all.preckdstage)[3,1], confint(lm_uni_DPP4.preckdstage)[3,1], confint(lm_uni_GLP1.preckdstage)[3,1] , confint(lm_uni_MFN.preckdstage)[3,1] , confint(lm_uni_SGLT2.preckdstage)[3,1] , confint(lm_uni_SU.preckdstage)[3,1] , confint(lm_uni_TZD.preckdstage)[3,1]),
        UCI = c(confint(lm_uni_all.preckdstage)[3,2], confint(lm_uni_DPP4.preckdstage)[3,2], confint(lm_uni_GLP1.preckdstage)[3,2] , confint(lm_uni_MFN.preckdstage)[3,2] , confint(lm_uni_SGLT2.preckdstage)[3,2] , confint(lm_uni_SU.preckdstage)[3,2] , confint(lm_uni_TZD.preckdstage)[3,2])
      )
    ),
    data.frame(
      cbind(
        variable = rep("CLD", 7),
        type = rep("Biological", 7),
        drug = drugs,
        coef = c(lm_uni_all.predrug_cld$coef[2], lm_uni_DPP4.predrug_cld$coef[2], lm_uni_GLP1.predrug_cld$coef[2] , lm_uni_MFN.predrug_cld$coef[2] , lm_uni_SGLT2.predrug_cld$coef[2] , lm_uni_SU.predrug_cld$coef[2] , lm_uni_TZD.predrug_cld$coef[2]),
        LCI = c(confint(lm_uni_all.predrug_cld)[2,1], confint(lm_uni_DPP4.predrug_cld)[2,1], confint(lm_uni_GLP1.predrug_cld)[2,1] , confint(lm_uni_MFN.predrug_cld)[2,1] , confint(lm_uni_SGLT2.predrug_cld)[2,1] , confint(lm_uni_SU.predrug_cld)[2,1] , confint(lm_uni_TZD.predrug_cld)[2,1]),
        UCI = c(confint(lm_uni_all.predrug_cld)[2,2], confint(lm_uni_DPP4.predrug_cld)[2,2], confint(lm_uni_GLP1.predrug_cld)[2,2] , confint(lm_uni_MFN.predrug_cld)[2,2] , confint(lm_uni_SGLT2.predrug_cld)[2,2] , confint(lm_uni_SU.predrug_cld)[2,2] , confint(lm_uni_TZD.predrug_cld)[2,2])
      )
    )
  ) %>%
  mutate(
    variable = factor(variable),
    drug = factor(drug, levels = rev(drugs)),
    type = factor(type, levels = c("Biological", "Behavioural")),
    coef = exp(as.numeric(coef)),
    LCI = exp(as.numeric(LCI)),
    UCI = exp(as.numeric(UCI))
  )


pdf("results/figures/univariate_analysis.pdf", width = 14, height = 10)

wrap_plots(
  
  coefficients %>%
    filter(type == "Behavioural") %>%
    ggplot() + 
    geom_vline(aes(xintercept = 1), colour = "black") +
    geom_errorbar(aes(x = coef, xmin = LCI, xmax = UCI, y = variable, colour = drug),position = position_dodge(width = 1)) +
    geom_pointrange(aes(x = coef, xmin = LCI, xmax = UCI, y = variable, colour = drug), position = position_dodge(width = 1)) +
    xlab("Odds Ratio") +
    coord_cartesian(xlim = c(0, 2.5)) +
    facet_grid(variable~., scales = "free_y") +
    scale_colour_manual(values = c("Pooled" = "black", "SGLT2" = "#E69F00", "GLP1" = "#56B4E9", "SU" = "#CC79A7", "DPP4" = "#0072B2", "TZD" = "#D55E00", "MFN" = "grey"), name = "Therapy", guide = guide_legend(reverse = TRUE)) +
    theme_bw() +
    theme(
      axis.title.y = element_blank(),
      legend.position = "bottom",
      strip.text = element_blank()
    ),
  
  
  coefficients %>%
    filter(type == "Biological") %>%
    ggplot() + 
    geom_vline(aes(xintercept = 1), colour = "black") +
    geom_errorbar(aes(x = coef, xmin = LCI, xmax = UCI, y = variable, colour = drug),position = position_dodge(width = 1)) +
    geom_pointrange(aes(x = coef, xmin = LCI, xmax = UCI, y = variable, colour = drug), position = position_dodge(width = 1)) +
    xlab("Odds Ratio") +
    coord_cartesian(xlim = c(0, 2.5)) +
    facet_grid(variable~., scales = "free_y") +
    
    scale_colour_manual(values = c("Pooled" = "black", "SGLT2" = "#E69F00", "GLP1" = "#56B4E9", "SU" = "#CC79A7", "DPP4" = "#0072B2", "TZD" = "#D55E00", "MFN" = "grey"), name = "Therapy", guide = guide_legend(reverse = TRUE)) +
    theme_bw() +
    theme(
      axis.title.y = element_blank(),
      legend.position = "bottom",
      strip.text = element_blank()
    ),
  
  ncol = 1, nrow = 2
  
) + 
  plot_layout(guides = "collect", height = c(3, 14), axis_titles = "collect", axes = "collect") +
  plot_annotation(tag_levels = list(c('Behavioural', "Biological"))) & 
  theme(
    legend.direction = "vertical",
    axis.text.y = element_text(size = 15),
    axis.text.x = element_text(size = 15),
    axis.title.x = element_text(size = 18)
  )

dev.off()


