####################
## Description: 
##  - In this file we:
##    - Investigate heterogeneity of therapy discontinuation
#################### 


# load working directory
setwd("Samples/T2D_Discontinuation")

# load functions
source("code/00.set_up_data_and_functions.R")

# load libraries
library(tidyverse)
library(MatchIt)
library(qpdf)
library(unihtee)
library(data.table)

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
  dataset = "3m.disc.dataset.dev",
  follow_up = "6-months",
  full_prescribing_history = TRUE
) %>%
  drop_na(-stopdrug_6m_6mFU, -stopdrug_12m_6mFU)


# load dataset
cprd_dataset.val <- set_up_data(
  raw_data = "20240308_t2d_1stinstance",
  diagnosis = FALSE,
  therapies = c("DPP4", "GLP1", "MFN", "SGLT2", "SU", "TZD"),
  dataset = "3m.disc.dataset.val",
  follow_up = "6-months",
  full_prescribing_history = TRUE
) %>%
  drop_na(-stopdrug_6m_6mFU, -stopdrug_12m_6mFU)



# load propensity scores
ps.only_dataset <- readRDS("results/PS_model/ps.dataset_lm_all.rds")


# join propensity scores into dataset
cprd_dataset.dev <- cprd_dataset.dev %>%
  left_join(
    ps.only_dataset %>%
      select(patid, dstartdate, prop.score.overlap, prop.score.IPW),
    by = c("patid", "dstartdate")
  )
cprd_dataset.val <- cprd_dataset.val %>%
  left_join(
    ps.only_dataset %>%
      select(patid, dstartdate, prop.score.overlap, prop.score.IPW),
    by = c("patid", "dstartdate")
  )


probabilities.only_dataset <- readRDS("results/Models/Predictions/04.model_predictions_3m_all_drugs.rds")


# join predictions into dataset
cprd_dataset.dev <- cprd_dataset.dev %>%
  left_join(
    probabilities.only_dataset,
    by = c("patid", "dstartdate")
  )
cprd_dataset.val <- cprd_dataset.val %>%
  left_join(
    probabilities.only_dataset,
    by = c("patid", "dstartdate")
  )



###############################################################################
###############################################################################
######################### Testing heterogeneity ###############################
###############################################################################
###############################################################################

## set up dataset
variables_needed <- c(
  # outcome
  "stopdrug_3m_6mFU",
  # treatment
  "drugclass",
  # variables
  "dstartdate_dm_dur", "dstartdate_age", "drugline", "numdrugs", "smoking_cat",
  "imd2015_10", "predrug_statins", "stopdrug_3m_3mFU_MFN_hist", "ethnicity_5cat",
  "gender", "predrug_bloodmed", "prehba1c", "preegfr", "prebmi", 
  # "prehdl", 
  # "predrug_cardio_event", "predrug_angina", "predrug_myocardialinfarction", 
  # "predrug_ihd", "predrug_pad", "predrug_revasc", "predrug_stroke", 
  # "predrug_heart_event", "predrug_heartfailure", "predrug_hypertension",
  # "predrug_micro_event", "predrug_retinopathy", "predrug_diabeticnephropathy",
  # "predrug_neuropathy", "preckdstage", "predrug_cld", 
  "predrug_frailty_proxy"
)

variables_needed_cat <- cprd_dataset.dev %>%
  select(all_of(variables_needed)) %>%
  select_if(is.factor) %>% 
  colnames()

#:---- GLP1 vs DPP4
heterogeneity_glp1_dpp4 <- cprd_dataset.dev %>%
  select(all_of(variables_needed)) %>%
  filter(drugclass %in% c("GLP1", "DPP4")) %>%
  mutate(drugclass = factor(drugclass)) %>%
  mutate_all(as.numeric) %>% mutate_at(all_of(variables_needed_cat), ~ . - 1)

## targeted maximum likelihood estimates and testing procedure
heterogeneity_glp1_dpp4_p_values <- unihtee(
  data = heterogeneity_glp1_dpp4 %>% as.data.table(),
  confounders = colnames(heterogeneity_glp1_dpp4)[-c(1,2)],
  modifiers = colnames(heterogeneity_glp1_dpp4)[-c(1,2)],
  exposure = "drugclass", outcome = "stopdrug_3m_6mFU", outcome_type = "binary", effect = "relative", estimator = "tmle"
)

#:---- GLP1 vs SU
heterogeneity_glp1_su <- cprd_dataset.dev %>%
  select(all_of(variables_needed)) %>%
  filter(drugclass %in% c("GLP1", "SU")) %>%
  mutate(drugclass = factor(drugclass)) %>%
  mutate_all(as.numeric) %>% mutate_at(all_of(variables_needed_cat), ~ . - 1)

## targeted maximum likelihood estimates and testing procedure
heterogeneity_glp1_su_p_values <- unihtee(
  data = heterogeneity_glp1_su %>% as.matrix() %>% as.data.table(),
  confounders = colnames(heterogeneity_glp1_su)[-c(1,2)],
  modifiers = colnames(heterogeneity_glp1_su)[-c(1,2)],
  exposure = "drugclass", outcome = "stopdrug_3m_6mFU", outcome_type = "binary", effect = "relative", estimator = "tmle"
)

#:---- GLP1 vs TZD
heterogeneity_glp1_tzd <- cprd_dataset.dev %>%
  select(all_of(variables_needed)) %>%
  filter(drugclass %in% c("GLP1", "TZD")) %>%
  mutate(drugclass = factor(drugclass)) %>%
  mutate_all(as.numeric) %>% mutate_at(all_of(variables_needed_cat), ~ . - 1)

## targeted maximum likelihood estimates and testing procedure
heterogeneity_glp1_tzd_p_values <- unihtee(
  data = heterogeneity_glp1_tzd %>% as.matrix() %>% as.data.table(),
  confounders = colnames(heterogeneity_glp1_tzd)[-c(1,2)],
  modifiers = colnames(heterogeneity_glp1_tzd)[-c(1,2)],
  exposure = "drugclass", outcome = "stopdrug_3m_6mFU", outcome_type = "binary", effect = "relative", estimator = "tmle"
)

#:---- SGLT2 vs DPP4
heterogeneity_sglt2_dpp4 <- cprd_dataset.dev %>%
  select(all_of(variables_needed)) %>%
  filter(drugclass %in% c("SGLT2", "DPP4")) %>%
  mutate(drugclass = factor(drugclass)) %>%
  mutate_all(as.numeric) %>% mutate_at(all_of(variables_needed_cat), ~ . - 1)

## targeted maximum likelihood estimates and testing procedure
heterogeneity_sglt2_dpp4_p_values <- unihtee(
  data = heterogeneity_sglt2_dpp4 %>% as.matrix() %>% as.data.table(),
  confounders = colnames(heterogeneity_sglt2_dpp4)[-c(1,2)],
  modifiers = colnames(heterogeneity_sglt2_dpp4)[-c(1,2)],
  exposure = "drugclass", outcome = "stopdrug_3m_6mFU", outcome_type = "binary", effect = "relative", estimator = "tmle"
)

#:---- SGLT2 vs GLP1
heterogeneity_sglt2_glp1 <- cprd_dataset.dev %>%
  select(all_of(variables_needed)) %>%
  filter(drugclass %in% c("SGLT2", "GLP1")) %>%
  mutate(drugclass = factor(drugclass)) %>%
  mutate_all(as.numeric) %>% mutate_at(all_of(variables_needed_cat), ~ . - 1)

## targeted maximum likelihood estimates and testing procedure
heterogeneity_sglt2_glp1_p_values <- unihtee(
  data = heterogeneity_sglt2_glp1 %>% as.matrix() %>% as.data.table(),
  confounders = colnames(heterogeneity_sglt2_glp1)[-c(1,2)],
  modifiers = colnames(heterogeneity_sglt2_glp1)[-c(1,2)],
  exposure = "drugclass", outcome = "stopdrug_3m_6mFU", outcome_type = "binary", effect = "relative", estimator = "tmle"
)

#:---- SGLT2 vs SU
heterogeneity_sglt2_su <- cprd_dataset.dev %>%
  select(all_of(variables_needed)) %>%
  filter(drugclass %in% c("SGLT2", "SU")) %>%
  mutate(drugclass = factor(drugclass)) %>%
  mutate_all(as.numeric) %>% mutate_at(all_of(variables_needed_cat), ~ . - 1)

## targeted maximum likelihood estimates and testing procedure
heterogeneity_sglt2_su_p_values <- unihtee(
  data = heterogeneity_sglt2_su %>% as.matrix() %>% as.data.table(),
  confounders = colnames(heterogeneity_sglt2_su)[-c(1,2)],
  modifiers = colnames(heterogeneity_sglt2_su)[-c(1,2)],
  exposure = "drugclass", outcome = "stopdrug_3m_6mFU", outcome_type = "binary", effect = "relative", estimator = "tmle"
)

#:---- SGLT2 vs TZD
heterogeneity_sglt2_tzd <- cprd_dataset.dev %>%
  select(all_of(variables_needed)) %>%
  filter(drugclass %in% c("SGLT2", "TZD")) %>%
  mutate(drugclass = factor(drugclass)) %>%
  mutate_all(as.numeric) %>% mutate_at(all_of(variables_needed_cat), ~ . - 1)

## targeted maximum likelihood estimates and testing procedure
heterogeneity_sglt2_tzd_p_values <- unihtee(
  data = heterogeneity_sglt2_tzd %>% as.matrix() %>% as.data.table(),
  confounders = colnames(heterogeneity_sglt2_tzd)[-c(1,2)],
  modifiers = colnames(heterogeneity_sglt2_tzd)[-c(1,2)],
  exposure = "drugclass", outcome = "stopdrug_3m_6mFU", outcome_type = "binary", effect = "relative", estimator = "tmle"
)

#:---- TZD vs DPP4
heterogeneity_tzd_dpp4 <- cprd_dataset.dev %>%
  select(all_of(variables_needed)) %>%
  filter(drugclass %in% c("TZD", "DPP4")) %>%
  mutate(drugclass = factor(drugclass)) %>%
  mutate_all(as.numeric) %>% mutate_at(all_of(variables_needed_cat), ~ . - 1)

## targeted maximum likelihood estimates and testing procedure
heterogeneity_tzd_dpp4_p_values <- unihtee(
  data = heterogeneity_tzd_dpp4 %>% as.matrix() %>% as.data.table(),
  confounders = colnames(heterogeneity_tzd_dpp4)[-c(1,2)],
  modifiers = colnames(heterogeneity_tzd_dpp4)[-c(1,2)],
  exposure = "drugclass", outcome = "stopdrug_3m_6mFU", outcome_type = "binary", effect = "relative", estimator = "tmle"
)

#:---- TZD vs SU
heterogeneity_tzd_su <- cprd_dataset.dev %>%
  select(all_of(variables_needed)) %>%
  filter(drugclass %in% c("TZD", "SU")) %>%
  mutate(drugclass = factor(drugclass)) %>%
  mutate_all(as.numeric) %>% mutate_at(all_of(variables_needed_cat), ~ . - 1)

## targeted maximum likelihood estimates and testing procedure
heterogeneity_tzd_su_p_values <- unihtee(
  data = heterogeneity_tzd_su %>% as.matrix() %>% as.data.table(),
  confounders = colnames(heterogeneity_tzd_su)[-c(1,2)],
  modifiers = colnames(heterogeneity_tzd_su)[-c(1,2)],
  exposure = "drugclass", outcome = "stopdrug_3m_6mFU", outcome_type = "binary", effect = "relative", estimator = "tmle"
)

#:---- SU vs DPP4
heterogeneity_su_dpp4 <- cprd_dataset.dev %>%
  select(all_of(variables_needed)) %>%
  filter(drugclass %in% c("SU", "DPP4")) %>%
  mutate(drugclass = factor(drugclass)) %>%
  mutate_all(as.numeric) %>% mutate_at(all_of(variables_needed_cat), ~ . - 1)

## targeted maximum likelihood estimates and testing procedure
heterogeneity_su_dpp4_p_values <- unihtee(
  data = heterogeneity_su_dpp4 %>% as.matrix() %>% as.data.table(),
  confounders = colnames(heterogeneity_su_dpp4)[-c(1,2)],
  modifiers = colnames(heterogeneity_su_dpp4)[-c(1,2)],
  exposure = "drugclass", outcome = "stopdrug_3m_6mFU", outcome_type = "binary", effect = "relative", estimator = "tmle"
)



pdf("results/figures/05.heterogeneity_analysis_p_values.pdf")
rbind(
  heterogeneity_glp1_dpp4_p_values %>% select(modifier, p_value_fdr) %>% mutate(comparison = "GLP1 vs DPP4"),
  heterogeneity_glp1_su_p_values %>% select(modifier, p_value_fdr) %>% mutate(comparison = "GLP1 vs SU"),
  heterogeneity_glp1_tzd_p_values %>% select(modifier, p_value_fdr) %>% mutate(comparison = "GLP1 vs TZD"),
  heterogeneity_sglt2_dpp4_p_values %>% select(modifier, p_value_fdr) %>% mutate(comparison = "SGLT2 vs DPP4"),
  heterogeneity_sglt2_glp1_p_values %>% select(modifier, p_value_fdr) %>% mutate(comparison = "SGLT2 vs GLP1"),
  heterogeneity_sglt2_su_p_values %>% select(modifier, p_value_fdr) %>% mutate(comparison = "SGLT2 vs SU"),
  heterogeneity_sglt2_tzd_p_values %>% select(modifier, p_value_fdr) %>% mutate(comparison = "SGLT2 vs TZD"),
  heterogeneity_tzd_dpp4_p_values %>% select(modifier, p_value_fdr) %>% mutate(comparison = "TZD vs DPP4"),
  heterogeneity_tzd_su_p_values %>% select(modifier, p_value_fdr) %>% mutate(comparison = "TZD vs SU"),
  heterogeneity_su_dpp4_p_values %>% select(modifier, p_value_fdr) %>% mutate(comparison = "SU vs DPP4")
) %>%
  as.data.frame() %>%
  ggplot(aes(x = modifier, y = -log(p_value_fdr), colour = comparison)) +
  geom_hline(yintercept = -log(0.5), colour = "black", linetype = "dashed") +
  geom_point() +
  scale_y_continuous(trans='log10') +
  coord_flip() +
  theme_bw()
dev.off()









###############################################################################
###############################################################################
####################### Per combination comparison ############################
###############################################################################
###############################################################################



#
#:------------ All drugs model
#

# variables to adjust
breakdown <- c(
  # Extra info
  "dstartdate_dm_dur", "dstartdate_age", "drugline", "numdrugs",
  "smoking_cat", "imd2015_10", "gender",
  # Biomarkers
  "prehba1c", "preegfr", "prebmi", "prealt",
  "prehdl"
)


#:--------------
## UBER model variable adjusted

# ATE not adjusted

ATE.var_adj.no_adj_all_drugs <- calc_ATE(cprd_dataset.dev, drugs = c("MFN", "DPP4", "GLP1", "SGLT2", "SU", "TZD"), pred.variable = "no_weight", n_bootstrap = 100)

saveRDS(ATE.var_adj.no_adj_all_drugs, "results/Heterogeneity/05.ATE.var_adj.no_adj_all_drugs.rds")

# ATE overlap matching

ATE.var_adj.overlap_match_all_drugs <- calc_ATE(cprd_dataset.dev, drugs = c("MFN", "DPP4", "GLP1", "SGLT2", "SU", "TZD"), pred.variable = "no_weight", weight.variable = "overlap", matching = TRUE, n_bootstrap = 100)

saveRDS(ATE.var_adj.overlap_match_all_drugs, "results/Heterogeneity/05.ATE.var_adj.overlap_match_all_drugs.rds")

# ATE IPW matching

ATE.var_adj.IPW_match_all_drugs <- calc_ATE(cprd_dataset.dev, drugs = c("MFN", "DPP4", "GLP1", "SGLT2", "SU", "TZD"), pred.variable = "no_weight", weight.variable = "IPW", matching = TRUE, n_bootstrap = 100)

saveRDS(ATE.var_adj.IPW_match_all_drugs, "results/Heterogeneity/05.ATE.var_adj.IPW_match_all_drugs.rds")







pdf("results/figures/05.plot_1.pdf", width = 10, height = 7)

ATE.var_adj.no_adj_all_drugs %>%
  filter(drug_1 == "MFN" | drug_2 == "MFN") %>%
  ggplot() +
  geom_abline(aes(intercept = 0, slope = 1), colour = "red") +
  geom_hline(aes(yintercept = 0), colour = "grey") +
  geom_vline(aes(xintercept = 0), colour = "grey") +
  geom_point(aes(y = obs, x = diff.pred)) +
  geom_errorbar(aes(y = obs, ymin = lci, ymax = uci, x = diff.pred)) +
  scale_x_continuous("Predicted Conditional Average Treatment Discontinuation Difference", labels = scales::percent, breaks = seq(-1, 1, 0.02)) +
  scale_y_continuous("Decile Average Treatment Discontinuation Difference", labels = scales::percent, breaks = seq(-1, 1, 0.02)) +
  ggtitle("Comparison of MFN against other therapies") +
  facet_wrap(~Type, scales = "free", nrow = 2) +
  theme_bw()

dev.off()



pdf("results/figures/05.plot_2.pdf", width = 12, height = 9)

ATE.var_adj.no_adj_all_drugs %>%
  filter(drug_1 != "MFN" & drug_2 != "MFN") %>%
  ggplot() +
  geom_abline(aes(intercept = 0, slope = 1), colour = "red") +
  geom_hline(aes(yintercept = 0), colour = "grey") +
  geom_vline(aes(xintercept = 0), colour = "grey") +
  geom_point(aes(y = obs, x = diff.pred)) +
  geom_errorbar(aes(y = obs, ymin = lci, ymax = uci, x = diff.pred)) +
  scale_x_continuous("Predicted Conditional Average Treatment Discontinuation Difference", labels = scales::percent, breaks = seq(-1, 1, 0.02)) +
  scale_y_continuous("Decile Average Treatment Discontinuation Difference (no adjustment)", labels = scales::percent, breaks = seq(-1, 1, 0.02)) +
  ggtitle("Therapy discontinuation heterogeneity (no adjustment)") +
  facet_wrap(~Type, scales = "free") +
  theme_bw()

ATE.var_adj.overlap_match_all_drugs %>%
  filter(drug_1 != "MFN" & drug_2 != "MFN") %>%
  ggplot() +
  geom_abline(aes(intercept = 0, slope = 1), colour = "red") +
  geom_hline(aes(yintercept = 0), colour = "grey") +
  geom_vline(aes(xintercept = 0), colour = "grey") +
  geom_point(aes(y = obs, x = diff.pred)) +
  geom_errorbar(aes(y = obs, ymin = lci, ymax = uci, x = diff.pred)) +
  scale_x_continuous("Predicted Conditional Average Treatment Discontinuation Difference", labels = scales::percent, breaks = seq(-1, 1, 0.02)) +
  scale_y_continuous("Decile Average Treatment Discontinuation Difference", labels = scales::percent, breaks = seq(-1, 1, 0.02)) +
  ggtitle("Therapy discontinuation heterogeneity (Overlap weight matching)") +
  facet_wrap(~Type, scales = "free") +
  theme_bw()

ATE.var_adj.IPW_match_all_drugs %>%
  filter(drug_1 != "MFN" & drug_2 != "MFN") %>%
  ggplot() +
  geom_abline(aes(intercept = 0, slope = 1), colour = "red") +
  geom_hline(aes(yintercept = 0), colour = "grey") +
  geom_vline(aes(xintercept = 0), colour = "grey") +
  geom_point(aes(y = obs, x = diff.pred)) +
  geom_errorbar(aes(y = obs, ymin = lci, ymax = uci, x = diff.pred)) +
  scale_x_continuous("Predicted Conditional Average Treatment Discontinuation Difference", labels = scales::percent, breaks = seq(-1, 1, 0.02)) +
  scale_y_continuous("Decile Average Treatment Discontinuation Difference", labels = scales::percent, breaks = seq(-1, 1, 0.02)) +
  ggtitle("Therapy discontinuation heterogeneity (Overlap weight matching)") +
  facet_wrap(~Type, scales = "free") +
  theme_bw()

dev.off()












#
#:------------ 2nd line drugs model
#

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
  dataset = "3m.disc.dataset.dev",
  follow_up = "6-months",
  full_prescribing_history = TRUE
) %>%
  drop_na(-stopdrug_6m_6mFU, -stopdrug_12m_6mFU)


# load dataset
cprd_dataset.val <- set_up_data(
  raw_data = "20240308_t2d_1stinstance",
  diagnosis = FALSE,
  therapies = c("DPP4", "GLP1", "MFN", "SGLT2", "SU", "TZD"),
  dataset = "3m.disc.dataset.val",
  follow_up = "6-months",
  full_prescribing_history = TRUE
) %>%
  drop_na(-stopdrug_6m_6mFU, -stopdrug_12m_6mFU)



# load propensity scores
ps.only_dataset <- readRDS("results/PS_model/ps.dataset_lm_all.rds")


# join propensity scores into dataset
cprd_dataset.dev <- cprd_dataset.dev %>%
  left_join(
    ps.only_dataset %>%
      select(patid, dstartdate, prop.score.overlap, prop.score.IPW),
    by = c("patid", "dstartdate")
  )
cprd_dataset.val <- cprd_dataset.val %>%
  left_join(
    ps.only_dataset %>%
      select(patid, dstartdate, prop.score.overlap, prop.score.IPW),
    by = c("patid", "dstartdate")
  )


probabilities.only_dataset <- readRDS("results/Models/Predictions/04.model_predictions_3m_2nd_line_drugs.rds")


# join predictions into dataset
cprd_dataset.dev <- cprd_dataset.dev %>%
  left_join(
    probabilities.only_dataset,
    by = c("patid", "dstartdate")
  )
cprd_dataset.val <- cprd_dataset.val %>%
  left_join(
    probabilities.only_dataset,
    by = c("patid", "dstartdate")
  )



# variables to adjust
breakdown <- c(
  # Extra info
  "dstartdate_dm_dur", "dstartdate_age", "drugline", "numdrugs",
  "smoking_cat", "imd2015_10", "gender",
  # Biomarkers
  "prehba1c", "preegfr", "prebmi", "prealt",
  "prehdl"
)


#:--------------
## UBER model variable adjusted

# ATE not adjusted

ATE.var_adj.no_adj_2nd_line <- calc_ATE(cprd_dataset.dev, drugs = c("DPP4", "GLP1", "SGLT2", "SU", "TZD"), pred.variable = "no_weight", n_bootstrap = 100)

saveRDS(ATE.var_adj.no_adj_2nd_line, "results/Heterogeneity/05.ATE.var_adj.no_adj_2nd_line.rds")

# ATE overlap matching

ATE.var_adj.overlap_match_2nd_line <- calc_ATE(cprd_dataset.dev, drugs = c("DPP4", "GLP1", "SGLT2", "SU", "TZD"), pred.variable = "no_weight", weight.variable = "overlap", matching = TRUE, n_bootstrap = 100)

saveRDS(ATE.var_adj.overlap_match_2nd_line, "results/Heterogeneity/05.ATE.var_adj.overlap_match_2nd_line.rds")

# ATE IPW matching

ATE.var_adj.IPW_match_2nd_line <- calc_ATE(cprd_dataset.dev, drugs = c("DPP4", "GLP1", "SGLT2", "SU", "TZD"), pred.variable = "no_weight", weight.variable = "IPW", matching = TRUE, n_bootstrap = 100)

saveRDS(ATE.var_adj.IPW_match_2nd_line, "results/Heterogeneity/05.ATE.var_adj.IPW_match_2nd_line.rds")





pdf("results/figures/05.plot_3.pdf", width = 12, height = 9)

ATE.var_adj.no_adj_2nd_line %>%
  ggplot() +
  geom_abline(aes(intercept = 0, slope = 1), colour = "red") +
  geom_hline(aes(yintercept = 0), colour = "grey") +
  geom_vline(aes(xintercept = 0), colour = "grey") +
  geom_point(aes(y = obs, x = diff.pred)) +
  geom_errorbar(aes(y = obs, ymin = lci, ymax = uci, x = diff.pred)) +
  scale_x_continuous("Predicted Conditional Average Treatment Discontinuation Difference", labels = scales::percent, breaks = seq(-1, 1, 0.02)) +
  scale_y_continuous("Decile Average Treatment Discontinuation Difference (no adjustment)", labels = scales::percent, breaks = seq(-1, 1, 0.02)) +
  ggtitle("2nd line Therapy discontinuation heterogeneity (no adjustment)") +
  facet_wrap(~Type, scales = "free") +
  theme_bw()

ATE.var_adj.no_adj_2nd_line %>%
  ggplot() +
  geom_abline(aes(intercept = 0, slope = 1), colour = "red") +
  geom_hline(aes(yintercept = 0), colour = "grey") +
  geom_vline(aes(xintercept = 0), colour = "grey") +
  geom_point(aes(y = obs, x = diff.pred)) +
  geom_errorbar(aes(y = obs, ymin = lci, ymax = uci, x = diff.pred)) +
  scale_x_continuous("Predicted Conditional Average Treatment Discontinuation Difference", labels = scales::percent, breaks = seq(-1, 1, 0.02)) +
  scale_y_continuous("Decile Average Treatment Discontinuation Difference (no adjustment)", labels = scales::percent, breaks = seq(-1, 1, 0.02)) +
  ggtitle("2nd line Therapy discontinuation heterogeneity (no adjustment)") +
  facet_wrap(~Type, scales = "fixed") +
  theme_bw()

ATE.var_adj.overlap_match_2nd_line %>%
  ggplot() +
  geom_abline(aes(intercept = 0, slope = 1), colour = "red") +
  geom_hline(aes(yintercept = 0), colour = "grey") +
  geom_vline(aes(xintercept = 0), colour = "grey") +
  geom_point(aes(y = obs, x = diff.pred)) +
  geom_errorbar(aes(y = obs, ymin = lci, ymax = uci, x = diff.pred)) +
  scale_x_continuous("Predicted Conditional Average Treatment Discontinuation Difference", labels = scales::percent, breaks = seq(-1, 1, 0.02)) +
  scale_y_continuous("Decile Average Treatment Discontinuation Difference", labels = scales::percent, breaks = seq(-1, 1, 0.02)) +
  ggtitle("2nd line Therapy discontinuation heterogeneity (Overlap weight matching)") +
  facet_wrap(~Type, scales = "free") +
  theme_bw()

ATE.var_adj.IPW_match_2nd_line %>%
  ggplot() +
  geom_abline(aes(intercept = 0, slope = 1), colour = "red") +
  geom_hline(aes(yintercept = 0), colour = "grey") +
  geom_vline(aes(xintercept = 0), colour = "grey") +
  geom_point(aes(y = obs, x = diff.pred)) +
  geom_errorbar(aes(y = obs, ymin = lci, ymax = uci, x = diff.pred)) +
  scale_x_continuous("Predicted Conditional Average Treatment Discontinuation Difference", labels = scales::percent, breaks = seq(-1, 1, 0.02)) +
  scale_y_continuous("Decile Average Treatment Discontinuation Difference", labels = scales::percent, breaks = seq(-1, 1, 0.02)) +
  ggtitle("2nd line Therapy discontinuation heterogeneity (Overlap weight matching)") +
  facet_wrap(~Type, scales = "free") +
  theme_bw()

dev.off()





#:-----------------------------------------------------------------
# Combining pdfs

qpdf::pdf_combine(input = c("results/figures/05.plot_1.pdf", "results/figures/05.plot_2.pdf", "results/figures/05.plot_3.pdf"),
                  output = "results/figures/05.heterogeneity_differences_analysis.pdf")


file.remove(c("results/figures/05.plot_1.pdf", "results/figures/05.plot_2.pdf", "results/figures/05.plot_3.pdf"))
