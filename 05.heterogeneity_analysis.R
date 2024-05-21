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
  drop_na(-stopdrug_6m_6mFU, -stopdrug_12m_6mFU)


# load dataset
cprd_dataset.val <- set_up_data(
  raw_data = "20240308_t2d_1stinstance",
  diagnosis = FALSE,
  therapies = c("DPP4", "GLP1", "MFN", "SGLT2", "SU", "TZD"),
  dataset = "3m.disc.dataset.val"
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


probabilities.only_dataset <- readRDS("results/Models/Predictions/model_predictions_3m_all_drugs.rds")


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
  "pretotalcholesterol"
)


#:--------------
## UBER model variable adjusted

# ATE not adjusted

ATE.var_adj.no_adj_all_drugs <- calc_ATE(cprd_dataset.dev, drugs = c("MFN", "DPP4", "GLP1", "SGLT2", "SU", "TZD"), pred.variable = "no_weight", n_bootstrap = 100)

saveRDS(ATE.var_adj.no_adj_all_drugs, "results/Heterogeneity/ATE.var_adj.no_adj_all_drugs.rds")

# ATE overlap matching

ATE.var_adj.overlap_match_all_drugs <- calc_ATE(cprd_dataset.dev, drugs = c("MFN", "DPP4", "GLP1", "SGLT2", "SU", "TZD"), pred.variable = "no_weight", weight.variable = "overlap", matching = TRUE, n_bootstrap = 100)

saveRDS(ATE.var_adj.overlap_match_all_drugs, "results/Heterogeneity/ATE.var_adj.overlap_match_all_drugs.rds")

# ATE IPW matching

ATE.var_adj.IPW_match_all_drugs <- calc_ATE(cprd_dataset.dev, drugs = c("MFN", "DPP4", "GLP1", "SGLT2", "SU", "TZD"), pred.variable = "no_weight", weight.variable = "IPW", matching = TRUE, n_bootstrap = 100)

saveRDS(ATE.var_adj.IPW_match_all_drugs, "results/Heterogeneity/ATE.var_adj.IPW_match_all_drugs.rds")







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
  dataset = "3m.disc.dataset.dev"
) %>%
  drop_na(-stopdrug_6m_6mFU, -stopdrug_12m_6mFU)


# load dataset
cprd_dataset.val <- set_up_data(
  raw_data = "20240308_t2d_1stinstance",
  diagnosis = FALSE,
  therapies = c("DPP4", "GLP1", "MFN", "SGLT2", "SU", "TZD"),
  dataset = "3m.disc.dataset.val"
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


probabilities.only_dataset <- readRDS("results/Models/Predictions/model_predictions_3m_2nd_line_drugs.rds")


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
  "pretotalcholesterol"
)


#:--------------
## UBER model variable adjusted

# ATE not adjusted

ATE.var_adj.no_adj_2nd_line <- calc_ATE(cprd_dataset.dev, drugs = c("DPP4", "GLP1", "SGLT2", "SU", "TZD"), pred.variable = "no_weight", n_bootstrap = 100)

saveRDS(ATE.var_adj.no_adj_2nd_line, "results/Heterogeneity/ATE.var_adj.no_adj_2nd_line.rds")

# ATE overlap matching

ATE.var_adj.overlap_match_2nd_line <- calc_ATE(cprd_dataset.dev, drugs = c("DPP4", "GLP1", "SGLT2", "SU", "TZD"), pred.variable = "no_weight", weight.variable = "overlap", matching = TRUE, n_bootstrap = 100)

saveRDS(ATE.var_adj.overlap_match_2nd_line, "results/Heterogeneity/ATE.var_adj.overlap_match_2nd_line.rds")

# ATE IPW matching

ATE.var_adj.IPW_match_2nd_line <- calc_ATE(cprd_dataset.dev, drugs = c("DPP4", "GLP1", "SGLT2", "SU", "TZD"), pred.variable = "no_weight", weight.variable = "IPW", matching = TRUE, n_bootstrap = 100)

saveRDS(ATE.var_adj.IPW_match_2nd_line, "results/Heterogeneity/ATE.var_adj.IPW_match_2nd_line.rds")





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
                  output = "results/figures/heterogeneity_differences_analysis.pdf")


file.remove(c("results/figures/05.plot_1.pdf", "results/figures/05.plot_2.pdf", "results/figures/05.plot_3.pdf"))
