####################
## Description: 
##  - In this file we:
##    - Fit a linear model all the variables into one model.
#################### 


# load working directory
setwd("Samples/T2D_Discontinuation")

# load functions
source("code/00.set_up_data_and_functions.R")

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
  therapies = c("MFN", "GLP1", "DPP4", "SGLT2", "TZD", "SU"),
  dataset = "3m.disc.dataset.dev",
  follow_up = "6-months",
  full_prescribing_history = TRUE
) %>%
  drop_na(-stopdrug_6m_6mFU, -stopdrug_12m_6mFU)

# load dataset
cprd_dataset.val <- set_up_data(
  raw_data = "20240308_t2d_1stinstance",
  diagnosis = FALSE,
  therapies = c("MFN", "GLP1", "DPP4", "SGLT2", "TZD", "SU"),
  dataset = "3m.disc.dataset.val",
  follow_up = "6-months",
  full_prescribing_history = TRUE
) %>%
  drop_na(-stopdrug_6m_6mFU, -stopdrug_12m_6mFU)


###########################################################
###########################################################

transformation_matrix <- NULL

transformation_matrix[["dstartdate_dm_dur"]] <- data.frame(cbind(mean = mean(cprd_dataset.dev$dstartdate_dm_dur, na.rm = TRUE), sd = sd(cprd_dataset.dev$dstartdate_dm_dur, na.rm = TRUE)))
transformation_matrix[["dstartdate_age"]] <- data.frame(cbind(mean = mean(cprd_dataset.dev$dstartdate_age, na.rm = TRUE), sd = sd(cprd_dataset.dev$dstartdate_age, na.rm = TRUE)))
transformation_matrix[["prehba1c"]] <- data.frame(cbind(mean = mean(cprd_dataset.dev$prehba1c, na.rm = TRUE), sd = sd(cprd_dataset.dev$prehba1c, na.rm = TRUE)))
transformation_matrix[["preegfr"]] <- data.frame(cbind(mean = mean(cprd_dataset.dev$preegfr, na.rm = TRUE), sd = sd(cprd_dataset.dev$preegfr, na.rm = TRUE)))
transformation_matrix[["prebmi"]] <- data.frame(cbind(mean = mean(cprd_dataset.dev$prebmi, na.rm = TRUE), sd = sd(cprd_dataset.dev$prebmi, na.rm = TRUE)))
transformation_matrix[["prealt"]] <- data.frame(cbind(mean = mean(cprd_dataset.dev$prealt, na.rm = TRUE), sd = sd(cprd_dataset.dev$prealt, na.rm = TRUE)))
transformation_matrix[["prehdl"]] <- data.frame(cbind(mean = mean(cprd_dataset.dev$prehdl, na.rm = TRUE), sd = sd(cprd_dataset.dev$prehdl, na.rm = TRUE)))


cprd_dataset.dev <- cprd_dataset.dev %>%
  mutate(
    imd2015_10 = as.numeric(imd2015_10),
    drugline = as.numeric(drugline)
  ) %>%
  mutate(
    dstartdate_dm_dur = (cprd_dataset.dev$dstartdate_dm_dur - transformation_matrix[["dstartdate_dm_dur"]][[1]])/transformation_matrix[["dstartdate_dm_dur"]][[2]],
    dstartdate_age = (cprd_dataset.dev$dstartdate_age - transformation_matrix[["dstartdate_age"]][[1]])/transformation_matrix[["dstartdate_age"]][[2]],
    prehba1c = (cprd_dataset.dev$prehba1c - transformation_matrix[["prehba1c"]][[1]])/transformation_matrix[["prehba1c"]][[2]],
    preegfr = (cprd_dataset.dev$preegfr - transformation_matrix[["preegfr"]][[1]])/transformation_matrix[["preegfr"]][[2]],
    prebmi = (cprd_dataset.dev$prebmi - transformation_matrix[["prebmi"]][[1]])/transformation_matrix[["prebmi"]][[2]],
    prealt = (cprd_dataset.dev$prealt - transformation_matrix[["prealt"]][[1]])/transformation_matrix[["prealt"]][[2]],
    prehdl = (cprd_dataset.dev$prehdl - transformation_matrix[["prehdl"]][[1]])/transformation_matrix[["prehdl"]][[2]]
  ) %>%
  mutate(ethnicity_5cat = factor(ethnicity_5cat, levels = c(0, 1, 2, 3, 4), labels = c("White", "Other", "Other" , "Other" , "Other"))) %>%
  mutate(smoking_cat = factor(smoking_cat, levels = c("Active smoker", "Ex-smoker", "Non-smoker"), labels = c("Active smoker", "Other", "Other"))) %>%
  mutate(alcohol_cat = factor(alcohol_cat, levels = c("Excess", "Harmful", "None", "Within limits"), labels = c("Excess/Harmful", "Excess/Harmful", "Within limits/None", "Within limits/None"))) %>%
  mutate(preckdstage = factor(preckdstage, levels = c("stage_0", "stage_1", "stage_2", "stage_3a", "stage_3b", "stage_4", "stage_5"), labels = c("stage 0/1/2", "stage 0/1/2", "stage 0/1/2", "stage 3/4/5", "stage 3/4/5", "stage 3/4/5", "stage 3/4/5")))


cprd_dataset.val <- cprd_dataset.val %>%
  mutate(
    imd2015_10 = as.numeric(imd2015_10),
    drugline = as.numeric(drugline)
  ) %>%
  mutate(
    dstartdate_dm_dur = (cprd_dataset.val$dstartdate_dm_dur - transformation_matrix[["dstartdate_dm_dur"]][[1]])/transformation_matrix[["dstartdate_dm_dur"]][[2]],
    dstartdate_age = (cprd_dataset.val$dstartdate_age - transformation_matrix[["dstartdate_age"]][[1]])/transformation_matrix[["dstartdate_age"]][[2]],
    prehba1c = (cprd_dataset.val$prehba1c - transformation_matrix[["prehba1c"]][[1]])/transformation_matrix[["prehba1c"]][[2]],
    preegfr = (cprd_dataset.val$preegfr - transformation_matrix[["preegfr"]][[1]])/transformation_matrix[["preegfr"]][[2]],
    prebmi = (cprd_dataset.val$prebmi - transformation_matrix[["prebmi"]][[1]])/transformation_matrix[["prebmi"]][[2]],
    prealt = (cprd_dataset.val$prealt - transformation_matrix[["prealt"]][[1]])/transformation_matrix[["prealt"]][[2]],
    prehdl = (cprd_dataset.val$prehdl - transformation_matrix[["prehdl"]][[1]])/transformation_matrix[["prehdl"]][[2]]
  ) %>%
  mutate(ethnicity_5cat = factor(ethnicity_5cat, levels = c(0, 1, 2, 3, 4), labels = c("White", "Other", "Other" , "Other" , "Other"))) %>%
  mutate(smoking_cat = factor(smoking_cat, levels = c("Active smoker", "Ex-smoker", "Non-smoker"), labels = c("Active smoker", "Other", "Other"))) %>%
  mutate(alcohol_cat = factor(alcohol_cat, levels = c("Excess", "Harmful", "None", "Within limits"), labels = c("Excess/Harmful", "Excess/Harmful", "Within limits/None", "Within limits/None"))) %>%
  mutate(preckdstage = factor(preckdstage, levels = c("stage_0", "stage_1", "stage_2", "stage_3a", "stage_3b", "stage_4", "stage_5"), labels = c("stage 0/1/2", "stage 0/1/2", "stage 0/1/2", "stage 3/4/5", "stage 3/4/5", "stage 3/4/5", "stage 3/4/5")))


roc_dataset <- NULL

formula_high_vars <- "~ dstartdate_age + gender + imd2015_10 + prebmi + dstartdate_dm_dur + prehba1c + drugline + predrug_frailty_proxy + ethnicity_5cat"

outcome <- "stopdrug_3m_6mFU"

lm.3m.mfn.high_vars <- glm(formula = as.formula(paste0(outcome, formula_high_vars)), data = cprd_dataset.dev %>% filter(drugclass == "MFN"), family = binomial())

roc_dataset <- rbind(
  roc_dataset,
  pROC::ci(pROC::roc(
  cprd_dataset.dev %>% filter(drugclass == "MFN") %>% select(all_of(outcome)) %>% unlist(),
  predict(lm.3m.mfn.high_vars, type = "response", newdata = cprd_dataset.dev %>% filter(drugclass == "MFN"))
    ), of = "thresholds") %>%
    as.data.frame() %>%
    mutate(model = "MFN", dataset = "dev", variables = "high"),
  pROC::ci(pROC::roc(
    cprd_dataset.val %>% filter(drugclass == "MFN") %>% select(all_of(outcome)) %>% unlist(),
    predict(lm.3m.mfn.high_vars, type = "response", newdata = cprd_dataset.val %>% filter(drugclass == "MFN"))
  ), of = "thresholds") %>%
    as.data.frame() %>%
    mutate(model = "MFN", dataset = "val", variables = "high")
)

lm.3m.otherdrugs.high_vars <- glm(formula = as.formula(paste0(outcome, formula_high_vars, "+stopdrug_3m_3mFU_MFN_hist")), data = cprd_dataset.dev %>% filter(drugclass != "MFN"), family = binomial())

roc_dataset <- rbind(
  roc_dataset,
  pROC::ci(pROC::roc(
    cprd_dataset.dev %>% filter(drugclass != "MFN") %>% select(all_of(outcome)) %>% unlist(),
    predict(lm.3m.mfn.high_vars, type = "response", newdata = cprd_dataset.dev %>% filter(drugclass != "MFN"))
  ), of = "thresholds") %>%
    as.data.frame() %>%
    mutate(model = "Other", dataset = "dev", variables = "high"),
  pROC::ci(pROC::roc(
    cprd_dataset.val %>% filter(drugclass != "MFN") %>% select(all_of(outcome)) %>% unlist(),
    predict(lm.3m.mfn.high_vars, type = "response", newdata = cprd_dataset.val %>% filter(drugclass != "MFN"))
  ), of = "thresholds") %>%
    as.data.frame() %>%
    mutate(model = "Other", dataset = "val", variables = "high")
)


formula_high_vars_extra <- "~ dstartdate_age + gender + imd2015_10 + prebmi + dstartdate_dm_dur + prehba1c + drugline + predrug_frailty_proxy + ethnicity_5cat + numdrugs + predrug_bloodmed + smoking_cat + predrug_statins + preegfr + prehdl"

lm.3m.mfn.high_vars_extra <- glm(formula = as.formula(paste0(outcome, formula_high_vars_extra)), data = cprd_dataset.dev %>% filter(drugclass == "MFN"), family = binomial())

roc_dataset <- rbind(
  roc_dataset,
  pROC::ci(pROC::roc(
    cprd_dataset.dev %>% filter(drugclass == "MFN") %>% select(all_of(outcome)) %>% unlist(),
    predict(lm.3m.mfn.high_vars_extra, type = "response", newdata = cprd_dataset.dev %>% filter(drugclass == "MFN"))
  ), of = "thresholds") %>%
    as.data.frame() %>%
    mutate(model = "MFN", dataset = "dev", variables = "high_extra"),
  pROC::ci(pROC::roc(
    cprd_dataset.val %>% filter(drugclass == "MFN") %>% select(all_of(outcome)) %>% unlist(),
    predict(lm.3m.mfn.high_vars_extra, type = "response", newdata = cprd_dataset.val %>% filter(drugclass == "MFN"))
  ), of = "thresholds") %>%
    as.data.frame() %>%
    mutate(model = "MFN", dataset = "val", variables = "high_extra")
)


lm.3m.otherdrugs.high_vars_extra <- glm(formula = as.formula(paste0(outcome, formula_high_vars_extra, "+stopdrug_3m_3mFU_MFN_hist")), data = cprd_dataset.dev %>% filter(drugclass != "MFN"), family = binomial())


roc_dataset <- rbind(
  roc_dataset,
  pROC::ci(pROC::roc(
    cprd_dataset.dev %>% filter(drugclass != "MFN") %>% select(all_of(outcome)) %>% unlist(),
    predict(lm.3m.otherdrugs.high_vars_extra, type = "response", newdata = cprd_dataset.dev %>% filter(drugclass != "MFN"))
  ), of = "thresholds") %>%
    as.data.frame() %>%
    mutate(model = "Other", dataset = "dev", variables = "high_extra"),
  pROC::ci(pROC::roc(
    cprd_dataset.val %>% filter(drugclass != "MFN") %>% select(all_of(outcome)) %>% unlist(),
    predict(lm.3m.otherdrugs.high_vars_extra, type = "response", newdata = cprd_dataset.val %>% filter(drugclass != "MFN"))
  ), of = "thresholds") %>%
    as.data.frame() %>%
    mutate(model = "Other", dataset = "val", variables = "high_extra")
)


saveRDS(roc_dataset, "results/Models/Predictions/04.roc_plotting_multivariate_disc.rds")



pdf("results/figures/04.roc_multivariate_disc.pdf", width = 8, height = 6)

roc_dataset %>%
  mutate(
    model = factor(model, levels = c("MFN", "Other")),
    dataset = factor(dataset, levels = c("dev", "val"), labels = c("Development", "Validation")),
    variables = factor(variables, levels = c("high", "high_extra"), labels = c("Highlights", "Highlights + extra"))
  ) %>%
  ggplot() +
  geom_abline(aes(intercept = 0, slope = 1), colour = "black", linetype = "dashed") +
  geom_path(aes(x = 1 - specificity.2.5., y = sensitivity.2.5., colour = model), alpha = 0.4) +
  geom_path(aes(x = 1 - specificity.50., y = sensitivity.50., colour = model)) +
  geom_path(aes(x = 1 - specificity.97.5., y = sensitivity.97.5., colour = model), alpha = 0.4) +
  scale_colour_manual(values = c("Other" = "black", "Pooled" = "black", "SGLT2" = "#E69F00", "GLP1" = "#56B4E9", "SU" = "#CC79A7", "DPP4" = "#0072B2", "TZD" = "#D55E00", "MFN" = "grey"), 
                    breaks = c("MFN", "Other"), labels = c("Metformin", "Other therapies"),
                    name = "", guide = guide_legend(reverse = TRUE)) +
  xlab("1 - Specificity") +
  ylab("Sensitivity") +
  xlim(0, 1) +
  ylim(0, 1) +
  facet_grid(dataset ~ variables) +
  theme_bw() +
  theme(
    legend.position = "bottom"
  )


dev.off()

###########################################################
###########################################################


# load dataset
cprd_dataset.dev <- set_up_data(
  raw_data = "20240308_t2d_1stinstance",
  diagnosis = FALSE,
  therapies = c("MFN", "GLP1", "DPP4", "SGLT2", "TZD", "SU"),
  dataset = "3m.disc.dataset.dev",
  follow_up = "6-months",
  full_prescribing_history = TRUE
) %>%
  drop_na(-stopdrug_6m_6mFU, -stopdrug_12m_6mFU)

# load dataset
cprd_dataset.val <- set_up_data(
  raw_data = "20240308_t2d_1stinstance",
  diagnosis = FALSE,
  therapies = c("MFN", "GLP1", "DPP4", "SGLT2", "TZD", "SU"),
  dataset = "3m.disc.dataset.val",
  follow_up = "6-months",
  full_prescribing_history = TRUE
) %>%
  drop_na(-stopdrug_6m_6mFU, -stopdrug_12m_6mFU)



#:-------------------------------
## Fitting models

calc_roc_dataset <- function(data, model, outcome_var, model_chr, outcome_chr, dataset_chr, drug_chr) {
  
  roc_dataset <- rbind(
    data.frame(
      model = model_chr,
      outcome = outcome_chr,
      dataset = dataset_chr,
      drug = drug_chr,
      N = nrow(data),
      events = data %>% rename("outcome_var" = outcome_var) %>% filter(outcome_var == 1) %>% nrow(),
      roc = as.numeric(roc(data %>% select(all_of(outcome_var)) %>% unlist(), predict(lm.3m.pooled.high_vars, type = "response", newdata = data), ci = TRUE)$ci[2]),
      roc_lci = as.numeric(roc(data %>% select(all_of(outcome_var)) %>% unlist(), predict(lm.3m.pooled.high_vars, type = "response", newdata = data), ci = TRUE)$ci[1]),
      roc_uci = as.numeric(roc(data %>% select(all_of(outcome_var)) %>% unlist(), predict(lm.3m.pooled.high_vars, type = "response", newdata = data), ci = TRUE)$ci[3])
    )
  )
  
}

roc_dataset <- NULL

### Prediction models based on highlighted variables

formula_high_vars <- "~ dstartdate_age + gender + imd2015_10 + prebmi + dstartdate_dm_dur + prehba1c + drugline + predrug_frailty_proxy + ethnicity_5cat"

outcome <- "stopdrug_3m_6mFU"

#### Pooled

lm.3m.pooled.high_vars <- glm(formula = as.formula(paste0(outcome, formula_high_vars)), data = cprd_dataset.dev, family = binomial())

roc_dataset <- rbind(
  roc_dataset,
  calc_roc_dataset(cprd_dataset.dev, lm.3m.pooled.high_vars, outcome, model_chr = "Highlighted vars", outcome_chr = "3-months", dataset_chr = "Development", drug_chr = "Pooled"),
  calc_roc_dataset(cprd_dataset.val, lm.3m.pooled.high_vars, outcome, model_chr = "Highlighted vars", outcome_chr = "3-months", dataset_chr = "Validation", drug_chr = "Pooled")
)

#### By drugclass

# MFN
lm.3m.MFN.high_vars <- glm(formula = as.formula(paste0(outcome, formula_high_vars)), data = cprd_dataset.dev %>% filter(drugclass == "MFN"), family = binomial())

roc_dataset <- rbind(
  roc_dataset,
  calc_roc_dataset(cprd_dataset.dev %>% filter(drugclass == "MFN"), lm.3m.MFN.high_vars, outcome, model_chr = "Highlighted vars", outcome_chr = "3-months", dataset_chr = "Development", drug_chr = "MFN"),
  calc_roc_dataset(cprd_dataset.val %>% filter(drugclass == "MFN"), lm.3m.MFN.high_vars, outcome, model_chr = "Highlighted vars", outcome_chr = "3-months", dataset_chr = "Validation", drug_chr = "MFN")
)

# GLP1
lm.3m.GLP1.high_vars <- glm(formula = as.formula(paste0(outcome, formula_high_vars)), data = cprd_dataset.dev %>% filter(drugclass == "GLP1"), family = binomial())

roc_dataset <- rbind(
  roc_dataset,
  calc_roc_dataset(cprd_dataset.dev %>% filter(drugclass == "GLP1"), lm.3m.GLP1.high_vars, outcome, model_chr = "Highlighted vars", outcome_chr = "3-months", dataset_chr = "Development", drug_chr = "GLP1"),
  calc_roc_dataset(cprd_dataset.val %>% filter(drugclass == "GLP1"), lm.3m.GLP1.high_vars, outcome, model_chr = "Highlighted vars", outcome_chr = "3-months", dataset_chr = "Validation", drug_chr = "GLP1")
)

# DPP4
lm.3m.DPP4.high_vars <- glm(formula = as.formula(paste0(outcome, formula_high_vars)), data = cprd_dataset.dev %>% filter(drugclass == "DPP4"), family = binomial())

roc_dataset <- rbind(
  roc_dataset,
  calc_roc_dataset(cprd_dataset.dev %>% filter(drugclass == "DPP4"), lm.3m.DPP4.high_vars, outcome, model_chr = "Highlighted vars", outcome_chr = "3-months", dataset_chr = "Development", drug_chr = "DPP4"),
  calc_roc_dataset(cprd_dataset.val %>% filter(drugclass == "DPP4"), lm.3m.DPP4.high_vars, outcome, model_chr = "Highlighted vars", outcome_chr = "3-months", dataset_chr = "Validation", drug_chr = "DPP4")
)

# SGLT2
lm.3m.SGLT2.high_vars <- glm(formula = as.formula(paste0(outcome, formula_high_vars)), data = cprd_dataset.dev %>% filter(drugclass == "SGLT2"), family = binomial())

roc_dataset <- rbind(
  roc_dataset,
  calc_roc_dataset(cprd_dataset.dev %>% filter(drugclass == "SGLT2"), lm.3m.SGLT2.high_vars, outcome, model_chr = "Highlighted vars", outcome_chr = "3-months", dataset_chr = "Development", drug_chr = "SGLT2"),
  calc_roc_dataset(cprd_dataset.val %>% filter(drugclass == "SGLT2"), lm.3m.SGLT2.high_vars, outcome, model_chr = "Highlighted vars", outcome_chr = "3-months", dataset_chr = "Validation", drug_chr = "SGLT2")
)

# TZD
lm.3m.TZD.high_vars <- glm(formula = as.formula(paste0(outcome, formula_high_vars)), data = cprd_dataset.dev %>% filter(drugclass == "TZD"), family = binomial())

roc_dataset <- rbind(
  roc_dataset,
  calc_roc_dataset(cprd_dataset.dev %>% filter(drugclass == "TZD"), lm.3m.TZD.high_vars, outcome, model_chr = "Highlighted vars", outcome_chr = "3-months", dataset_chr = "Development", drug_chr = "TZD"),
  calc_roc_dataset(cprd_dataset.val %>% filter(drugclass == "TZD"), lm.3m.TZD.high_vars, outcome, model_chr = "Highlighted vars", outcome_chr = "3-months", dataset_chr = "Validation", drug_chr = "TZD")
)

# SU
lm.3m.SU.high_vars <- glm(formula = as.formula(paste0(outcome, formula_high_vars)), data = cprd_dataset.dev %>% filter(drugclass == "SU"), family = binomial())

roc_dataset <- rbind(
  roc_dataset,
  calc_roc_dataset(cprd_dataset.dev %>% filter(drugclass == "SU"), lm.3m.SU.high_vars, outcome, model_chr = "Highlighted vars", outcome_chr = "3-months", dataset_chr = "Development", drug_chr = "SU"),
  calc_roc_dataset(cprd_dataset.val %>% filter(drugclass == "SU"), lm.3m.SU.high_vars, outcome, model_chr = "Highlighted vars", outcome_chr = "3-months", dataset_chr = "Validation", drug_chr = "SU")
)

### Prediction models based on highlighted variable + additional

formula_high_vars_extra <- "~ dstartdate_age + gender + imd2015_10 + prebmi + dstartdate_dm_dur + prehba1c + drugline + predrug_frailty_proxy + ethnicity_5cat + numdrugs + predrug_bloodmed + smoking_cat + predrug_statins + preegfr + prehdl"

#### Pooled

lm.3m.pooled.high_vars_extra <- glm(formula = as.formula(paste0(outcome, formula_high_vars_extra, "+ drugclass:stopdrug_3m_3mFU_MFN_hist")), data = cprd_dataset.dev, family = binomial())

roc_dataset <- rbind(
  roc_dataset,
  calc_roc_dataset(cprd_dataset.dev, lm.3m.pooled.high_vars_extra, outcome, model_chr = "Highlighted vars + extra", outcome_chr = "3-months", dataset_chr = "Development", drug_chr = "Pooled"),
  calc_roc_dataset(cprd_dataset.val, lm.3m.pooled.high_vars_extra, outcome, model_chr = "Highlighted vars + extra", outcome_chr = "3-months", dataset_chr = "Validation", drug_chr = "Pooled")
)

#### By drugclass

# MFN
lm.3m.MFN.high_vars_extra <- glm(formula = as.formula(paste0(outcome, formula_high_vars_extra)), data = cprd_dataset.dev %>% filter(drugclass == "MFN"), family = binomial())

roc_dataset <- rbind(
  roc_dataset,
  calc_roc_dataset(cprd_dataset.dev %>% filter(drugclass == "MFN"), lm.3m.MFN.high_vars_extra, outcome, model_chr = "Highlighted vars + extra", outcome_chr = "3-months", dataset_chr = "Development", drug_chr = "MFN"),
  calc_roc_dataset(cprd_dataset.val %>% filter(drugclass == "MFN"), lm.3m.MFN.high_vars_extra, outcome, model_chr = "Highlighted vars + extra", outcome_chr = "3-months", dataset_chr = "Validation", drug_chr = "MFN")
)

# GLP1
lm.3m.GLP1.high_vars_extra <- glm(formula = as.formula(paste0(outcome, formula_high_vars_extra, "+ stopdrug_3m_3mFU_MFN_hist")), data = cprd_dataset.dev %>% filter(drugclass == "GLP1"), family = binomial())

roc_dataset <- rbind(
  roc_dataset,
  calc_roc_dataset(cprd_dataset.dev %>% filter(drugclass == "GLP1"), lm.3m.GLP1.high_vars_extra, outcome, model_chr = "Highlighted vars + extra", outcome_chr = "3-months", dataset_chr = "Development", drug_chr = "GLP1"),
  calc_roc_dataset(cprd_dataset.val %>% filter(drugclass == "GLP1"), lm.3m.GLP1.high_vars_extra, outcome, model_chr = "Highlighted vars + extra", outcome_chr = "3-months", dataset_chr = "Validation", drug_chr = "GLP1")
)

# DPP4
lm.3m.DPP4.high_vars_extra <- glm(formula = as.formula(paste0(outcome, formula_high_vars_extra, "+ stopdrug_3m_3mFU_MFN_hist")), data = cprd_dataset.dev %>% filter(drugclass == "DPP4"), family = binomial())

roc_dataset <- rbind(
  roc_dataset,
  calc_roc_dataset(cprd_dataset.dev %>% filter(drugclass == "DPP4"), lm.3m.DPP4.high_vars_extra, outcome, model_chr = "Highlighted vars + extra", outcome_chr = "3-months", dataset_chr = "Development", drug_chr = "DPP4"),
  calc_roc_dataset(cprd_dataset.val %>% filter(drugclass == "DPP4"), lm.3m.DPP4.high_vars_extra, outcome, model_chr = "Highlighted vars + extra", outcome_chr = "3-months", dataset_chr = "Validation", drug_chr = "DPP4")
)

# SGLT2
lm.3m.SGLT2.high_vars_extra <- glm(formula = as.formula(paste0(outcome, formula_high_vars_extra, "+ stopdrug_3m_3mFU_MFN_hist")), data = cprd_dataset.dev %>% filter(drugclass == "SGLT2"), family = binomial())

roc_dataset <- rbind(
  roc_dataset,
  calc_roc_dataset(cprd_dataset.dev %>% filter(drugclass == "SGLT2"), lm.3m.SGLT2.high_vars_extra, outcome, model_chr = "Highlighted vars + extra", outcome_chr = "3-months", dataset_chr = "Development", drug_chr = "SGLT2"),
  calc_roc_dataset(cprd_dataset.val %>% filter(drugclass == "SGLT2"), lm.3m.SGLT2.high_vars_extra, outcome, model_chr = "Highlighted vars + extra", outcome_chr = "3-months", dataset_chr = "Validation", drug_chr = "SGLT2")
)

# TZD
lm.3m.TZD.high_vars_extra <- glm(formula = as.formula(paste0(outcome, formula_high_vars_extra, "+ stopdrug_3m_3mFU_MFN_hist")), data = cprd_dataset.dev %>% filter(drugclass == "TZD"), family = binomial())

roc_dataset <- rbind(
  roc_dataset,
  calc_roc_dataset(cprd_dataset.dev %>% filter(drugclass == "TZD"), lm.3m.TZD.high_vars_extra, outcome, model_chr = "Highlighted vars + extra", outcome_chr = "3-months", dataset_chr = "Development", drug_chr = "TZD"),
  calc_roc_dataset(cprd_dataset.val %>% filter(drugclass == "TZD"), lm.3m.TZD.high_vars_extra, outcome, model_chr = "Highlighted vars + extra", outcome_chr = "3-months", dataset_chr = "Validation", drug_chr = "TZD")
)

# SU
lm.3m.SU.high_vars_extra <- glm(formula = as.formula(paste0(outcome, formula_high_vars_extra, "+ stopdrug_3m_3mFU_MFN_hist")), data = cprd_dataset.dev %>% filter(drugclass == "SU"), family = binomial())

roc_dataset <- rbind(
  roc_dataset,
  calc_roc_dataset(cprd_dataset.dev %>% filter(drugclass == "SU"), lm.3m.SU.high_vars_extra, outcome, model_chr = "Highlighted vars + extra", outcome_chr = "3-months", dataset_chr = "Development", drug_chr = "SU"),
  calc_roc_dataset(cprd_dataset.val %>% filter(drugclass == "SU"), lm.3m.SU.high_vars_extra, outcome, model_chr = "Highlighted vars + extra", outcome_chr = "3-months", dataset_chr = "Validation", drug_chr = "SU")
)




###############################################################################
###############################################################################
############################### Read Data In ##################################
###############################################################################
###############################################################################

# load dataset
cprd_dataset.dev <- set_up_data(
  raw_data = "20240308_t2d_1stinstance",
  diagnosis = FALSE,
  therapies = c("MFN", "GLP1", "DPP4", "SGLT2", "TZD", "SU"),
  dataset = "6m.disc.dataset.dev",
  follow_up = "6-months",
  full_prescribing_history = TRUE
) %>%
  drop_na(-stopdrug_12m_6mFU)

# load dataset
cprd_dataset.val <- set_up_data(
  raw_data = "20240308_t2d_1stinstance",
  diagnosis = FALSE,
  therapies = c("MFN", "GLP1", "DPP4", "SGLT2", "TZD", "SU"),
  dataset = "6m.disc.dataset.val",
  follow_up = "6-months",
  full_prescribing_history = TRUE
) %>%
  drop_na(-stopdrug_12m_6mFU)


### Prediction models based on highlighted variables

formula_high_vars <- "~ dstartdate_age + gender + imd2015_10 + prebmi + dstartdate_dm_dur + prehba1c + drugline + predrug_frailty_proxy + ethnicity_5cat"

outcome <- "stopdrug_6m_6mFU"

#### Pooled

lm.6m.pooled.high_vars <- glm(formula = as.formula(paste0(outcome, formula_high_vars)), data = cprd_dataset.dev, family = binomial())

roc_dataset <- rbind(
  roc_dataset,
  calc_roc_dataset(cprd_dataset.dev, lm.6m.pooled.high_vars, outcome, model_chr = "Highlighted vars", outcome_chr = "6-months", dataset_chr = "Development", drug_chr = "Pooled"),
  calc_roc_dataset(cprd_dataset.val, lm.6m.pooled.high_vars, outcome, model_chr = "Highlighted vars", outcome_chr = "6-months", dataset_chr = "Validation", drug_chr = "Pooled")
)

#### By drugclass

# MFN
lm.6m.MFN.high_vars <- glm(formula = as.formula(paste0(outcome, formula_high_vars)), data = cprd_dataset.dev %>% filter(drugclass == "MFN"), family = binomial())

roc_dataset <- rbind(
  roc_dataset,
  calc_roc_dataset(cprd_dataset.dev %>% filter(drugclass == "MFN"), lm.6m.MFN.high_vars, outcome, model_chr = "Highlighted vars", outcome_chr = "6-months", dataset_chr = "Development", drug_chr = "MFN"),
  calc_roc_dataset(cprd_dataset.val %>% filter(drugclass == "MFN"), lm.6m.MFN.high_vars, outcome, model_chr = "Highlighted vars", outcome_chr = "6-months", dataset_chr = "Validation", drug_chr = "MFN")
)

# GLP1
lm.6m.GLP1.high_vars <- glm(formula = as.formula(paste0(outcome, formula_high_vars)), data = cprd_dataset.dev %>% filter(drugclass == "GLP1"), family = binomial())

roc_dataset <- rbind(
  roc_dataset,
  calc_roc_dataset(cprd_dataset.dev %>% filter(drugclass == "GLP1"), lm.6m.GLP1.high_vars, outcome, model_chr = "Highlighted vars", outcome_chr = "6-months", dataset_chr = "Development", drug_chr = "GLP1"),
  calc_roc_dataset(cprd_dataset.val %>% filter(drugclass == "GLP1"), lm.6m.GLP1.high_vars, outcome, model_chr = "Highlighted vars", outcome_chr = "6-months", dataset_chr = "Validation", drug_chr = "GLP1")
)

# DPP4
lm.6m.DPP4.high_vars <- glm(formula = as.formula(paste0(outcome, formula_high_vars)), data = cprd_dataset.dev %>% filter(drugclass == "DPP4"), family = binomial())

roc_dataset <- rbind(
  roc_dataset,
  calc_roc_dataset(cprd_dataset.dev %>% filter(drugclass == "DPP4"), lm.6m.DPP4.high_vars, outcome, model_chr = "Highlighted vars", outcome_chr = "6-months", dataset_chr = "Development", drug_chr = "DPP4"),
  calc_roc_dataset(cprd_dataset.val %>% filter(drugclass == "DPP4"), lm.6m.DPP4.high_vars, outcome, model_chr = "Highlighted vars", outcome_chr = "6-months", dataset_chr = "Validation", drug_chr = "DPP4")
)

# SGLT2
lm.6m.SGLT2.high_vars <- glm(formula = as.formula(paste0(outcome, formula_high_vars)), data = cprd_dataset.dev %>% filter(drugclass == "SGLT2"), family = binomial())

roc_dataset <- rbind(
  roc_dataset,
  calc_roc_dataset(cprd_dataset.dev %>% filter(drugclass == "SGLT2"), lm.6m.SGLT2.high_vars, outcome, model_chr = "Highlighted vars", outcome_chr = "6-months", dataset_chr = "Development", drug_chr = "SGLT2"),
  calc_roc_dataset(cprd_dataset.val %>% filter(drugclass == "SGLT2"), lm.6m.SGLT2.high_vars, outcome, model_chr = "Highlighted vars", outcome_chr = "6-months", dataset_chr = "Validation", drug_chr = "SGLT2")
)

# TZD
lm.6m.TZD.high_vars <- glm(formula = as.formula(paste0(outcome, formula_high_vars)), data = cprd_dataset.dev %>% filter(drugclass == "TZD"), family = binomial())

roc_dataset <- rbind(
  roc_dataset,
  calc_roc_dataset(cprd_dataset.dev %>% filter(drugclass == "TZD"), lm.6m.TZD.high_vars, outcome, model_chr = "Highlighted vars", outcome_chr = "6-months", dataset_chr = "Development", drug_chr = "TZD"),
  calc_roc_dataset(cprd_dataset.val %>% filter(drugclass == "TZD"), lm.6m.TZD.high_vars, outcome, model_chr = "Highlighted vars", outcome_chr = "6-months", dataset_chr = "Validation", drug_chr = "TZD")
)

# SU
lm.6m.SU.high_vars <- glm(formula = as.formula(paste0(outcome, formula_high_vars)), data = cprd_dataset.dev %>% filter(drugclass == "SU"), family = binomial())

roc_dataset <- rbind(
  roc_dataset,
  calc_roc_dataset(cprd_dataset.dev %>% filter(drugclass == "SU"), lm.6m.SU.high_vars, outcome, model_chr = "Highlighted vars", outcome_chr = "6-months", dataset_chr = "Development", drug_chr = "SU"),
  calc_roc_dataset(cprd_dataset.val %>% filter(drugclass == "SU"), lm.6m.SU.high_vars, outcome, model_chr = "Highlighted vars", outcome_chr = "6-months", dataset_chr = "Validation", drug_chr = "SU")
)


### Prediction models based on highlighted variable + additional

formula_high_vars_extra <- "~ dstartdate_age + gender + imd2015_10 + prebmi + dstartdate_dm_dur + prehba1c + drugline + predrug_frailty_proxy + ethnicity_5cat + numdrugs + predrug_bloodmed + smoking_cat + predrug_statins + preegfr + prehdl"

#### Pooled

lm.6m.pooled.high_vars_extra <- glm(formula = as.formula(paste0(outcome, formula_high_vars_extra, "+ drugclass:stopdrug_3m_3mFU_MFN_hist")), data = cprd_dataset.dev, family = binomial())

roc_dataset <- rbind(
  roc_dataset,
  calc_roc_dataset(cprd_dataset.dev, lm.6m.pooled.high_vars_extra, outcome, model_chr = "Highlighted vars + extra", outcome_chr = "6-months", dataset_chr = "Development", drug_chr = "Pooled"),
  calc_roc_dataset(cprd_dataset.val, lm.6m.pooled.high_vars_extra, outcome, model_chr = "Highlighted vars + extra", outcome_chr = "6-months", dataset_chr = "Validation", drug_chr = "Pooled")
)

#### By drugclass

# MFN
lm.6m.MFN.high_vars_extra <- glm(formula = as.formula(paste0(outcome, formula_high_vars_extra)), data = cprd_dataset.dev %>% filter(drugclass == "MFN"), family = binomial())

roc_dataset <- rbind(
  roc_dataset,
  calc_roc_dataset(cprd_dataset.dev %>% filter(drugclass == "MFN"), lm.6m.MFN.high_vars_extra, outcome, model_chr = "Highlighted vars + extra", outcome_chr = "6-months", dataset_chr = "Development", drug_chr = "MFN"),
  calc_roc_dataset(cprd_dataset.val %>% filter(drugclass == "MFN"), lm.6m.MFN.high_vars_extra, outcome, model_chr = "Highlighted vars + extra", outcome_chr = "6-months", dataset_chr = "Validation", drug_chr = "MFN")
)

# GLP1
lm.6m.GLP1.high_vars_extra <- glm(formula = as.formula(paste0(outcome, formula_high_vars_extra, "+ stopdrug_3m_3mFU_MFN_hist")), data = cprd_dataset.dev %>% filter(drugclass == "GLP1"), family = binomial())

roc_dataset <- rbind(
  roc_dataset,
  calc_roc_dataset(cprd_dataset.dev %>% filter(drugclass == "GLP1"), lm.6m.GLP1.high_vars_extra, outcome, model_chr = "Highlighted vars + extra", outcome_chr = "6-months", dataset_chr = "Development", drug_chr = "GLP1"),
  calc_roc_dataset(cprd_dataset.val %>% filter(drugclass == "GLP1"), lm.6m.GLP1.high_vars_extra, outcome, model_chr = "Highlighted vars + extra", outcome_chr = "6-months", dataset_chr = "Validation", drug_chr = "GLP1")
)

# DPP4
lm.6m.DPP4.high_vars_extra <- glm(formula = as.formula(paste0(outcome, formula_high_vars_extra, "+ stopdrug_3m_3mFU_MFN_hist")), data = cprd_dataset.dev %>% filter(drugclass == "DPP4"), family = binomial())

roc_dataset <- rbind(
  roc_dataset,
  calc_roc_dataset(cprd_dataset.dev %>% filter(drugclass == "DPP4"), lm.6m.DPP4.high_vars_extra, outcome, model_chr = "Highlighted vars + extra", outcome_chr = "6-months", dataset_chr = "Development", drug_chr = "DPP4"),
  calc_roc_dataset(cprd_dataset.val %>% filter(drugclass == "DPP4"), lm.6m.DPP4.high_vars_extra, outcome, model_chr = "Highlighted vars + extra", outcome_chr = "6-months", dataset_chr = "Validation", drug_chr = "DPP4")
)

# SGLT2
lm.6m.SGLT2.high_vars_extra <- glm(formula = as.formula(paste0(outcome, formula_high_vars_extra, "+ stopdrug_3m_3mFU_MFN_hist")), data = cprd_dataset.dev %>% filter(drugclass == "SGLT2"), family = binomial())

roc_dataset <- rbind(
  roc_dataset,
  calc_roc_dataset(cprd_dataset.dev %>% filter(drugclass == "SGLT2"), lm.6m.SGLT2.high_vars_extra, outcome, model_chr = "Highlighted vars + extra", outcome_chr = "6-months", dataset_chr = "Development", drug_chr = "SGLT2"),
  calc_roc_dataset(cprd_dataset.val %>% filter(drugclass == "SGLT2"), lm.6m.SGLT2.high_vars_extra, outcome, model_chr = "Highlighted vars + extra", outcome_chr = "6-months", dataset_chr = "Validation", drug_chr = "SGLT2")
)

# TZD
lm.6m.TZD.high_vars_extra <- glm(formula = as.formula(paste0(outcome, formula_high_vars_extra, "+ stopdrug_3m_3mFU_MFN_hist")), data = cprd_dataset.dev %>% filter(drugclass == "TZD"), family = binomial())

roc_dataset <- rbind(
  roc_dataset,
  calc_roc_dataset(cprd_dataset.dev %>% filter(drugclass == "TZD"), lm.6m.TZD.high_vars_extra, outcome, model_chr = "Highlighted vars + extra", outcome_chr = "6-months", dataset_chr = "Development", drug_chr = "TZD"),
  calc_roc_dataset(cprd_dataset.val %>% filter(drugclass == "TZD"), lm.6m.TZD.high_vars_extra, outcome, model_chr = "Highlighted vars + extra", outcome_chr = "6-months", dataset_chr = "Validation", drug_chr = "TZD")
)

# SU
lm.6m.SU.high_vars_extra <- glm(formula = as.formula(paste0(outcome, formula_high_vars_extra, "+ stopdrug_3m_3mFU_MFN_hist")), data = cprd_dataset.dev %>% filter(drugclass == "SU"), family = binomial())

roc_dataset <- rbind(
  roc_dataset,
  calc_roc_dataset(cprd_dataset.dev %>% filter(drugclass == "SU"), lm.6m.SU.high_vars_extra, outcome, model_chr = "Highlighted vars + extra", outcome_chr = "6-months", dataset_chr = "Development", drug_chr = "SU"),
  calc_roc_dataset(cprd_dataset.val %>% filter(drugclass == "SU"), lm.6m.SU.high_vars_extra, outcome, model_chr = "Highlighted vars + extra", outcome_chr = "6-months", dataset_chr = "Validation", drug_chr = "SU")
)



###############################################################################
###############################################################################
############################### Read Data In ##################################
###############################################################################
###############################################################################

# load dataset
cprd_dataset.dev <- set_up_data(
  raw_data = "20240308_t2d_1stinstance",
  diagnosis = FALSE,
  therapies = c("MFN", "GLP1", "DPP4", "SGLT2", "TZD", "SU"),
  dataset = "12m.disc.dataset.dev",
  follow_up = "6-months",
  full_prescribing_history = TRUE
) %>%
  drop_na()

# load dataset
cprd_dataset.val <- set_up_data(
  raw_data = "20240308_t2d_1stinstance",
  diagnosis = FALSE,
  therapies = c("MFN", "GLP1", "DPP4", "SGLT2", "TZD", "SU"),
  dataset = "12m.disc.dataset.val",
  follow_up = "6-months",
  full_prescribing_history = TRUE
) %>%
  drop_na()


### Prediction models based on highlighted variables

formula_high_vars <- "~ dstartdate_age + gender + imd2015_10 + prebmi + dstartdate_dm_dur + prehba1c + drugline + predrug_frailty_proxy + ethnicity_5cat"

outcome <- "stopdrug_12m_6mFU"

#### Pooled

lm.12m.pooled.high_vars <- glm(formula = as.formula(paste0(outcome, formula_high_vars)), data = cprd_dataset.dev, family = binomial())

roc_dataset <- rbind(
  roc_dataset,
  calc_roc_dataset(cprd_dataset.dev, lm.12m.pooled.high_vars, outcome, model_chr = "Highlighted vars", outcome_chr = "12-months", dataset_chr = "Development", drug_chr = "Pooled"),
  calc_roc_dataset(cprd_dataset.val, lm.12m.pooled.high_vars, outcome, model_chr = "Highlighted vars", outcome_chr = "12-months", dataset_chr = "Validation", drug_chr = "Pooled")
)

#### By drugclass

# MFN
lm.12m.MFN.high_vars <- glm(formula = as.formula(paste0(outcome, formula_high_vars)), data = cprd_dataset.dev %>% filter(drugclass == "MFN"), family = binomial())

roc_dataset <- rbind(
  roc_dataset,
  calc_roc_dataset(cprd_dataset.dev %>% filter(drugclass == "MFN"), lm.12m.MFN.high_vars, outcome, model_chr = "Highlighted vars", outcome_chr = "12-months", dataset_chr = "Development", drug_chr = "MFN"),
  calc_roc_dataset(cprd_dataset.val %>% filter(drugclass == "MFN"), lm.12m.MFN.high_vars, outcome, model_chr = "Highlighted vars", outcome_chr = "12-months", dataset_chr = "Validation", drug_chr = "MFN")
)

# GLP1
lm.12m.GLP1.high_vars <- glm(formula = as.formula(paste0(outcome, formula_high_vars)), data = cprd_dataset.dev %>% filter(drugclass == "GLP1"), family = binomial())

roc_dataset <- rbind(
  roc_dataset,
  calc_roc_dataset(cprd_dataset.dev %>% filter(drugclass == "GLP1"), lm.12m.GLP1.high_vars, outcome, model_chr = "Highlighted vars", outcome_chr = "12-months", dataset_chr = "Development", drug_chr = "GLP1"),
  calc_roc_dataset(cprd_dataset.val %>% filter(drugclass == "GLP1"), lm.12m.GLP1.high_vars, outcome, model_chr = "Highlighted vars", outcome_chr = "12-months", dataset_chr = "Validation", drug_chr = "GLP1")
)

# DPP4
lm.12m.DPP4.high_vars <- glm(formula = as.formula(paste0(outcome, formula_high_vars)), data = cprd_dataset.dev %>% filter(drugclass == "DPP4"), family = binomial())

roc_dataset <- rbind(
  roc_dataset,
  calc_roc_dataset(cprd_dataset.dev %>% filter(drugclass == "DPP4"), lm.12m.DPP4.high_vars, outcome, model_chr = "Highlighted vars", outcome_chr = "12-months", dataset_chr = "Development", drug_chr = "DPP4"),
  calc_roc_dataset(cprd_dataset.val %>% filter(drugclass == "DPP4"), lm.12m.DPP4.high_vars, outcome, model_chr = "Highlighted vars", outcome_chr = "12-months", dataset_chr = "Validation", drug_chr = "DPP4")
)

# SGLT2
lm.12m.SGLT2.high_vars <- glm(formula = as.formula(paste0(outcome, formula_high_vars)), data = cprd_dataset.dev %>% filter(drugclass == "SGLT2"), family = binomial())

roc_dataset <- rbind(
  roc_dataset,
  calc_roc_dataset(cprd_dataset.dev %>% filter(drugclass == "SGLT2"), lm.12m.SGLT2.high_vars, outcome, model_chr = "Highlighted vars", outcome_chr = "12-months", dataset_chr = "Development", drug_chr = "SGLT2"),
  calc_roc_dataset(cprd_dataset.val %>% filter(drugclass == "SGLT2"), lm.12m.SGLT2.high_vars, outcome, model_chr = "Highlighted vars", outcome_chr = "12-months", dataset_chr = "Validation", drug_chr = "SGLT2")
)

# TZD
lm.12m.TZD.high_vars <- glm(formula = as.formula(paste0(outcome, formula_high_vars)), data = cprd_dataset.dev %>% filter(drugclass == "TZD"), family = binomial())

roc_dataset <- rbind(
  roc_dataset,
  calc_roc_dataset(cprd_dataset.dev %>% filter(drugclass == "TZD"), lm.12m.TZD.high_vars, outcome, model_chr = "Highlighted vars", outcome_chr = "12-months", dataset_chr = "Development", drug_chr = "TZD"),
  calc_roc_dataset(cprd_dataset.val %>% filter(drugclass == "TZD"), lm.12m.TZD.high_vars, outcome, model_chr = "Highlighted vars", outcome_chr = "12-months", dataset_chr = "Validation", drug_chr = "TZD")
)

# SU
lm.12m.SU.high_vars <- glm(formula = as.formula(paste0(outcome, formula_high_vars)), data = cprd_dataset.dev %>% filter(drugclass == "SU"), family = binomial())

roc_dataset <- rbind(
  roc_dataset,
  calc_roc_dataset(cprd_dataset.dev %>% filter(drugclass == "SU"), lm.12m.SU.high_vars, outcome, model_chr = "Highlighted vars", outcome_chr = "12-months", dataset_chr = "Development", drug_chr = "SU"),
  calc_roc_dataset(cprd_dataset.val %>% filter(drugclass == "SU"), lm.12m.SU.high_vars, outcome, model_chr = "Highlighted vars", outcome_chr = "12-months", dataset_chr = "Validation", drug_chr = "SU")
)


### Prediction models based on highlighted variable + additional

formula_high_vars_extra <- "~ dstartdate_age + gender + imd2015_10 + prebmi + dstartdate_dm_dur + prehba1c + drugline + predrug_frailty_proxy + ethnicity_5cat + numdrugs + predrug_bloodmed + smoking_cat + predrug_statins + preegfr + prehdl"

#### Pooled

lm.12m.pooled.high_vars_extra <- glm(formula = as.formula(paste0(outcome, formula_high_vars_extra, "+ drugclass:stopdrug_3m_3mFU_MFN_hist")), data = cprd_dataset.dev, family = binomial())

roc_dataset <- rbind(
  roc_dataset,
  calc_roc_dataset(cprd_dataset.dev, lm.12m.pooled.high_vars_extra, outcome, model_chr = "Highlighted vars + extra", outcome_chr = "12-months", dataset_chr = "Development", drug_chr = "Pooled"),
  calc_roc_dataset(cprd_dataset.val, lm.12m.pooled.high_vars_extra, outcome, model_chr = "Highlighted vars + extra", outcome_chr = "12-months", dataset_chr = "Validation", drug_chr = "Pooled")
)

#### By drugclass

# MFN
lm.12m.MFN.high_vars_extra <- glm(formula = as.formula(paste0(outcome, formula_high_vars_extra)), data = cprd_dataset.dev %>% filter(drugclass == "MFN"), family = binomial())

roc_dataset <- rbind(
  roc_dataset,
  calc_roc_dataset(cprd_dataset.dev %>% filter(drugclass == "MFN"), lm.12m.MFN.high_vars_extra, outcome, model_chr = "Highlighted vars + extra", outcome_chr = "12-months", dataset_chr = "Development", drug_chr = "MFN"),
  calc_roc_dataset(cprd_dataset.val %>% filter(drugclass == "MFN"), lm.12m.MFN.high_vars_extra, outcome, model_chr = "Highlighted vars + extra", outcome_chr = "12-months", dataset_chr = "Validation", drug_chr = "MFN")
)

# GLP1
lm.12m.GLP1.high_vars_extra <- glm(formula = as.formula(paste0(outcome, formula_high_vars_extra, "+ stopdrug_3m_3mFU_MFN_hist")), data = cprd_dataset.dev %>% filter(drugclass == "GLP1"), family = binomial())

roc_dataset <- rbind(
  roc_dataset,
  calc_roc_dataset(cprd_dataset.dev %>% filter(drugclass == "GLP1"), lm.12m.GLP1.high_vars_extra, outcome, model_chr = "Highlighted vars + extra", outcome_chr = "12-months", dataset_chr = "Development", drug_chr = "GLP1"),
  calc_roc_dataset(cprd_dataset.val %>% filter(drugclass == "GLP1"), lm.12m.GLP1.high_vars_extra, outcome, model_chr = "Highlighted vars + extra", outcome_chr = "12-months", dataset_chr = "Validation", drug_chr = "GLP1")
)

# DPP4
lm.12m.DPP4.high_vars_extra <- glm(formula = as.formula(paste0(outcome, formula_high_vars_extra, "+ stopdrug_3m_3mFU_MFN_hist")), data = cprd_dataset.dev %>% filter(drugclass == "DPP4"), family = binomial())

roc_dataset <- rbind(
  roc_dataset,
  calc_roc_dataset(cprd_dataset.dev %>% filter(drugclass == "DPP4"), lm.12m.DPP4.high_vars_extra, outcome, model_chr = "Highlighted vars + extra", outcome_chr = "12-months", dataset_chr = "Development", drug_chr = "DPP4"),
  calc_roc_dataset(cprd_dataset.val %>% filter(drugclass == "DPP4"), lm.12m.DPP4.high_vars_extra, outcome, model_chr = "Highlighted vars + extra", outcome_chr = "12-months", dataset_chr = "Validation", drug_chr = "DPP4")
)

# SGLT2
lm.12m.SGLT2.high_vars_extra <- glm(formula = as.formula(paste0(outcome, formula_high_vars_extra, "+ stopdrug_3m_3mFU_MFN_hist")), data = cprd_dataset.dev %>% filter(drugclass == "SGLT2"), family = binomial())

roc_dataset <- rbind(
  roc_dataset,
  calc_roc_dataset(cprd_dataset.dev %>% filter(drugclass == "SGLT2"), lm.12m.SGLT2.high_vars_extra, outcome, model_chr = "Highlighted vars + extra", outcome_chr = "12-months", dataset_chr = "Development", drug_chr = "SGLT2"),
  calc_roc_dataset(cprd_dataset.val %>% filter(drugclass == "SGLT2"), lm.12m.SGLT2.high_vars_extra, outcome, model_chr = "Highlighted vars + extra", outcome_chr = "12-months", dataset_chr = "Validation", drug_chr = "SGLT2")
)

# TZD
lm.12m.TZD.high_vars_extra <- glm(formula = as.formula(paste0(outcome, formula_high_vars_extra, "+ stopdrug_3m_3mFU_MFN_hist")), data = cprd_dataset.dev %>% filter(drugclass == "TZD"), family = binomial())

roc_dataset <- rbind(
  roc_dataset,
  calc_roc_dataset(cprd_dataset.dev %>% filter(drugclass == "TZD"), lm.12m.TZD.high_vars_extra, outcome, model_chr = "Highlighted vars + extra", outcome_chr = "12-months", dataset_chr = "Development", drug_chr = "TZD"),
  calc_roc_dataset(cprd_dataset.val %>% filter(drugclass == "TZD"), lm.12m.TZD.high_vars_extra, outcome, model_chr = "Highlighted vars + extra", outcome_chr = "12-months", dataset_chr = "Validation", drug_chr = "TZD")
)

# SU
lm.12m.SU.high_vars_extra <- glm(formula = as.formula(paste0(outcome, formula_high_vars_extra, "+ stopdrug_3m_3mFU_MFN_hist")), data = cprd_dataset.dev %>% filter(drugclass == "SU"), family = binomial())

roc_dataset <- rbind(
  roc_dataset,
  calc_roc_dataset(cprd_dataset.dev %>% filter(drugclass == "SU"), lm.12m.SU.high_vars_extra, outcome, model_chr = "Highlighted vars + extra", outcome_chr = "12-months", dataset_chr = "Development", drug_chr = "SU"),
  calc_roc_dataset(cprd_dataset.val %>% filter(drugclass == "SU"), lm.12m.SU.high_vars_extra, outcome, model_chr = "Highlighted vars + extra", outcome_chr = "12-months", dataset_chr = "Validation", drug_chr = "SU")
)


saveRDS(roc_dataset, "results/Models/Predictions/04.roc_multivariate_disc.rds")



pdf("results/figures/04.auc_multivariate_disc.pdf", width = 12, height = 6)

roc_dataset %>%
  mutate(
    model = factor(model, levels = c("Highlighted vars", "Highlighted vars + extra"), labels = c("Highlighted features", "Highlighted features + routine clinical features + biomarkers")),
    outcome = factor(outcome, levels = c("12-months", "6-months", "3-months")),
    dataset = factor(dataset, levels = c("Development", "Validation"), labels = c("Development dataset", "Validation dataset")),
    drug = factor(drug, levels = rev(c("Pooled", "MFN", "GLP1", "DPP4", "SGLT2", "TZD", "SU"))),
  ) %>%
  ggplot(aes(y = outcome, x = roc, xmin = roc_lci, xmax = roc_uci, colour = drug)) +
  geom_point(position = position_dodge(width = 1)) + 
  geom_errorbar(position = position_dodge(width = 1)) +
  scale_x_continuous("AUROC", limits = c(0.475, 0.65), breaks = seq(0.475, 0.65, 0.025)) +
  scale_colour_manual(values = c("Pooled" = "black", "SGLT2" = "#E69F00", "GLP1" = "#56B4E9", "SU" = "#CC79A7", "DPP4" = "#0072B2", "TZD" = "#D55E00", "MFN" = "grey"), breaks = c("Pooled", "MFN", "GLP1", "DPP4", "SGLT2", "TZD", "SU"), labels = c("Pooled", "Metformin", "GLP-1RA", "DPP4i", "SGLT2i", "TZD", "SU"), name = "Therapy", guide = guide_legend(reverse = TRUE)) +
  guides(colour = guide_legend(nrow = 2, byrow = TRUE)) +
  facet_grid(dataset ~ model) +
  theme_bw() +
  theme(
    legend.position = "bottom",
    axis.title.y = element_blank(),
    axis.title.x = element_text(size = 15),
    axis.text = element_text(size = 13),
    strip.text = element_text(size = 13),
    legend.title = element_text(size = 13),
    legend.text = element_text(size = 13),
    panel.spacing = unit(0.3, "cm", data = NULL)
  )

dev.off()

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
  dataset = "3m.disc.dataset.val",
  follow_up = "6-months",
  full_prescribing_history = TRUE
) %>%
  drop_na(-stopdrug_6m_6mFU, -stopdrug_12m_6mFU)


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
                            prehdl*drugclass + 
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


# model_selection <- step(
#   glm(stopdrug_3m_6mFU ~ 1, data = cprd_dataset.dev, family = binomial()),
#   scope = formula_glm,
#   direction = "both"
#   )

#:----- What was found:
# Step:  AIC=99534.68
# stopdrug_3m_6mFU ~ drugclass + prehdl + gender +
#   stopdrug_3m_3mFU_MFN_hist + drugline + predrug_bloodmed +
#   prebmi + predrug_statins + dstartdate_age + ethnicity_5cat +
#   predrug_cardio_event + prealt + numdrugs + prehba1c + smoking_cat +
#   imd2015_10 + dstartdate_dm_dur + predrug_frailty_proxy +
#   drugclass:gender + drugclass:prehdl + drugclass:prebmi +
#   drugclass:dstartdate_age + drugclass:prealt + drugclass:prehba1c +
#   drugclass:dstartdate_dm_dur



lm_uber_model.no_weight_adjust <- glm(as.formula(formula_glm), data = cprd_dataset.dev, family = binomial())


roc(cprd_dataset.dev %>% select(stopdrug_3m_6mFU) %>% unlist(), predict(lm_uber_model.no_weight_adjust, type = "response"), ci = TRUE)



roc(cprd_dataset.val %>% select(stopdrug_3m_6mFU) %>% unlist(), predict(lm_uber_model.no_weight_adjust, newdata = cprd_dataset.val, type = "response"), ci = TRUE)



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

saveRDS(probabilities.only_dataset, "results/Models/Predictions/04.model_predictions_3m_all_drugs.rds")


cprd_dataset.dev <- cprd_dataset.dev %>%
  mutate(
    pred.no_weight = predict(lm_uber_model.no_weight_adjust, type = "response")
  )

cprd_dataset.val <- cprd_dataset.val %>%
  mutate(
    pred.no_weight = predict(lm_uber_model.no_weight_adjust, newdata = cprd_dataset.val, type = "response")
  )

pdf("results/figures/04.uber_model_calibration_3m_all_drugs.pdf", width = 4, height = 4)

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
  dataset = "3m.disc.dataset.dev",
  follow_up = "6-months",
  full_prescribing_history = TRUE
) %>%
  drop_na(-stopdrug_6m_6mFU, -stopdrug_12m_6mFU)


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
  dataset = "3m.disc.dataset.val",
  follow_up = "6-months",
  full_prescribing_history = TRUE
) %>%
  drop_na(-stopdrug_6m_6mFU, -stopdrug_12m_6mFU)


# modify some variables for fitting
cprd_dataset.val <- cprd_dataset.val %>%
  mutate(
    ethnicity_5cat = factor(ethnicity_5cat, levels = c(0, 1, 2, 3, 4), labels = c("White", "Other", "Other" , "Other" , "Other")),
    smoking_cat = factor(smoking_cat, levels = c("Active smoker", "Ex-smoker", "Non-smoker"), labels = c("Active smoker", "Other", "Other")),
    alcohol_cat = factor(alcohol_cat, levels = c("Excess", "Harmful", "None", "Within limits"), labels = c("Excess/Harmful", "Excess/Harmful", "Within limits/None", "Within limits/None")),
    preckdstage = factor(preckdstage, levels = c("stage_0", "stage_1", "stage_2", "stage_3a", "stage_3b", "stage_4", "stage_5"), labels = c("stage 0/1/2", "stage 0/1/2", "stage 0/1/2", "stage 3/4/5", "stage 3/4/5", "stage 3/4/5", "stage 3/4/5"))
  )


lm_uber_model.no_weight_adjust <- glm(as.formula(formula_glm), data = cprd_dataset.dev, family = binomial())



# model_selection <- step(
#   glm(stopdrug_3m_6mFU ~ 1, data = cprd_dataset.dev, family = binomial()),
#   scope = formula_glm,
#   direction = "both"
#   )

#:----- What was found:
# Step:  AIC=84036.86
# stopdrug_3m_6mFU ~ prehdl + drugclass + stopdrug_3m_3mFU_MFN_hist +
#   ethnicity_5cat + predrug_bloodmed + gender + dstartdate_dm_dur +
#   predrug_statins + dstartdate_age + drugline + numdrugs +
#   prehba1c + predrug_cardio_event + preckdstage + prealt +
#   prebmi + imd2015_10 + predrug_heart_event + drugclass:dstartdate_dm_dur +
#   drugclass:dstartdate_age + drugclass:gender + drugclass:prehba1c +
#   drugclass:prealt + drugclass:prebmi




roc(cprd_dataset.dev %>% select(stopdrug_3m_6mFU) %>% unlist(), predict(lm_uber_model.no_weight_adjust, type = "response"), ci = TRUE)


roc(cprd_dataset.val %>% select(stopdrug_3m_6mFU) %>% unlist(), predict(lm_uber_model.no_weight_adjust, newdata = cprd_dataset.val, type = "response"), ci = TRUE)


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

saveRDS(probabilities.only_dataset, "results/Models/Predictions/04.model_predictions_3m_2nd_line_drugs.rds")


cprd_dataset.dev <- cprd_dataset.dev %>%
  mutate(
    pred.no_weight = predict(lm_uber_model.no_weight_adjust, type = "response")
  )

cprd_dataset.val <- cprd_dataset.val %>%
  mutate(
    pred.no_weight = predict(lm_uber_model.no_weight_adjust, newdata = cprd_dataset.val, type = "response")
  )

pdf("results/figures/04.uber_model_calibration_3m_2nd_line_drugs.pdf", width = 4, height = 4)

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




