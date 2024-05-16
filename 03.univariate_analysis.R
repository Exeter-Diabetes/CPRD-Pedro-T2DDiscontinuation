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
  dataset = "3m.disc.dataset"
) %>%
  drop_na(-stopdrug_6m_6mFU, -stopdrug_12m_6mFU)


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
########################### Univariate analysis ###############################
###############################################################################
###############################################################################


univariate_analysis <- function(data, drugs, outcome = "stopdrug_3m_6mFU", variable, variable_name, type, prop.score.variable = NULL) {
  
  # set up formula
  formula_chr <- paste(outcome, "~", variable)
  
  # output measures
  coef <- NULL; LCI <- NULL; UCI <- NULL
  
  # set up propensity score
  if (!is.null(prop.score.variable)) {
    data <- data %>%
      rename("prop.score.model" = prop.score.variable)
  }
  
  if (!is.null(prop.score.variable)) {
    
    # Pooled model
    lm_all <- glm(formula(formula_chr), data = data, weights = prop.score.model, family = quasibinomial())
    
    coef = c(coef, lm_all$coef[2]); LCI = c(LCI, confint(lm_all)[2,1]); UCI = c(UCI, confint(lm_all)[2,2])
    
    for (drug_i in drugs) {
      
      # Drug model
      lm_drug <- glm(formula(formula_chr), data = data %>% filter(drugclass == drug_i), weights = prop.score.model, family = quasibinomial())
      
      coef = c(coef, lm_drug$coef[2]); LCI = c(LCI, confint(lm_drug)[2,1]); UCI = c(UCI, confint(lm_drug)[2,2])
      
    }
    
  } else {
    
    # Pooled model
    lm_all <- glm(formula(formula_chr), data = data, family = binomial())
    
    coef = c(coef, lm_all$coef[2]); LCI = c(LCI, confint(lm_all)[2,1]); UCI = c(UCI, confint(lm_all)[2,2])
    
    for (drug_i in drugs) {
      
      # Drug model
      lm_drug <- glm(formula(formula_chr), data = data %>% filter(drugclass == drug_i), family = binomial())
      
      coef = c(coef, lm_drug$coef[2]); LCI = c(LCI, confint(lm_drug)[2,1]); UCI = c(UCI, confint(lm_drug)[2,2])
      
    }
    
  }
  
  output <- data.frame(
    cbind(
      variable = rep(variable_name, length(drugs) + 1),
      type = rep(type, length(drugs) + 1),
      drug = c("Pooled", drugs),
      coef = as.numeric(coef),
      LCI = as.numeric(LCI),
      UCI = as.numeric(UCI)
    )
  )
 
  return(output)
   
}



###############################################################################
###############################################################################
################# Univariate analysis - 3m discontinuation ####################
###############################################################################
###############################################################################


drugs = c("MFN", "DPP4", "GLP1", "SGLT2", "SU", "TZD")


# dstartdate_dm_dur

cprd_dataset <- cprd_dataset %>%
  mutate_at(c("dstartdate_dm_dur"), ~(scale(.) %>% as.vector))

dstartdate_dm_dur_sum <- univariate_analysis(
  cprd_dataset, drugs = drugs, outcome = "stopdrug_3m_6mFU", 
  variable = "dstartdate_dm_dur", variable_name = "Diabetes duration (per SD)", type = "Clinical features & biomarkers"
)


# dstartdate_age

cprd_dataset <- cprd_dataset %>%
  mutate_at(c("dstartdate_age"), ~(scale(.) %>% as.vector))

dstartdate_age_sum <- univariate_analysis(
  cprd_dataset, drugs = drugs, outcome = "stopdrug_3m_6mFU", 
  variable = "dstartdate_age", variable_name = "Age (per SD)", type = "Clinical features & biomarkers"
)


# prehba1c

cprd_dataset <- cprd_dataset %>%
  mutate_at(c("prehba1c"), ~(scale(.) %>% as.vector))

prehba1c_sum <- univariate_analysis(
  cprd_dataset, drugs = drugs, outcome = "stopdrug_3m_6mFU", 
  variable = "prehba1c", variable_name = "HbA1c (per SD)", type = "Clinical features & biomarkers"
)


# preegfr

cprd <- cprd_dataset %>%
  mutate_at(c("preegfr"), ~(scale(.) %>% as.vector))

preegfr_sum <- univariate_analysis(
  cprd_dataset, drugs = drugs, outcome = "stopdrug_3m_6mFU", 
  variable = "preegfr", variable_name = "eGFR (per SD)", type = "Clinical features & biomarkers"
)


# prebmi

cprd_dataset <- cprd_dataset %>%
  mutate_at(c("prebmi"), ~(scale(.) %>% as.vector))

prebmi_sum <- univariate_analysis(
  cprd_dataset, drugs = drugs, outcome = "stopdrug_3m_6mFU", 
  variable = "prebmi", variable_name = "BMI (per SD)", type = "Clinical features & biomarkers"
)

# prealt

cprd <- cprd_dataset %>%
  mutate_at(c("prealt"), ~(scale(.) %>% as.vector))

prealt_sum <- univariate_analysis(
  cprd_dataset, drugs = drugs, outcome = "stopdrug_3m_6mFU", 
  variable = "prealt", variable_name = "ALT (per SD)", type = "Clinical features & biomarkers"
)


# pretotalcholesterol

cprd <- cprd_dataset %>%
  mutate_at(c("pretotalcholesterol"), ~(scale(.) %>% as.vector))

pretotalcholesterol_sum <- univariate_analysis(
  cprd_dataset, drugs = drugs, outcome = "stopdrug_3m_6mFU", 
  variable = "pretotalcholesterol", variable_name = "Total cholesterol (per SD)", type = "Clinical features & biomarkers"
)


# gender

gender_sum <- univariate_analysis(
  cprd_dataset, drugs = drugs, outcome = "stopdrug_3m_6mFU", 
  variable = "gender", variable_name = "Sex (ref Male)", type = "Clinical features & biomarkers"
)


# predrug_bloodmed

predrug_bloodmed_sum <- univariate_analysis(
  cprd_dataset, drugs = drugs, outcome = "stopdrug_3m_6mFU", 
  variable = "predrug_bloodmed", variable_name = "Blood pressure medication", type = "Behavioural"
)


# predrug_statins

predrug_statins_sum <- univariate_analysis(
  cprd_dataset, drugs = drugs, outcome = "stopdrug_3m_6mFU", 
  variable = "predrug_statins", variable_name = "Statins", type = "Behavioural"
)


# stopdrug_3m_3mFU_MFN_hist

stopdrug_3m_3mFU_MFN_hist_sum <- univariate_analysis(
  cprd_dataset, drugs = c("DPP4", "GLP1", "SGLT2", "SU", "TZD"), outcome = "stopdrug_3m_6mFU", 
  variable = "stopdrug_3m_3mFU_MFN_hist", variable_name = "MFN discontinuation", type = "Behavioural"
)


# ethnicity

cprd_dataset <- cprd_dataset %>%
  mutate(ethnicity_5cat = factor(ethnicity_5cat, levels = c(0, 1, 2, 3, 4), labels = c("White", "Other", "Other" , "Other" , "Other")))

ethnicity_5cat_sum <- univariate_analysis(
  cprd_dataset, drugs = drugs, outcome = "stopdrug_3m_6mFU", 
  variable = "ethnicity_5cat", variable_name = "Ethnicity (ref White)", type = "Behavioural"
)


# smoking_cat

cprd_dataset <- cprd_dataset %>%
  mutate(smoking_cat = factor(smoking_cat, levels = c("Active smoker", "Ex-smoker", "Non-smoker"), labels = c("Active smoker", "Other", "Other")))

smoking_cat_sum <- univariate_analysis(
  cprd_dataset, drugs = drugs, outcome = "stopdrug_3m_6mFU", 
  variable = "smoking_cat", variable_name = "Smoking (ref Active Smoker)", type = "Behavioural"
)


# alcohol_cat

cprd_dataset <- cprd_dataset %>%
  mutate(alcohol_cat = factor(alcohol_cat, levels = c("Excess", "Harmful", "None", "Within limits"), labels = c("Excess/Harmful", "Excess/Harmful", "Within limits/None", "Within limits/None")))

alcohol_cat_sum <- univariate_analysis(
  cprd_dataset, drugs = drugs, outcome = "stopdrug_3m_6mFU", 
  variable = "alcohol_cat", variable_name = "Alcohol (ref Excess)", type = "Behavioural"
)


# predrug_cardio_event

predrug_cardio_event_sum <- univariate_analysis(
  cprd_dataset, drugs = drugs, outcome = "stopdrug_3m_6mFU", 
  variable = "predrug_cardio_event", variable_name = "Cardiovascular", type = "Comorbidity event"
)


# predrug_heart_event

predrug_heart_event_sum <- univariate_analysis(
  cprd_dataset, drugs = drugs, outcome = "stopdrug_3m_6mFU", 
  variable = "predrug_heart_event", variable_name = "Heart", type = "Comorbidity event"
)


# predrug_micro_event

predrug_micro_event_sum <- univariate_analysis(
  cprd_dataset, drugs = drugs, outcome = "stopdrug_3m_6mFU", 
  variable = "predrug_micro_event", variable_name = "Microvascular", type = "Comorbidity event"
)


# preckdstage

cprd_dataset <- cprd_dataset %>%
  mutate(preckdstage = factor(preckdstage, levels = c("stage_0", "stage_1", "stage_2", "stage_3a", "stage_3b", "stage_4", "stage_5"), labels = c("stage 0/1/2", "stage 0/1/2", "stage 0/1/2", "stage 3/4/5", "stage 3/4/5", "stage 3/4/5", "stage 3/4/5")))

preckdstage_sum <- univariate_analysis(
  cprd_dataset, drugs = drugs, outcome = "stopdrug_3m_6mFU", 
  variable = "preckdstage", variable_name = "CKD stage 3/4/5 (ref stage 0/1/2)", type = "Comorbidity event"
)


# predrug_cld

preckdstage_sum <- univariate_analysis(
  cprd_dataset, drugs = drugs, outcome = "stopdrug_3m_6mFU", 
  variable = "preckdstage", variable_name = "CKD stage 3/4/5 (ref stage 0/1/2)", type = "Comorbidity event"
)


# predrug_frailty_proxy

predrug_frailty_proxy_sum <- univariate_analysis(
  cprd_dataset, drugs = drugs, outcome = "stopdrug_3m_6mFU", 
  variable = "predrug_frailty_proxy", variable_name = "Frailty proxy", type = "Comorbidity event"
)



coefficients <- rbind(
  dstartdate_dm_dur_sum,
  dstartdate_age_sum,
  prehba1c_sum,
  preegfr_sum,
  prebmi_sum,
  prealt_sum,
  pretotalcholesterol_sum,
  gender_sum,
  predrug_bloodmed_sum,
  predrug_statins_sum,
  stopdrug_3m_3mFU_MFN_hist_sum,
  ethnicity_5cat_sum,
  smoking_cat_sum,
  alcohol_cat_sum,
  predrug_cardio_event_sum,
  predrug_heart_event_sum,
  predrug_micro_event_sum,
  preckdstage_sum,
  preckdstage_sum,
  predrug_frailty_proxy_sum
) %>%
  data.frame() %>%
  mutate(
    variable = factor(variable),
    drug = factor(drug, levels = rev(c("Pooled", drugs))),
    type = factor(type, levels = c("Clinical features & biomarkers", "Behavioural", "Comorbidity event")),
    coef = exp(as.numeric(coef)),
    LCI = exp(as.numeric(LCI)),
    UCI = exp(as.numeric(UCI))
  )


pdf("results/figures/univariate_analysis_3m.pdf", width = 12, height = 11)

wrap_plots(
  
  coefficients %>%
    filter(type == "Behavioural") %>%
    ggplot() + 
    geom_vline(aes(xintercept = 1), colour = "black") +
    geom_errorbar(aes(x = coef, xmin = LCI, xmax = UCI, y = variable, colour = drug),position = position_dodge(width = 1)) +
    geom_pointrange(aes(x = coef, xmin = LCI, xmax = UCI, y = variable, colour = drug), position = position_dodge(width = 1)) +
    xlab("Odds Ratio") +
    scale_x_log10(breaks = c(0.25, 0.5, 1, 2, 4), limits = c(0.4, 4.1)) +
    facet_grid(variable~., scales = "free_y") +
    scale_colour_manual(values = c("Pooled" = "black", "SGLT2" = "#E69F00", "GLP1" = "#56B4E9", "SU" = "#CC79A7", "DPP4" = "#0072B2", "TZD" = "#D55E00", "MFN" = "grey"), name = "Therapy", guide = guide_legend(reverse = TRUE)) +
    theme_bw() +
    theme(
      axis.title.y = element_blank(),
      legend.position = "bottom",
      strip.text = element_blank()
    ),
  
  
  coefficients %>%
    filter(type == "Clinical features & biomarkers") %>%
    ggplot() + 
    geom_vline(aes(xintercept = 1), colour = "black") +
    geom_errorbar(aes(x = coef, xmin = LCI, xmax = UCI, y = variable, colour = drug),position = position_dodge(width = 1)) +
    geom_pointrange(aes(x = coef, xmin = LCI, xmax = UCI, y = variable, colour = drug), position = position_dodge(width = 1)) +
    xlab("Odds Ratio") +
    scale_x_log10(breaks = c(0.25, 0.5, 1, 2, 4), limits = c(0.4, 4.1)) +
    facet_grid(variable~., scales = "free_y") +
    
    scale_colour_manual(values = c("Pooled" = "black", "SGLT2" = "#E69F00", "GLP1" = "#56B4E9", "SU" = "#CC79A7", "DPP4" = "#0072B2", "TZD" = "#D55E00", "MFN" = "grey"), name = "Therapy", guide = guide_legend(reverse = TRUE)) +
    theme_bw() +
    theme(
      axis.title.y = element_blank(),
      legend.position = "bottom",
      strip.text = element_blank()
    ),
  
  coefficients %>%
    filter(type == "Comorbidity event") %>%
    ggplot() + 
    geom_vline(aes(xintercept = 1), colour = "black") +
    geom_errorbar(aes(x = coef, xmin = LCI, xmax = UCI, y = variable, colour = drug),position = position_dodge(width = 1)) +
    geom_pointrange(aes(x = coef, xmin = LCI, xmax = UCI, y = variable, colour = drug), position = position_dodge(width = 1)) +
    xlab("Odds Ratio") +
    scale_x_log10(breaks = c(0.25, 0.5, 1, 2, 4), limits = c(0.4, 4.1)) +
    facet_grid(variable~., scales = "free_y") +
    
    scale_colour_manual(values = c("Pooled" = "black", "SGLT2" = "#E69F00", "GLP1" = "#56B4E9", "SU" = "#CC79A7", "DPP4" = "#0072B2", "TZD" = "#D55E00", "MFN" = "grey"), name = "Therapy", guide = guide_legend(reverse = TRUE)) +
    theme_bw() +
    theme(
      axis.title.y = element_blank(),
      legend.position = "bottom",
      strip.text = element_blank()
    ),
  
  ncol = 1, nrow = 3
  
) + 
  plot_layout(guides = "collect", height = c(6, 8, 6), axis_titles = "collect", axes = "collect") +
  plot_annotation(tag_levels = list(c("Behavioural", "Clinical features & biomarkers", "Comorbidity event"))) & 
  theme(
    legend.direction = "vertical",
    axis.text.y = element_text(size = 15),
    axis.text.x = element_text(size = 15),
    axis.title.x = element_text(size = 18)
  )

dev.off()



###############################################################################
###############################################################################
################# Univariate analysis - 6m discontinuation ####################
###############################################################################
###############################################################################


# load dataset
cprd_dataset <- cprd_dataset %>%
  drop_na(stopdrug_6m_6mFU)


drugs = c("MFN", "DPP4", "GLP1", "SGLT2", "SU", "TZD")


# dstartdate_dm_dur

cprd_dataset <- cprd_dataset %>%
  mutate_at(c("dstartdate_dm_dur"), ~(scale(.) %>% as.vector))

dstartdate_dm_dur_sum <- univariate_analysis(
  cprd_dataset, drugs = drugs, outcome = "stopdrug_6m_6mFU", 
  variable = "dstartdate_dm_dur", variable_name = "Diabetes duration (per SD)", type = "Clinical features & biomarkers"
)


# dstartdate_age

cprd_dataset <- cprd_dataset %>%
  mutate_at(c("dstartdate_age"), ~(scale(.) %>% as.vector))

dstartdate_age_sum <- univariate_analysis(
  cprd_dataset, drugs = drugs, outcome = "stopdrug_6m_6mFU", 
  variable = "dstartdate_age", variable_name = "Age (per SD)", type = "Clinical features & biomarkers"
)


# prehba1c

cprd_dataset <- cprd_dataset %>%
  mutate_at(c("prehba1c"), ~(scale(.) %>% as.vector))

prehba1c_sum <- univariate_analysis(
  cprd_dataset, drugs = drugs, outcome = "stopdrug_6m_6mFU", 
  variable = "prehba1c", variable_name = "HbA1c (per SD)", type = "Clinical features & biomarkers"
)


# preegfr

cprd <- cprd_dataset %>%
  mutate_at(c("preegfr"), ~(scale(.) %>% as.vector))

preegfr_sum <- univariate_analysis(
  cprd_dataset, drugs = drugs, outcome = "stopdrug_6m_6mFU", 
  variable = "preegfr", variable_name = "eGFR (per SD)", type = "Clinical features & biomarkers"
)


# prebmi

cprd_dataset <- cprd_dataset %>%
  mutate_at(c("prebmi"), ~(scale(.) %>% as.vector))

prebmi_sum <- univariate_analysis(
  cprd_dataset, drugs = drugs, outcome = "stopdrug_6m_6mFU", 
  variable = "prebmi", variable_name = "BMI (per SD)", type = "Clinical features & biomarkers"
)

# prealt

cprd <- cprd_dataset %>%
  mutate_at(c("prealt"), ~(scale(.) %>% as.vector))

prealt_sum <- univariate_analysis(
  cprd_dataset, drugs = drugs, outcome = "stopdrug_6m_6mFU", 
  variable = "prealt", variable_name = "ALT (per SD)", type = "Clinical features & biomarkers"
)


# pretotalcholesterol

cprd <- cprd_dataset %>%
  mutate_at(c("pretotalcholesterol"), ~(scale(.) %>% as.vector))

pretotalcholesterol_sum <- univariate_analysis(
  cprd_dataset, drugs = drugs, outcome = "stopdrug_6m_6mFU", 
  variable = "pretotalcholesterol", variable_name = "Total cholesterol (per SD)", type = "Clinical features & biomarkers"
)


# gender

gender_sum <- univariate_analysis(
  cprd_dataset, drugs = drugs, outcome = "stopdrug_6m_6mFU", 
  variable = "gender", variable_name = "Sex (ref Male)", type = "Clinical features & biomarkers"
)


# predrug_bloodmed

predrug_bloodmed_sum <- univariate_analysis(
  cprd_dataset, drugs = drugs, outcome = "stopdrug_6m_6mFU", 
  variable = "predrug_bloodmed", variable_name = "Blood pressure medication", type = "Behavioural"
)


# predrug_statins

predrug_statins_sum <- univariate_analysis(
  cprd_dataset, drugs = drugs, outcome = "stopdrug_6m_6mFU", 
  variable = "predrug_statins", variable_name = "Statins", type = "Behavioural"
)


# stopdrug_3m_3mFU_MFN_hist

stopdrug_3m_3mFU_MFN_hist_sum <- univariate_analysis(
  cprd_dataset, drugs = c("DPP4", "GLP1", "SGLT2", "SU", "TZD"), outcome = "stopdrug_6m_6mFU", 
  variable = "stopdrug_3m_3mFU_MFN_hist", variable_name = "MFN discontinuation", type = "Behavioural"
)


# ethnicity

cprd_dataset <- cprd_dataset %>%
  mutate(ethnicity_5cat = factor(ethnicity_5cat, levels = c(0, 1, 2, 3, 4), labels = c("White", "Other", "Other" , "Other" , "Other")))

ethnicity_5cat_sum <- univariate_analysis(
  cprd_dataset, drugs = drugs, outcome = "stopdrug_6m_6mFU", 
  variable = "ethnicity_5cat", variable_name = "Ethnicity (ref White)", type = "Behavioural"
)


# smoking_cat

cprd_dataset <- cprd_dataset %>%
  mutate(smoking_cat = factor(smoking_cat, levels = c("Active smoker", "Ex-smoker", "Non-smoker"), labels = c("Active smoker", "Other", "Other")))

smoking_cat_sum <- univariate_analysis(
  cprd_dataset, drugs = drugs, outcome = "stopdrug_6m_6mFU", 
  variable = "smoking_cat", variable_name = "Smoking (ref Active Smoker)", type = "Behavioural"
)


# alcohol_cat

cprd_dataset <- cprd_dataset %>%
  mutate(alcohol_cat = factor(alcohol_cat, levels = c("Excess", "Harmful", "None", "Within limits"), labels = c("Excess/Harmful", "Excess/Harmful", "Within limits/None", "Within limits/None")))

alcohol_cat_sum <- univariate_analysis(
  cprd_dataset, drugs = drugs, outcome = "stopdrug_6m_6mFU", 
  variable = "alcohol_cat", variable_name = "Alcohol (ref Excess)", type = "Behavioural"
)


# predrug_cardio_event

predrug_cardio_event_sum <- univariate_analysis(
  cprd_dataset, drugs = drugs, outcome = "stopdrug_6m_6mFU", 
  variable = "predrug_cardio_event", variable_name = "Cardiovascular", type = "Comorbidity event"
)


# predrug_heart_event

predrug_heart_event_sum <- univariate_analysis(
  cprd_dataset, drugs = drugs, outcome = "stopdrug_6m_6mFU", 
  variable = "predrug_heart_event", variable_name = "Heart", type = "Comorbidity event"
)


# predrug_micro_event

predrug_micro_event_sum <- univariate_analysis(
  cprd_dataset, drugs = drugs, outcome = "stopdrug_6m_6mFU", 
  variable = "predrug_micro_event", variable_name = "Microvascular", type = "Comorbidity event"
)


# preckdstage

cprd_dataset <- cprd_dataset %>%
  mutate(preckdstage = factor(preckdstage, levels = c("stage_0", "stage_1", "stage_2", "stage_3a", "stage_3b", "stage_4", "stage_5"), labels = c("stage 0/1/2", "stage 0/1/2", "stage 0/1/2", "stage 3/4/5", "stage 3/4/5", "stage 3/4/5", "stage 3/4/5")))

preckdstage_sum <- univariate_analysis(
  cprd_dataset, drugs = drugs, outcome = "stopdrug_6m_6mFU", 
  variable = "preckdstage", variable_name = "CKD stage 3/4/5 (ref stage 0/1/2)", type = "Comorbidity event"
)


# predrug_cld

preckdstage_sum <- univariate_analysis(
  cprd_dataset, drugs = drugs, outcome = "stopdrug_6m_6mFU", 
  variable = "preckdstage", variable_name = "CKD stage 3/4/5 (ref stage 0/1/2)", type = "Comorbidity event"
)


# predrug_frailty_proxy

predrug_frailty_proxy_sum <- univariate_analysis(
  cprd_dataset, drugs = drugs, outcome = "stopdrug_6m_6mFU", 
  variable = "predrug_frailty_proxy", variable_name = "Frailty proxy", type = "Comorbidity event"
)



coefficients <- rbind(
  dstartdate_dm_dur_sum,
  dstartdate_age_sum,
  prehba1c_sum,
  preegfr_sum,
  prebmi_sum,
  prealt_sum,
  pretotalcholesterol_sum,
  gender_sum,
  predrug_bloodmed_sum,
  predrug_statins_sum,
  stopdrug_3m_3mFU_MFN_hist_sum,
  ethnicity_5cat_sum,
  smoking_cat_sum,
  alcohol_cat_sum,
  predrug_cardio_event_sum,
  predrug_heart_event_sum,
  predrug_micro_event_sum,
  preckdstage_sum,
  preckdstage_sum,
  predrug_frailty_proxy_sum
) %>%
  data.frame() %>%
  mutate(
    variable = factor(variable),
    drug = factor(drug, levels = rev(c("Pooled", drugs))),
    type = factor(type, levels = c("Clinical features & biomarkers", "Behavioural", "Comorbidity event")),
    coef = exp(as.numeric(coef)),
    LCI = exp(as.numeric(LCI)),
    UCI = exp(as.numeric(UCI))
  )


pdf("results/figures/univariate_analysis_6m.pdf", width = 12, height = 11)

wrap_plots(
  
  coefficients %>%
    filter(type == "Behavioural") %>%
    ggplot() + 
    geom_vline(aes(xintercept = 1), colour = "black") +
    geom_errorbar(aes(x = coef, xmin = LCI, xmax = UCI, y = variable, colour = drug),position = position_dodge(width = 1)) +
    geom_pointrange(aes(x = coef, xmin = LCI, xmax = UCI, y = variable, colour = drug), position = position_dodge(width = 1)) +
    xlab("Odds Ratio") +
    scale_x_log10(breaks = c(0.25, 0.5, 1, 2, 4), limits = c(0.4, 4.1)) +
    facet_grid(variable~., scales = "free_y") +
    scale_colour_manual(values = c("Pooled" = "black", "SGLT2" = "#E69F00", "GLP1" = "#56B4E9", "SU" = "#CC79A7", "DPP4" = "#0072B2", "TZD" = "#D55E00", "MFN" = "grey"), name = "Therapy", guide = guide_legend(reverse = TRUE)) +
    theme_bw() +
    theme(
      axis.title.y = element_blank(),
      legend.position = "bottom",
      strip.text = element_blank()
    ),
  
  
  coefficients %>%
    filter(type == "Clinical features & biomarkers") %>%
    ggplot() + 
    geom_vline(aes(xintercept = 1), colour = "black") +
    geom_errorbar(aes(x = coef, xmin = LCI, xmax = UCI, y = variable, colour = drug),position = position_dodge(width = 1)) +
    geom_pointrange(aes(x = coef, xmin = LCI, xmax = UCI, y = variable, colour = drug), position = position_dodge(width = 1)) +
    xlab("Odds Ratio") +
    scale_x_log10(breaks = c(0.25, 0.5, 1, 2, 4), limits = c(0.4, 4.1)) +
    facet_grid(variable~., scales = "free_y") +
    
    scale_colour_manual(values = c("Pooled" = "black", "SGLT2" = "#E69F00", "GLP1" = "#56B4E9", "SU" = "#CC79A7", "DPP4" = "#0072B2", "TZD" = "#D55E00", "MFN" = "grey"), name = "Therapy", guide = guide_legend(reverse = TRUE)) +
    theme_bw() +
    theme(
      axis.title.y = element_blank(),
      legend.position = "bottom",
      strip.text = element_blank()
    ),
  
  coefficients %>%
    filter(type == "Comorbidity event") %>%
    ggplot() + 
    geom_vline(aes(xintercept = 1), colour = "black") +
    geom_errorbar(aes(x = coef, xmin = LCI, xmax = UCI, y = variable, colour = drug),position = position_dodge(width = 1)) +
    geom_pointrange(aes(x = coef, xmin = LCI, xmax = UCI, y = variable, colour = drug), position = position_dodge(width = 1)) +
    xlab("Odds Ratio") +
    scale_x_log10(breaks = c(0.25, 0.5, 1, 2, 4), limits = c(0.4, 4.1)) +
    facet_grid(variable~., scales = "free_y") +
    
    scale_colour_manual(values = c("Pooled" = "black", "SGLT2" = "#E69F00", "GLP1" = "#56B4E9", "SU" = "#CC79A7", "DPP4" = "#0072B2", "TZD" = "#D55E00", "MFN" = "grey"), name = "Therapy", guide = guide_legend(reverse = TRUE)) +
    theme_bw() +
    theme(
      axis.title.y = element_blank(),
      legend.position = "bottom",
      strip.text = element_blank()
    ),
  
  ncol = 1, nrow = 3
  
) + 
  plot_layout(guides = "collect", height = c(6, 8, 6), axis_titles = "collect", axes = "collect") +
  plot_annotation(tag_levels = list(c("Behavioural", "Clinical features & biomarkers", "Comorbidity event"))) & 
  theme(
    legend.direction = "vertical",
    axis.text.y = element_text(size = 15),
    axis.text.x = element_text(size = 15),
    axis.title.x = element_text(size = 18)
  )

dev.off()




###############################################################################
###############################################################################
################# Univariate analysis - 12m discontinuation ####################
###############################################################################
###############################################################################


# load dataset
cprd_dataset <- cprd_dataset %>%
  drop_na(stopdrug_12m_6mFU)


drugs = c("MFN", "DPP4", "GLP1", "SGLT2", "SU", "TZD")


# dstartdate_dm_dur

cprd_dataset <- cprd_dataset %>%
  mutate_at(c("dstartdate_dm_dur"), ~(scale(.) %>% as.vector))

dstartdate_dm_dur_sum <- univariate_analysis(
  cprd_dataset, drugs = drugs, outcome = "stopdrug_12m_6mFU", 
  variable = "dstartdate_dm_dur", variable_name = "Diabetes duration (per SD)", type = "Clinical features & biomarkers"
)


# dstartdate_age

cprd_dataset <- cprd_dataset %>%
  mutate_at(c("dstartdate_age"), ~(scale(.) %>% as.vector))

dstartdate_age_sum <- univariate_analysis(
  cprd_dataset, drugs = drugs, outcome = "stopdrug_12m_6mFU", 
  variable = "dstartdate_age", variable_name = "Age (per SD)", type = "Clinical features & biomarkers"
)


# prehba1c

cprd_dataset <- cprd_dataset %>%
  mutate_at(c("prehba1c"), ~(scale(.) %>% as.vector))

prehba1c_sum <- univariate_analysis(
  cprd_dataset, drugs = drugs, outcome = "stopdrug_12m_6mFU", 
  variable = "prehba1c", variable_name = "HbA1c (per SD)", type = "Clinical features & biomarkers"
)


# preegfr

cprd <- cprd_dataset %>%
  mutate_at(c("preegfr"), ~(scale(.) %>% as.vector))

preegfr_sum <- univariate_analysis(
  cprd_dataset, drugs = drugs, outcome = "stopdrug_12m_6mFU", 
  variable = "preegfr", variable_name = "eGFR (per SD)", type = "Clinical features & biomarkers"
)


# prebmi

cprd_dataset <- cprd_dataset %>%
  mutate_at(c("prebmi"), ~(scale(.) %>% as.vector))

prebmi_sum <- univariate_analysis(
  cprd_dataset, drugs = drugs, outcome = "stopdrug_12m_6mFU", 
  variable = "prebmi", variable_name = "BMI (per SD)", type = "Clinical features & biomarkers"
)

# prealt

cprd <- cprd_dataset %>%
  mutate_at(c("prealt"), ~(scale(.) %>% as.vector))

prealt_sum <- univariate_analysis(
  cprd_dataset, drugs = drugs, outcome = "stopdrug_12m_6mFU", 
  variable = "prealt", variable_name = "ALT (per SD)", type = "Clinical features & biomarkers"
)


# pretotalcholesterol

cprd <- cprd_dataset %>%
  mutate_at(c("pretotalcholesterol"), ~(scale(.) %>% as.vector))

pretotalcholesterol_sum <- univariate_analysis(
  cprd_dataset, drugs = drugs, outcome = "stopdrug_12m_6mFU", 
  variable = "pretotalcholesterol", variable_name = "Total cholesterol (per SD)", type = "Clinical features & biomarkers"
)


# gender

gender_sum <- univariate_analysis(
  cprd_dataset, drugs = drugs, outcome = "stopdrug_12m_6mFU", 
  variable = "gender", variable_name = "Sex (ref Male)", type = "Clinical features & biomarkers"
)


# predrug_bloodmed

predrug_bloodmed_sum <- univariate_analysis(
  cprd_dataset, drugs = drugs, outcome = "stopdrug_12m_6mFU", 
  variable = "predrug_bloodmed", variable_name = "Blood pressure medication", type = "Behavioural"
)


# predrug_statins

predrug_statins_sum <- univariate_analysis(
  cprd_dataset, drugs = drugs, outcome = "stopdrug_12m_6mFU", 
  variable = "predrug_statins", variable_name = "Statins", type = "Behavioural"
)


# stopdrug_3m_3mFU_MFN_hist

stopdrug_3m_3mFU_MFN_hist_sum <- univariate_analysis(
  cprd_dataset, drugs = c("DPP4", "GLP1", "SGLT2", "SU", "TZD"), outcome = "stopdrug_12m_6mFU", 
  variable = "stopdrug_3m_3mFU_MFN_hist", variable_name = "MFN discontinuation", type = "Behavioural"
)


# ethnicity

cprd_dataset <- cprd_dataset %>%
  mutate(ethnicity_5cat = factor(ethnicity_5cat, levels = c(0, 1, 2, 3, 4), labels = c("White", "Other", "Other" , "Other" , "Other")))

ethnicity_5cat_sum <- univariate_analysis(
  cprd_dataset, drugs = drugs, outcome = "stopdrug_3m_6mFU", 
  variable = "ethnicity_5cat", variable_name = "Ethnicity (ref White)", type = "Behavioural"
)


# smoking_cat

cprd_dataset <- cprd_dataset %>%
  mutate(smoking_cat = factor(smoking_cat, levels = c("Active smoker", "Ex-smoker", "Non-smoker"), labels = c("Active smoker", "Other", "Other")))

smoking_cat_sum <- univariate_analysis(
  cprd_dataset, drugs = drugs, outcome = "stopdrug_12m_6mFU", 
  variable = "smoking_cat", variable_name = "Smoking (ref Active Smoker)", type = "Behavioural"
)


# alcohol_cat

cprd_dataset <- cprd_dataset %>%
  mutate(alcohol_cat = factor(alcohol_cat, levels = c("Excess", "Harmful", "None", "Within limits"), labels = c("Excess/Harmful", "Excess/Harmful", "Within limits/None", "Within limits/None")))

alcohol_cat_sum <- univariate_analysis(
  cprd_dataset, drugs = drugs, outcome = "stopdrug_12m_6mFU", 
  variable = "alcohol_cat", variable_name = "Alcohol (ref Excess)", type = "Behavioural"
)


# predrug_cardio_event

predrug_cardio_event_sum <- univariate_analysis(
  cprd_dataset, drugs = drugs, outcome = "stopdrug_12m_6mFU", 
  variable = "predrug_cardio_event", variable_name = "Cardiovascular", type = "Comorbidity event"
)


# predrug_heart_event

predrug_heart_event_sum <- univariate_analysis(
  cprd_dataset, drugs = drugs, outcome = "stopdrug_12m_6mFU", 
  variable = "predrug_heart_event", variable_name = "Heart", type = "Comorbidity event"
)


# predrug_micro_event

predrug_micro_event_sum <- univariate_analysis(
  cprd_dataset, drugs = drugs, outcome = "stopdrug_12m_6mFU", 
  variable = "predrug_micro_event", variable_name = "Microvascular", type = "Comorbidity event"
)


# preckdstage

cprd_dataset <- cprd_dataset %>%
  mutate(preckdstage = factor(preckdstage, levels = c("stage_0", "stage_1", "stage_2", "stage_3a", "stage_3b", "stage_4", "stage_5"), labels = c("stage 0/1/2", "stage 0/1/2", "stage 0/1/2", "stage 3/4/5", "stage 3/4/5", "stage 3/4/5", "stage 3/4/5")))

preckdstage_sum <- univariate_analysis(
  cprd_dataset, drugs = drugs, outcome = "stopdrug_12m_6mFU", 
  variable = "preckdstage", variable_name = "CKD stage 3/4/5 (ref stage 0/1/2)", type = "Comorbidity event"
)


# predrug_cld

preckdstage_sum <- univariate_analysis(
  cprd_dataset, drugs = drugs, outcome = "stopdrug_12m_6mFU", 
  variable = "preckdstage", variable_name = "CKD stage 3/4/5 (ref stage 0/1/2)", type = "Comorbidity event"
)


# predrug_frailty_proxy

predrug_frailty_proxy_sum <- univariate_analysis(
  cprd_dataset, drugs = drugs, outcome = "stopdrug_12m_6mFU", 
  variable = "predrug_frailty_proxy", variable_name = "Frailty proxy", type = "Comorbidity event"
)



coefficients <- rbind(
  dstartdate_dm_dur_sum,
  dstartdate_age_sum,
  prehba1c_sum,
  preegfr_sum,
  prebmi_sum,
  prealt_sum,
  pretotalcholesterol_sum,
  gender_sum,
  predrug_bloodmed_sum,
  predrug_statins_sum,
  stopdrug_3m_3mFU_MFN_hist_sum,
  ethnicity_5cat_sum,
  smoking_cat_sum,
  alcohol_cat_sum,
  predrug_cardio_event_sum,
  predrug_heart_event_sum,
  predrug_micro_event_sum,
  preckdstage_sum,
  preckdstage_sum,
  predrug_frailty_proxy_sum
) %>%
  data.frame() %>%
  mutate(
    variable = factor(variable),
    drug = factor(drug, levels = rev(c("Pooled", drugs))),
    type = factor(type, levels = c("Clinical features & biomarkers", "Behavioural", "Comorbidity event")),
    coef = exp(as.numeric(coef)),
    LCI = exp(as.numeric(LCI)),
    UCI = exp(as.numeric(UCI))
  )


pdf("results/figures/univariate_analysis_12m.pdf", width = 12, height = 11)

wrap_plots(
  
  coefficients %>%
    filter(type == "Behavioural") %>%
    ggplot() + 
    geom_vline(aes(xintercept = 1), colour = "black") +
    geom_errorbar(aes(x = coef, xmin = LCI, xmax = UCI, y = variable, colour = drug),position = position_dodge(width = 1)) +
    geom_pointrange(aes(x = coef, xmin = LCI, xmax = UCI, y = variable, colour = drug), position = position_dodge(width = 1)) +
    xlab("Odds Ratio") +
    scale_x_log10(breaks = c(0.25, 0.5, 1, 2, 4), limits = c(0.4, 4.1)) +
    facet_grid(variable~., scales = "free_y") +
    scale_colour_manual(values = c("Pooled" = "black", "SGLT2" = "#E69F00", "GLP1" = "#56B4E9", "SU" = "#CC79A7", "DPP4" = "#0072B2", "TZD" = "#D55E00", "MFN" = "grey"), name = "Therapy", guide = guide_legend(reverse = TRUE)) +
    theme_bw() +
    theme(
      axis.title.y = element_blank(),
      legend.position = "bottom",
      strip.text = element_blank()
    ),
  
  
  coefficients %>%
    filter(type == "Clinical features & biomarkers") %>%
    ggplot() + 
    geom_vline(aes(xintercept = 1), colour = "black") +
    geom_errorbar(aes(x = coef, xmin = LCI, xmax = UCI, y = variable, colour = drug),position = position_dodge(width = 1)) +
    geom_pointrange(aes(x = coef, xmin = LCI, xmax = UCI, y = variable, colour = drug), position = position_dodge(width = 1)) +
    xlab("Odds Ratio") +
    scale_x_log10(breaks = c(0.25, 0.5, 1, 2, 4), limits = c(0.4, 4.1)) +
    facet_grid(variable~., scales = "free_y") +
    
    scale_colour_manual(values = c("Pooled" = "black", "SGLT2" = "#E69F00", "GLP1" = "#56B4E9", "SU" = "#CC79A7", "DPP4" = "#0072B2", "TZD" = "#D55E00", "MFN" = "grey"), name = "Therapy", guide = guide_legend(reverse = TRUE)) +
    theme_bw() +
    theme(
      axis.title.y = element_blank(),
      legend.position = "bottom",
      strip.text = element_blank()
    ),
  
  coefficients %>%
    filter(type == "Comorbidity event") %>%
    ggplot() + 
    geom_vline(aes(xintercept = 1), colour = "black") +
    geom_errorbar(aes(x = coef, xmin = LCI, xmax = UCI, y = variable, colour = drug),position = position_dodge(width = 1)) +
    geom_pointrange(aes(x = coef, xmin = LCI, xmax = UCI, y = variable, colour = drug), position = position_dodge(width = 1)) +
    xlab("Odds Ratio") +
    scale_x_log10(breaks = c(0.25, 0.5, 1, 2, 4), limits = c(0.4, 4.1)) +
    facet_grid(variable~., scales = "free_y") +
    
    scale_colour_manual(values = c("Pooled" = "black", "SGLT2" = "#E69F00", "GLP1" = "#56B4E9", "SU" = "#CC79A7", "DPP4" = "#0072B2", "TZD" = "#D55E00", "MFN" = "grey"), name = "Therapy", guide = guide_legend(reverse = TRUE)) +
    theme_bw() +
    theme(
      axis.title.y = element_blank(),
      legend.position = "bottom",
      strip.text = element_blank()
    ),
  
  ncol = 1, nrow = 3
  
) + 
  plot_layout(guides = "collect", height = c(6, 8, 6), axis_titles = "collect", axes = "collect") +
  plot_annotation(tag_levels = list(c("Behavioural", "Clinical features & biomarkers", "Comorbidity event"))) & 
  theme(
    legend.direction = "vertical",
    axis.text.y = element_text(size = 15),
    axis.text.x = element_text(size = 15),
    axis.title.x = element_text(size = 18)
  )

dev.off()
