####################
## Description: 
##  - In this file we:
##    - Fit a BART model all the variables into one model.
#################### 

## increase memory usage to 100GB of RAM (needs to be run before library(bartMachine))
options(java.parameters = "-Xmx100g")

# load working directory
setwd("Samples/T2D_Discontinuation")

# load functions
source("code/00.set_up_data_and_functions.R")

# load libraries
library(tidyverse)
library(pROC)
library(predtools)
library(bartMachine)
library(BART)
library(patchwork)



###############################################################################
###############################################################################
############################### Read Data In ##################################
###############################################################################
###############################################################################

# load dataset
cprd_dataset.3m <- set_up_data(
  raw_data = "20240814_t2d_1stinstance",
  diagnosis = FALSE,
  therapies = c("GLP1", "DPP4", "SGLT2", "TZD", "SU"),
  dataset = "3m.disc.dataset",
  follow_up = "6-months",
  full_prescribing_history = TRUE
) %>%
  drop_na(-prehdl, -stopdrug_6m_6mFU, -stopdrug_12m_6mFU)

# load model
bartmachine_full_model_3m <- readRDS("results/Models/bartmachine/bartmachine_full_model_3m.rds")



# Predict percentage for 3m
## SGLT2
if (class(try(
  
  predictions_SGLT2_3m <- readRDS("results/Models/bartmachine/predictions_SGLT2_3m.rds")
  
  , silent = TRUE)) == "try-error") {
  
  predictions_SGLT2_3m <- predict(bartmachine_full_model_3m,
                                  cprd_dataset.3m %>%
                                    select(drugclass, 
                                           dstartdate_age, 
                                           gender, 
                                           imd2015_10, 
                                           prebmi, 
                                           dstartdate_dm_dur, 
                                           prehba1c, 
                                           drugline, 
                                           predrug_frailty_proxy,
                                           ethnicity_5cat,
                                           numdrugs,
                                           predrug_bloodmed,
                                           smoking_cat,
                                           predrug_statins,
                                           preegfr,
                                           stopdrug_3m_3mFU_MFN_hist) %>%
                                    mutate(drugclass = "SGLT2"))
  
  saveRDS(predictions_SGLT2_3m, "results/Models/bartmachine/predictions_SGLT2_3m.rds")
  
}


## GLP1
if (class(try(
  
  predictions_GLP1_3m <- readRDS("results/Models/bartmachine/predictions_GLP1_3m.rds")
  
  , silent = TRUE)) == "try-error") {
  
  predictions_GLP1_3m <- predict(bartmachine_full_model_3m,
                                  cprd_dataset.3m %>%
                                    select(drugclass, 
                                           dstartdate_age, 
                                           gender, 
                                           imd2015_10, 
                                           prebmi, 
                                           dstartdate_dm_dur, 
                                           prehba1c, 
                                           drugline, 
                                           predrug_frailty_proxy,
                                           ethnicity_5cat,
                                           numdrugs,
                                           predrug_bloodmed,
                                           smoking_cat,
                                           predrug_statins,
                                           preegfr,
                                           stopdrug_3m_3mFU_MFN_hist) %>%
                                    mutate(drugclass = "GLP1"))
  
  saveRDS(predictions_GLP1_3m, "results/Models/bartmachine/predictions_GLP1_3m.rds")
  
}


## DPP4
if (class(try(
  
  predictions_DPP4_3m <- readRDS("results/Models/bartmachine/predictions_DPP4_3m.rds")
  
  , silent = TRUE)) == "try-error") {
  
  predictions_DPP4_3m <- predict(bartmachine_full_model_3m,
                                 cprd_dataset.3m %>%
                                   select(drugclass, 
                                          dstartdate_age, 
                                          gender, 
                                          imd2015_10, 
                                          prebmi, 
                                          dstartdate_dm_dur, 
                                          prehba1c, 
                                          drugline, 
                                          predrug_frailty_proxy,
                                          ethnicity_5cat,
                                          numdrugs,
                                          predrug_bloodmed,
                                          smoking_cat,
                                          predrug_statins,
                                          preegfr,
                                          stopdrug_3m_3mFU_MFN_hist) %>%
                                   mutate(drugclass = "DPP4"))
  
  saveRDS(predictions_DPP4_3m, "results/Models/bartmachine/predictions_DPP4_3m.rds")
  
}


## TZD
if (class(try(
  
  predictions_TZD_3m <- readRDS("results/Models/bartmachine/predictions_TZD_3m.rds")
  
  , silent = TRUE)) == "try-error") {
  
  predictions_TZD_3m <- predict(bartmachine_full_model_3m,
                                 cprd_dataset.3m %>%
                                   select(drugclass, 
                                          dstartdate_age, 
                                          gender, 
                                          imd2015_10, 
                                          prebmi, 
                                          dstartdate_dm_dur, 
                                          prehba1c, 
                                          drugline, 
                                          predrug_frailty_proxy,
                                          ethnicity_5cat,
                                          numdrugs,
                                          predrug_bloodmed,
                                          smoking_cat,
                                          predrug_statins,
                                          preegfr,
                                          stopdrug_3m_3mFU_MFN_hist) %>%
                                   mutate(drugclass = "TZD"))
  
  saveRDS(predictions_TZD_3m, "results/Models/bartmachine/predictions_TZD_3m.rds")
  
}


## SU
if (class(try(
  
  predictions_SU_3m <- readRDS("results/Models/bartmachine/predictions_SU_3m.rds")
  
  , silent = TRUE)) == "try-error") {
  
  predictions_SU_3m <- predict(bartmachine_full_model_3m,
                                cprd_dataset.3m %>%
                                  select(drugclass, 
                                         dstartdate_age, 
                                         gender, 
                                         imd2015_10, 
                                         prebmi, 
                                         dstartdate_dm_dur, 
                                         prehba1c, 
                                         drugline, 
                                         predrug_frailty_proxy,
                                         ethnicity_5cat,
                                         numdrugs,
                                         predrug_bloodmed,
                                         smoking_cat,
                                         predrug_statins,
                                         preegfr,
                                         stopdrug_3m_3mFU_MFN_hist) %>%
                                  mutate(drugclass = "SU"))
  
  saveRDS(predictions_SU_3m, "results/Models/bartmachine/predictions_SU_3m.rds")
  
}


###############################################################################
breakpoints <- c(-0.1, -0.05, -0.025, 0, 0.025, 0.05, 0.1)

interim_dataset <- cprd_dataset.3m %>%
  as.data.frame() %>%
  mutate(
    pred.SGLT2 = predictions_SGLT2_3m,
    pred.GLP1 = predictions_GLP1_3m,
    pred.DPP4 = predictions_DPP4_3m,
    pred.SU = predictions_SU_3m,
    pred.TZD = predictions_TZD_3m
    )

# GLP1/DPP4
interim_dataset %>% filter(drugclass %in% c("GLP1", "DPP4")) %>%
  mutate(benefit = pred.GLP1 - pred.DPP4) %>% select(benefit) %>% unlist() %>% max(na.rm = TRUE)
# GLP1/SU
interim_dataset %>% filter(drugclass %in% c("GLP1", "SU")) %>%
  mutate(benefit = pred.GLP1 - pred.SU) %>% select(benefit) %>% unlist() %>% max(na.rm = TRUE)
# GLP1/TZD
interim_dataset %>% filter(drugclass %in% c("GLP1", "TZD")) %>%
  mutate(benefit = pred.GLP1 - pred.TZD) %>% select(benefit) %>% unlist() %>% max(na.rm = TRUE)
# SGLT2/DPP4
interim_dataset %>% filter(drugclass %in% c("SGLT2", "DPP4")) %>%
  mutate(benefit = pred.SGLT2 - pred.DPP4) %>% select(benefit) %>% unlist() %>% max(na.rm = TRUE)
# SGLT2/GLP1
interim_dataset %>% filter(drugclass %in% c("SGLT2", "GLP1")) %>%
  mutate(benefit = pred.SGLT2 - pred.GLP1) %>% select(benefit) %>% unlist() %>% max(na.rm = TRUE)
# SGLT2/SU
interim_dataset %>% filter(drugclass %in% c("SGLT2", "SU")) %>%
  mutate(benefit = pred.SGLT2 - pred.SU) %>% select(benefit) %>% unlist() %>% max(na.rm = TRUE)
# SGLT2/TZD
interim_dataset %>% filter(drugclass %in% c("SGLT2", "TZD")) %>%
  mutate(benefit = pred.SGLT2 - pred.TZD) %>% select(benefit) %>% unlist() %>% max(na.rm = TRUE)
# SU/DPP4
interim_dataset %>% filter(drugclass %in% c("SU", "DPP4")) %>%
  mutate(benefit = pred.SU - pred.DPP4) %>% select(benefit) %>% unlist() %>% max(na.rm = TRUE)
# TZD/DPP4
interim_dataset %>% filter(drugclass %in% c("TZD", "DPP4")) %>%
  mutate(benefit = pred.TZD - pred.DPP4) %>% select(benefit) %>% unlist() %>% max(na.rm = TRUE)
# TZD/SU
interim_dataset %>% filter(drugclass %in% c("TZD", "SU")) %>%
  mutate(benefit = pred.TZD - pred.SU) %>% select(benefit) %>% unlist() %>% max(na.rm = TRUE)




## Histogram of heterogeneities
#GLP1/DPP4
calc_CATE <- function(data, drugs, title) {
  
  potential_drugs <- c("GLP1", "DPP4", "SGLT2", "TZD", "SU"); drugs_names <- c("GLP-1RA", "DPP4i", "SGLT2i", "TZD", "SU")
  potential_colours <- c("#56B4E9", "#0072B2", "#E69F00", "#D55E00", "#CC79A7")
  positions <- c(which(potential_drugs %in% drugs[1]), which(potential_drugs %in% drugs[2]))
  current_drugs_names <- drugs_names[positions]
  current_drugs_colours <- potential_colours[positions]
  drug_1_name = paste0("pred.", drugs[1]); drug_2_name = paste0("pred.", drugs[2]); 
  
  interim_dataset <- data %>%
    filter(drugclass %in% drugs) %>%
    rename(
      "pred_drug_A" = drug_1_name, "pred_drug_B" = drug_2_name
    ) %>%
    mutate(benefit = pred_drug_A - pred_drug_B) %>%
    mutate(above = ifelse(benefit < 0, paste("Favours", current_drugs_names[1]), paste("Favours", current_drugs_names[2]))) %>%
    mutate(above = factor(above, levels = c(paste("Favours", current_drugs_names[1]), paste("Favours", current_drugs_names[2])))) %>%
    select(benefit, above) %>%
    mutate(title = paste(current_drugs_names[1], "-", current_drugs_names[2])) %>%
    as.data.frame()
  
  return(interim_dataset)
  
}

pdf("results/figures/08.CATE_histogram_drugs_3m.pdf", width = 13, height = 6)
calc_CATE(interim_dataset, drugs = c("SGLT2", "DPP4")) %>%
  rbind(
    calc_CATE(interim_dataset, drugs = c("SGLT2", "GLP1")),
    calc_CATE(interim_dataset, drugs = c("SGLT2", "SU")),
    calc_CATE(interim_dataset, drugs = c("SGLT2", "TZD")),
    calc_CATE(interim_dataset, drugs = c("SU", "DPP4")),
    calc_CATE(interim_dataset, drugs = c("GLP1", "DPP4")),
    calc_CATE(interim_dataset, drugs = c("GLP1", "SU")),
    calc_CATE(interim_dataset, drugs = c("GLP1", "TZD")),
    calc_CATE(interim_dataset, drugs = c("TZD", "SU")),
    calc_CATE(interim_dataset, drugs = c("TZD", "DPP4"))
  ) %>%
  mutate(title = factor(title, levels = c("SGLT2i - DPP4i", "SGLT2i - GLP-1RA", "SGLT2i - SU", "SGLT2i - TZD", "SU - DPP4i", "GLP-1RA - DPP4i", "GLP-1RA - SU", "GLP-1RA - TZD", "TZD - SU", "TZD - DPP4i"))) %>%
  ggplot(aes(x = benefit, fill = above)) +
  geom_histogram(position = "identity", alpha = 0.5, color = "black", breaks = seq(-0.25, 0.4, 0.025)) +
  geom_vline(aes(xintercept = 0), linetype = "dashed") +
  labs(x = "Predicted risk of discontinuation difference (%)", y = "Number of people") +
  scale_fill_manual(
    values = c("Favours SGLT2i" = "#E69F00", "Favours GLP-1RA" = "#56B4E9", "Favours SU" = "#CC79A7", "Favours DPP4i" = "#0072B2", "Favours TZD" = "#D55E00"),
    breaks = c("Favours GLP-1RA", "Favours DPP4i", "Favours SGLT2i", "Favours TZD", "Favours SU")) +
  scale_x_continuous(labels = scales::percent, breaks = seq(-0.3, 0.5, 0.1)) +
  theme_bw() +
  facet_wrap(~title, nrow = 2) +
  theme(legend.title = element_blank(),
        legend.direction = "horizontal",
        legend.position = "bottom",
        legend.box = "horizontal")

dev.off()




# vars <- c(
#   # Outcome
#   "stopdrug_3m_6mFU",
#   "stopdrug_6m_6mFU",
#   "stopdrug_12m_6mFU",
#   "drugclass",
#   # Extra info
#   "dstartdate_age", "dstartdate_dm_dur", "gender", "smoking_cat", 
#   "ethnicity_5cat", "imd2015_10",
#   # Biomarkers
#   "prebmi", "prehba1c", "preegfr", 
#   "predrug_statins", "predrug_bloodmed", "predrug_frailty_proxy", "stopdrug_3m_3mFU_MFN_hist",
#   "numdrugs", "drugline"
#   # Comorbidities
# )
# 
# vars_cat <- c(
#   # Outcome
#   "stopdrug_3m_6mFU",
#   "stopdrug_6m_6mFU",
#   "stopdrug_12m_6mFU",
#   "drugclass",
#   # Extra info
#   "drugline", "numdrugs", "smoking_cat", "imd2015_10",
#   "predrug_statins", "stopdrug_3m_3mFU_MFN_hist", "ethnicity_5cat", "gender", "predrug_bloodmed",
#   # Comorbidities
#   ## Frailty proxy
#   "predrug_frailty_proxy"
# )
# 
# 
# # SGLT2/DPP4
# table_characteristics_drugs <- CreateTableOne(
#   vars = vars,
#   factorVars = vars_cat,
#   includeNA = TRUE,
#   strata = c("group.benefit"),
#   data = interim_dataset %>% filter(drugclass %in% c("SGLT2", "DPP4")) %>%
#     mutate(benefit = pred.SGLT2 - pred.DPP4) %>% 
#     mutate(group.benefit = ntile(benefit, 10)) %>%
#     mutate(group.benefit = ifelse(group.benefit < 2, "Most disc on DPP4", ifelse(group.benefit >= 2 & group.benefit < 10, "Average", "Most disc on SGLT2"))) %>%
#     mutate(group.benefit = factor(group.benefit, levels = c("Most disc on DPP4", "Average" ,"Most disc on SGLT2"))),
#   test = FALSE
# )
# 
# table_characteristics_drugs_print <- print(table_characteristics_drugs, exact = "stage", quote = FALSE, noSpaces = TRUE, printToggle = FALSE, contDigits = 1)
# 
# write.csv(table_characteristics_drugs_print, file = "results/tables/08.01.heterogeneity_SGLT2_DPP4.csv")
# 
# # GLP1/TZD
# table_characteristics_drugs <- CreateTableOne(
#   vars = vars,
#   factorVars = vars_cat,
#   includeNA = TRUE,
#   strata = c("group.benefit"),
#   data = interim_dataset %>% filter(drugclass %in% c("GLP1", "TZD")) %>%
#     mutate(benefit = pred.GLP1 - pred.TZD) %>% 
#     mutate(group.benefit = ntile(benefit, 10)) %>%
#     mutate(group.benefit = ifelse(group.benefit < 2, "Most disc on TZD", ifelse(group.benefit >= 2 & group.benefit < 10, "Average", "Most disc on GLP1"))) %>%
#     mutate(group.benefit = factor(group.benefit, levels = c("Most disc on TZD", "Average" ,"Most disc on GLP1"))),
#   test = FALSE
# )
# 
# table_characteristics_drugs_print <- print(table_characteristics_drugs, exact = "stage", quote = FALSE, noSpaces = TRUE, printToggle = FALSE, contDigits = 1)
# 
# write.csv(table_characteristics_drugs_print, file = "results/tables/08.01.heterogeneity_GLP1_TZD.csv")

  

# ATE by deciles
ATE.deciles.adj <- calc_ATE(interim_dataset, ntiles = 10, drugs = c("SGLT2", "GLP1", "TZD", "SU", "DPP4"), pred.variable = "", n_bootstrap = 50,
                               breakdown = c("dstartdate_age", 
                                             "gender", 
                                             "imd2015_10", 
                                             "prebmi", 
                                             "dstartdate_dm_dur", 
                                             "prehba1c", 
                                             "drugline", 
                                             "predrug_frailty_proxy",
                                             "ethnicity_5cat",
                                             "numdrugs",
                                             "predrug_bloodmed",
                                             "smoking_cat",
                                             "predrug_statins",
                                             "preegfr",
                                             "stopdrug_3m_3mFU_MFN_hist"))



pdf("results/figures/08.drug_vs_drug_calibration_3m.pdf", width = 13, height = 6)

ATE.deciles.adj %>%
  group_by(Type) %>%
  mutate(
    total = sum(N),
    drug_1 = ifelse(drug_1 == "SGLT2", "SGLT2i", ifelse(drug_1 == "GLP1", "GLP-1RA", ifelse(drug_1 == "DPP4", "DPP4i", drug_1))),
    drug_2 = ifelse(drug_2 == "SGLT2", "SGLT2i", ifelse(drug_2 == "GLP1", "GLP-1RA", ifelse(drug_2 == "DPP4", "DPP4i", drug_2))),
    title = paste0(drug_1, " - ", drug_2, "\nn = ", format(total, big.mark = ","))
  ) %>%
  ungroup() %>%
  mutate(
    title = factor(title, levels = c("SGLT2i - DPP4i\nn = 109,101", "SGLT2i - GLP-1RA\nn = 60,641", "SGLT2i - SU\nn = 79,801", "SGLT2i - TZD\nn = 50,204", "SU - DPP4i\nn = 97,590", "GLP-1RA - DPP4i\nn = 78,430", "GLP-1RA - SU\nn = 49,130", "GLP-1RA - TZD\nn = 19,533", "TZD - SU\nn = 38,693", "TZD - DPP4i\nn = 67,993"))
  ) %>%
  ggplot(aes(x = diff.pred, y = obs, ymin = lci, ymax = uci)) +
  geom_vline(aes(xintercept = 0), linetype = "dashed", colour = "black") +
  geom_hline(aes(yintercept = 0), linetype = "dashed", colour = "black") +
  geom_abline(aes(intercept = 0, slope = 1), colour = "red") +
  geom_errorbar() +
  geom_point() +
  labs(x = "Predicted risk of discontinuation difference (%)", y = "Observed discontinuation difference (%)") +
  # ggtitle("T2D 3-month discontinuation calibration (adjusted, 50 bootstraps)") +
  scale_x_continuous(labels = scales::percent, breaks = seq(-0.3, 0.3, 0.05)) +
  scale_y_continuous(labels = scales::percent, breaks = seq(-0.3, 0.3, 0.05)) +
  theme_bw() +
  facet_wrap(~title, nrow = 2) +
  theme(
    panel.spacing = unit(0.3, "cm", data = NULL)
  )
  
dev.off()


# ATE not adjusted
ATE.overall.no_adj <- calc_ATE(interim_dataset, break_points = breakpoints, drugs = c("SGLT2", "GLP1", "TZD", "SU", "DPP4"), pred.variable = "", n_bootstrap = 10) %>%
  mutate(label = ifelse(group == "(-1,-0.1]", paste0("Benefit >10% (n=", N, ", events = ", events, ")"), ifelse(
    group == "(-0.1,-0.05]", paste0("Benefit 5-10% (n=", N, ", events = ", events, ")"), ifelse(
      group == "(-0.05,-0.025]", paste0("Benefit 2.5-5% (n=", N, ", events = ", events, ")"), ifelse(
        group == "(-0.025,0]", paste0("Benefit 0-2.5% (n=", N, ", events = ", events, ")"), ifelse(
          group == "(0,0.025]", paste0("Benefit 0-2.5% (n=", N, ", events = ", events, ")"), ifelse(
            group == "(0.025,0.05]", paste0("Benefit 2.5-5% (n=", N, ", events = ", events, ")"), ifelse(
              group == "(0.05,0.1]", paste0("Benefit 5-10% (n=", N, ", events = ", events, ")"), ifelse(
                group == "(0.1,1]", paste0("Benefit >10% (n=", N, ", events = ", events, ")"), group
              )
            )
          )
        )
      )
    )
  )))

# Predicted benefit
absolute.overall.no_adj <- calc_predicted_risk(interim_dataset, break_points = breakpoints, drugs = c("SGLT2", "GLP1", "TZD", "SU", "DPP4"), pred.variable = "") %>%
  mutate(label = ifelse(group == "(-1,-0.1]", paste0("Benefit >10% (n=", N, ", events = ", events, ")"), ifelse(
    group == "(-0.1,-0.05]", paste0("Benefit 5-10% (n=", N, ", events = ", events, ")"), ifelse(
      group == "(-0.05,-0.025]", paste0("Benefit 2.5-5% (n=", N, ", events = ", events, ")"), ifelse(
        group == "(-0.025,0]", paste0("Benefit 0-2.5% (n=", N, ", events = ", events, ")"), ifelse(
          group == "(0,0.025]", paste0("Benefit 0-2.5% (n=", N, ", events = ", events, ")"), ifelse(
            group == "(0.025,0.05]", paste0("Benefit 2.5-5% (n=", N, ", events = ", events, ")"), ifelse(
              group == "(0.05,0.1]", paste0("Benefit 5-10% (n=", N, ", events = ", events, ")"), ifelse(
                group == "(0.1,1]", paste0("Benefit >10% (n=", N, ", events = ", events, ")"), group
              )
            )
          )
        )
      )
    )
  )))


# Odds ratios 
odds_ratios.overall.no_adj <- calc_odds_ratios(interim_dataset, drugs = c("SGLT2", "GLP1", "TZD", "SU", "DPP4"), pred.variable = "", break_points = breakpoints) %>%
  mutate(label = ifelse(group == "(-1,-0.1]", paste0("Benefit >10% (n=", N, ", events = ", events, ")"), ifelse(
    group == "(-0.1,-0.05]", paste0("Benefit 5-10% (n=", N, ", events = ", events, ")"), ifelse(
      group == "(-0.05,-0.025]", paste0("Benefit 2.5-5% (n=", N, ", events = ", events, ")"), ifelse(
        group == "(-0.025,0]", paste0("Benefit 0-2.5% (n=", N, ", events = ", events, ")"), ifelse(
          group == "(0,0.025]", paste0("Benefit 0-2.5% (n=", N, ", events = ", events, ")"), ifelse(
            group == "(0.025,0.05]", paste0("Benefit 2.5-5% (n=", N, ", events = ", events, ")"), ifelse(
              group == "(0.05,0.1]", paste0("Benefit 5-10% (n=", N, ", events = ", events, ")"), ifelse(
                group == "(0.1,1]", paste0("Benefit >10% (n=", N, ", events = ", events, ")"), group
              )
            )
          )
        )
      )
    )
  )))






pdf("results/figures/08.absolute_relative_development_3m.pdf", width = 20, height = 7)


##: Per drug combination: SGLT2vsDPP4
drug1 = "SGLT2"; drug2 = "DPP4"

title <- ATE.overall.no_adj[1,]
title[,] <- NA

interim_1 <- title %>% as.data.frame() %>%mutate(group = paste0("Predicted benefit on ", drug1, ":"), label = paste0("Predicted benefit on ", drug1, ":")) %>%
  rbind(
    ATE.overall.no_adj %>%
      filter(drug_1 == drug1 & drug_2 == drug2 & diff.pred < 0),
    title %>% as.data.frame() %>% mutate(group = paste0("Predicted benefit on ", drug2, ":"), label = paste0("Predicted benefit on ", drug2, ":")),
    ATE.overall.no_adj %>%
      filter(drug_1 == drug1 & drug_2 == drug2 & diff.pred > 0)
  )

title <- odds_ratios.overall.no_adj[1,]
title[,] <- NA


interim_2 <- title %>% as.data.frame() %>%mutate(group = paste0("Predicted benefit on ", drug1, ":"), label = paste0("Predicted benefit on ", drug1, ":")) %>%
  rbind(
    odds_ratios.overall.no_adj %>%
      filter(drug_1 == drug1 & drug_2 == drug2 & obs > 0),
    title %>% as.data.frame() %>% mutate(group = paste0("Predicted benefit on ", drug2, ":"), label = paste0("Predicted benefit on ", drug2, ":")),
    odds_ratios.overall.no_adj %>%
      filter(drug_1 == drug1 & drug_2 == drug2 & obs < 0)
  )

title <- absolute.overall.no_adj[1,]
title[,] <- NA

interim_3 <- title %>% as.data.frame() %>% mutate(group = paste0("Predicted benefit on ", drug1, ":"), label = paste0("Predicted benefit on ", drug1, ":")) %>%
  mutate(drug_type = drug1) %>%
  rbind(
    absolute.overall.no_adj %>%
      filter(drug_1 == drug1 & drug_2 == drug2 & diff.pred < 0),
    title %>% as.data.frame() %>% mutate(group = paste0("Predicted benefit on ", drug2, ":"), label = paste0("Predicted benefit on ", drug2, ":")) %>%
      mutate(drug_type = drug2),
    absolute.overall.no_adj %>%
      filter(drug_1 == drug1 & drug_2 == drug2 & diff.pred > 0)
  )


patchwork::wrap_plots(
  
  interim_1 %>%
    ggplot(aes(y = group, x = obs, xmin = lci, xmax = uci)) +
    geom_vline(aes(xintercept = 0), colour = "grey") +
    geom_point() +
    geom_errorbar() + 
    scale_y_discrete(
      limits= rev(interim_1$group),
      labels = rev(interim_1$label %>% unlist())
    ) +
    scale_x_continuous("Predicted Treatment Discontinuation Difference", labels = scales::percent, breaks = seq(-1, 1, 0.02)) +
    coord_cartesian(xlim = c(-0.15, 0.1)) +
    ggtitle(paste0("Absolute benefit (below 0 = benefit on ", drug1,")")) +
    theme_bw(),
  
  interim_3 %>%
    ggplot(aes(y = group, x = obs, xmin = lci, xmax = uci, colour = drug_type)) +
    geom_point(position = position_dodge(width = 0.5)) +
    geom_errorbar(position = position_dodge(width = 0.5)) +
    scale_y_discrete(
      limits= rev(interim_1$group),
      labels = rev(interim_1$label %>% unlist())
    ) +
    scale_x_continuous("Predicted Discontinuation", labels = scales::percent) +
    scale_colour_manual("Therapy", values = c("GLP1" = "#56B4E9", "SGLT2" = "#E69F00", "SU" = "#CC79A7", "DPP4" = "#0072B2", "TZD" = "#D55E00")) +
    coord_cartesian(xlim = c(0.05, 0.3)) +
    ggtitle("Absolute risk") +
    theme_bw() +
    theme(
      legend.position = "bottom"
    ),
  
  
  interim_2 %>%
    mutate(
      obs = exp(obs),
      lci = exp(lci),
      uci = exp(uci)
    ) %>%
    ggplot() + 
    geom_vline(aes(xintercept = 1), colour = "grey") +
    geom_point(aes(y = group, x = obs)) + 
    geom_errorbar(aes(y = group, xmin = lci, xmax = uci, x = obs)) + 
    scale_y_discrete(
      limits= rev(interim_1$group),
      labels = rev(interim_1$label %>% unlist())
    ) +
    xlab("Odds Ratio") +
    scale_x_log10(breaks = c(0.25, 0.5, 1, 2, 4)) +
    coord_cartesian(xlim = c(0.2, 4.1)) +
    ggtitle(paste0("Relative benefit (below 1 = benefit on ", drug1, ")")) +
    theme_bw()
  
) + patchwork::plot_layout(axes = "collect") +
  patchwork::plot_annotation(
    title = paste0(drug1, " vs ", drug2)
  ) &
  theme(
    axis.title.y = element_blank(),
    axis.title.x = element_text(size = 13),
    axis.text = element_text(size = 13)
  )


##: Per drug combination: SGLT2vsGLP1
drug1 = "SGLT2"; drug2 = "GLP1"

title <- ATE.overall.no_adj[1,]
title[,] <- NA

interim_1 <- title %>% as.data.frame() %>%mutate(group = paste0("Predicted benefit on ", drug1, ":"), label = paste0("Predicted benefit on ", drug1, ":")) %>%
  rbind(
    ATE.overall.no_adj %>%
      filter(drug_1 == drug1 & drug_2 == drug2 & diff.pred < 0),
    title %>% as.data.frame() %>% mutate(group = paste0("Predicted benefit on ", drug2, ":"), label = paste0("Predicted benefit on ", drug2, ":")),
    ATE.overall.no_adj %>%
      filter(drug_1 == drug1 & drug_2 == drug2 & diff.pred > 0)
  )

title <- odds_ratios.overall.no_adj[1,]
title[,] <- NA


interim_2 <- title %>% as.data.frame() %>%mutate(group = paste0("Predicted benefit on ", drug1, ":"), label = paste0("Predicted benefit on ", drug1, ":")) %>%
  rbind(
    odds_ratios.overall.no_adj %>%
      filter(drug_1 == drug1 & drug_2 == drug2 & obs > 0),
    title %>% as.data.frame() %>% mutate(group = paste0("Predicted benefit on ", drug2, ":"), label = paste0("Predicted benefit on ", drug2, ":")),
    odds_ratios.overall.no_adj %>%
      filter(drug_1 == drug1 & drug_2 == drug2 & obs < 0)
  )

title <- absolute.overall.no_adj[1,]
title[,] <- NA

interim_3 <- title %>% as.data.frame() %>% mutate(group = paste0("Predicted benefit on ", drug1, ":"), label = paste0("Predicted benefit on ", drug1, ":")) %>%
  mutate(drug_type = drug1) %>%
  rbind(
    absolute.overall.no_adj %>%
      filter(drug_1 == drug1 & drug_2 == drug2 & diff.pred < 0),
    title %>% as.data.frame() %>% mutate(group = paste0("Predicted benefit on ", drug2, ":"), label = paste0("Predicted benefit on ", drug2, ":")) %>%
      mutate(drug_type = drug2),
    absolute.overall.no_adj %>%
      filter(drug_1 == drug1 & drug_2 == drug2 & diff.pred > 0)
  )


patchwork::wrap_plots(
  
  interim_1 %>%
    ggplot(aes(y = group, x = obs, xmin = lci, xmax = uci)) +
    geom_vline(aes(xintercept = 0), colour = "grey") +
    geom_point() +
    geom_errorbar() + 
    scale_y_discrete(
      limits= rev(interim_1$group),
      labels = rev(interim_1$label %>% unlist())
    ) +
    scale_x_continuous("Predicted Treatment Discontinuation Difference", labels = scales::percent, breaks = seq(-1, 1, 0.02)) +
    coord_cartesian(xlim = c(-0.15, 0.1)) +
    ggtitle(paste0("Absolute benefit (below 0 = benefit on ", drug1,")")) +
    theme_bw(),
  
  interim_3 %>%
    ggplot(aes(y = group, x = obs, xmin = lci, xmax = uci, colour = drug_type)) +
    geom_point(position = position_dodge(width = 0.5)) +
    geom_errorbar(position = position_dodge(width = 0.5)) +
    scale_y_discrete(
      limits= rev(interim_1$group),
      labels = rev(interim_1$label %>% unlist())
    ) +
    scale_x_continuous("Predicted Discontinuation", labels = scales::percent) +
    scale_colour_manual("Therapy", values = c("GLP1" = "#56B4E9", "SGLT2" = "#E69F00", "SU" = "#CC79A7", "DPP4" = "#0072B2", "TZD" = "#D55E00")) +
    coord_cartesian(xlim = c(0.05, 0.3)) +
    ggtitle("Absolute risk") +
    theme_bw() +
    theme(
      legend.position = "bottom"
    ),
  
  
  interim_2 %>%
    mutate(
      obs = exp(obs),
      lci = exp(lci),
      uci = exp(uci)
    ) %>%
    ggplot() + 
    geom_vline(aes(xintercept = 1), colour = "grey") +
    geom_point(aes(y = group, x = obs)) + 
    geom_errorbar(aes(y = group, xmin = lci, xmax = uci, x = obs)) + 
    scale_y_discrete(
      limits= rev(interim_1$group),
      labels = rev(interim_1$label %>% unlist())
    ) +
    xlab("Odds Ratio") +
    scale_x_log10(breaks = c(0.25, 0.5, 1, 2, 4)) +
    coord_cartesian(xlim = c(0.2, 4.1)) +
    ggtitle(paste0("Relative benefit (below 1 = benefit on ", drug1, ")")) +
    theme_bw()
  
) + patchwork::plot_layout(axes = "collect") +
  patchwork::plot_annotation(
    title = paste0(drug1, " vs ", drug2)
  ) &
  theme(
    axis.title.y = element_blank(),
    axis.title.x = element_text(size = 13),
    axis.text = element_text(size = 13)
  )



##: Per drug combination: SGLT2vsSU
drug1 = "SGLT2"; drug2 = "SU"

title <- ATE.overall.no_adj[1,]
title[,] <- NA

interim_1 <- title %>% as.data.frame() %>%mutate(group = paste0("Predicted benefit on ", drug1, ":"), label = paste0("Predicted benefit on ", drug1, ":")) %>%
  rbind(
    ATE.overall.no_adj %>%
      filter(drug_1 == drug1 & drug_2 == drug2 & diff.pred < 0),
    title %>% as.data.frame() %>% mutate(group = paste0("Predicted benefit on ", drug2, ":"), label = paste0("Predicted benefit on ", drug2, ":")),
    ATE.overall.no_adj %>%
      filter(drug_1 == drug1 & drug_2 == drug2 & diff.pred > 0)
  )

title <- odds_ratios.overall.no_adj[1,]
title[,] <- NA


interim_2 <- title %>% as.data.frame() %>%mutate(group = paste0("Predicted benefit on ", drug1, ":"), label = paste0("Predicted benefit on ", drug1, ":")) %>%
  rbind(
    odds_ratios.overall.no_adj %>%
      filter(drug_1 == drug1 & drug_2 == drug2 & obs > 0),
    title %>% as.data.frame() %>% mutate(group = paste0("Predicted benefit on ", drug2, ":"), label = paste0("Predicted benefit on ", drug2, ":")),
    odds_ratios.overall.no_adj %>%
      filter(drug_1 == drug1 & drug_2 == drug2 & obs < 0)
  )

title <- absolute.overall.no_adj[1,]
title[,] <- NA

interim_3 <- title %>% as.data.frame() %>% mutate(group = paste0("Predicted benefit on ", drug1, ":"), label = paste0("Predicted benefit on ", drug1, ":")) %>%
  mutate(drug_type = drug1) %>%
  rbind(
    absolute.overall.no_adj %>%
      filter(drug_1 == drug1 & drug_2 == drug2 & diff.pred < 0),
    title %>% as.data.frame() %>% mutate(group = paste0("Predicted benefit on ", drug2, ":"), label = paste0("Predicted benefit on ", drug2, ":")) %>%
      mutate(drug_type = drug2),
    absolute.overall.no_adj %>%
      filter(drug_1 == drug1 & drug_2 == drug2 & diff.pred > 0)
  )

patchwork::wrap_plots(
  
  interim_1 %>%
    ggplot(aes(y = group, x = obs, xmin = lci, xmax = uci)) +
    geom_vline(aes(xintercept = 0), colour = "grey") +
    geom_point() +
    geom_errorbar() + 
    scale_y_discrete(
      limits= rev(interim_1$group),
      labels = rev(interim_1$label %>% unlist())
    ) +
    scale_x_continuous("Predicted Treatment Discontinuation Difference", labels = scales::percent, breaks = seq(-1, 1, 0.02)) +
    coord_cartesian(xlim = c(-0.15, 0.1)) +
    ggtitle(paste0("Absolute benefit (below 0 = benefit on ", drug1,")")) +
    theme_bw(),
  
  interim_3 %>%
    ggplot(aes(y = group, x = obs, xmin = lci, xmax = uci, colour = drug_type)) +
    geom_point(position = position_dodge(width = 0.5)) +
    geom_errorbar(position = position_dodge(width = 0.5)) +
    scale_y_discrete(
      limits= rev(interim_1$group),
      labels = rev(interim_1$label %>% unlist())
    ) +
    scale_x_continuous("Predicted Discontinuation", labels = scales::percent) +
    scale_colour_manual("Therapy", values = c("GLP1" = "#56B4E9", "SGLT2" = "#E69F00", "SU" = "#CC79A7", "DPP4" = "#0072B2", "TZD" = "#D55E00")) +
    coord_cartesian(xlim = c(0.05, 0.3)) +
    ggtitle("Absolute risk") +
    theme_bw() +
    theme(
      legend.position = "bottom"
    ),
  
  
  interim_2 %>%
    mutate(
      obs = exp(obs),
      lci = exp(lci),
      uci = exp(uci)
    ) %>%
    ggplot() + 
    geom_vline(aes(xintercept = 1), colour = "grey") +
    geom_point(aes(y = group, x = obs)) + 
    geom_errorbar(aes(y = group, xmin = lci, xmax = uci, x = obs)) + 
    scale_y_discrete(
      limits= rev(interim_1$group),
      labels = rev(interim_1$label %>% unlist())
    ) +
    xlab("Odds Ratio") +
    scale_x_log10(breaks = c(0.25, 0.5, 1, 2, 4)) +
    coord_cartesian(xlim = c(0.2, 4.1)) +
    ggtitle(paste0("Relative benefit (below 1 = benefit on ", drug1, ")")) +
    theme_bw()
  
) + patchwork::plot_layout(axes = "collect") +
  patchwork::plot_annotation(
    title = paste0(drug1, " vs ", drug2)
  ) &
  theme(
    axis.title.y = element_blank(),
    axis.title.x = element_text(size = 13),
    axis.text = element_text(size = 13)
  )

##: Per drug combination: SGLT2vsTZD
drug1 = "SGLT2"; drug2 = "TZD"

title <- ATE.overall.no_adj[1,]
title[,] <- NA

interim_1 <- title %>% as.data.frame() %>%mutate(group = paste0("Predicted benefit on ", drug1, ":"), label = paste0("Predicted benefit on ", drug1, ":")) %>%
  rbind(
    ATE.overall.no_adj %>%
      filter(drug_1 == drug1 & drug_2 == drug2 & diff.pred < 0),
    title %>% as.data.frame() %>% mutate(group = paste0("Predicted benefit on ", drug2, ":"), label = paste0("Predicted benefit on ", drug2, ":")),
    ATE.overall.no_adj %>%
      filter(drug_1 == drug1 & drug_2 == drug2 & diff.pred > 0)
  )

title <- odds_ratios.overall.no_adj[1,]
title[,] <- NA


interim_2 <- title %>% as.data.frame() %>%mutate(group = paste0("Predicted benefit on ", drug1, ":"), label = paste0("Predicted benefit on ", drug1, ":")) %>%
  rbind(
    odds_ratios.overall.no_adj %>%
      filter(drug_1 == drug1 & drug_2 == drug2 & obs > 0),
    title %>% as.data.frame() %>% mutate(group = paste0("Predicted benefit on ", drug2, ":"), label = paste0("Predicted benefit on ", drug2, ":")),
    odds_ratios.overall.no_adj %>%
      filter(drug_1 == drug1 & drug_2 == drug2 & obs < 0)
  )

title <- absolute.overall.no_adj[1,]
title[,] <- NA

interim_3 <- title %>% as.data.frame() %>% mutate(group = paste0("Predicted benefit on ", drug1, ":"), label = paste0("Predicted benefit on ", drug1, ":")) %>%
  mutate(drug_type = drug1) %>%
  rbind(
    absolute.overall.no_adj %>%
      filter(drug_1 == drug1 & drug_2 == drug2 & diff.pred < 0),
    title %>% as.data.frame() %>% mutate(group = paste0("Predicted benefit on ", drug2, ":"), label = paste0("Predicted benefit on ", drug2, ":")) %>%
      mutate(drug_type = drug2),
    absolute.overall.no_adj %>%
      filter(drug_1 == drug1 & drug_2 == drug2 & diff.pred > 0)
  )


patchwork::wrap_plots(
  
  interim_1 %>%
    ggplot(aes(y = group, x = obs, xmin = lci, xmax = uci)) +
    geom_vline(aes(xintercept = 0), colour = "grey") +
    geom_point() +
    geom_errorbar() + 
    scale_y_discrete(
      limits= rev(interim_1$group),
      labels = rev(interim_1$label %>% unlist())
    ) +
    scale_x_continuous("Predicted Treatment Discontinuation Difference", labels = scales::percent, breaks = seq(-1, 1, 0.02)) +
    coord_cartesian(xlim = c(-0.15, 0.1)) +
    ggtitle(paste0("Absolute benefit (below 0 = benefit on ", drug1,")")) +
    theme_bw(),
  
  interim_3 %>%
    ggplot(aes(y = group, x = obs, xmin = lci, xmax = uci, colour = drug_type)) +
    geom_point(position = position_dodge(width = 0.5)) +
    geom_errorbar(position = position_dodge(width = 0.5)) +
    scale_y_discrete(
      limits= rev(interim_1$group),
      labels = rev(interim_1$label %>% unlist())
    ) +
    scale_x_continuous("Predicted Discontinuation", labels = scales::percent) +
    scale_colour_manual("Therapy", values = c("GLP1" = "#56B4E9", "SGLT2" = "#E69F00", "SU" = "#CC79A7", "DPP4" = "#0072B2", "TZD" = "#D55E00")) +
    coord_cartesian(xlim = c(0.05, 0.3)) +
    ggtitle("Absolute risk") +
    theme_bw() +
    theme(
      legend.position = "bottom"
    ),
  
  
  interim_2 %>%
    mutate(
      obs = exp(obs),
      lci = exp(lci),
      uci = exp(uci)
    ) %>%
    ggplot() + 
    geom_vline(aes(xintercept = 1), colour = "grey") +
    geom_point(aes(y = group, x = obs)) + 
    geom_errorbar(aes(y = group, xmin = lci, xmax = uci, x = obs)) + 
    scale_y_discrete(
      limits= rev(interim_1$group),
      labels = rev(interim_1$label %>% unlist())
    ) +
    xlab("Odds Ratio") +
    scale_x_log10(breaks = c(0.25, 0.5, 1, 2, 4)) +
    coord_cartesian(xlim = c(0.2, 4.1)) +
    ggtitle(paste0("Relative benefit (below 1 = benefit on ", drug1, ")")) +
    theme_bw()
  
) + patchwork::plot_layout(axes = "collect") +
  patchwork::plot_annotation(
    title = paste0(drug1, " vs ", drug2)
  ) &
  theme(
    axis.title.y = element_blank(),
    axis.title.x = element_text(size = 13),
    axis.text = element_text(size = 13)
  )


##: Per drug combination: GLP1vsDPP4
drug1 = "GLP1"; drug2 = "DPP4"

title <- ATE.overall.no_adj[1,]
title[,] <- NA

interim_1 <- title %>% as.data.frame() %>%mutate(group = paste0("Predicted benefit on ", drug1, ":"), label = paste0("Predicted benefit on ", drug1, ":")) %>%
  rbind(
    ATE.overall.no_adj %>%
      filter(drug_1 == drug1 & drug_2 == drug2 & diff.pred < 0),
    title %>% as.data.frame() %>% mutate(group = paste0("Predicted benefit on ", drug2, ":"), label = paste0("Predicted benefit on ", drug2, ":")),
    ATE.overall.no_adj %>%
      filter(drug_1 == drug1 & drug_2 == drug2 & diff.pred > 0)
  )

title <- odds_ratios.overall.no_adj[1,]
title[,] <- NA


interim_2 <- title %>% as.data.frame() %>%mutate(group = paste0("Predicted benefit on ", drug1, ":"), label = paste0("Predicted benefit on ", drug1, ":")) %>%
  rbind(
    odds_ratios.overall.no_adj %>%
      filter(drug_1 == drug1 & drug_2 == drug2 & obs > 0),
    title %>% as.data.frame() %>% mutate(group = paste0("Predicted benefit on ", drug2, ":"), label = paste0("Predicted benefit on ", drug2, ":")),
    odds_ratios.overall.no_adj %>%
      filter(drug_1 == drug1 & drug_2 == drug2 & obs < 0)
  )

title <- absolute.overall.no_adj[1,]
title[,] <- NA

interim_3 <- title %>% as.data.frame() %>% mutate(group = paste0("Predicted benefit on ", drug1, ":"), label = paste0("Predicted benefit on ", drug1, ":")) %>%
  mutate(drug_type = drug1) %>%
  rbind(
    absolute.overall.no_adj %>%
      filter(drug_1 == drug1 & drug_2 == drug2 & diff.pred < 0),
    title %>% as.data.frame() %>% mutate(group = paste0("Predicted benefit on ", drug2, ":"), label = paste0("Predicted benefit on ", drug2, ":")) %>%
      mutate(drug_type = drug2),
    absolute.overall.no_adj %>%
      filter(drug_1 == drug1 & drug_2 == drug2 & diff.pred > 0)
  )


patchwork::wrap_plots(
  
  interim_1 %>%
    ggplot(aes(y = group, x = obs, xmin = lci, xmax = uci)) +
    geom_vline(aes(xintercept = 0), colour = "grey") +
    geom_point() +
    geom_errorbar() + 
    scale_y_discrete(
      limits= rev(interim_1$group),
      labels = rev(interim_1$label %>% unlist())
    ) +
    scale_x_continuous("Predicted Treatment Discontinuation Difference", labels = scales::percent, breaks = seq(-1, 1, 0.02)) +
    coord_cartesian(xlim = c(-0.15, 0.1)) +
    ggtitle(paste0("Absolute benefit (below 0 = benefit on ", drug1,")")) +
    theme_bw(),
  
  interim_3 %>%
    ggplot(aes(y = group, x = obs, xmin = lci, xmax = uci, colour = drug_type)) +
    geom_point(position = position_dodge(width = 0.5)) +
    geom_errorbar(position = position_dodge(width = 0.5)) +
    scale_y_discrete(
      limits= rev(interim_1$group),
      labels = rev(interim_1$label %>% unlist())
    ) +
    scale_x_continuous("Predicted Discontinuation", labels = scales::percent) +
    scale_colour_manual("Therapy", values = c("GLP1" = "#56B4E9", "SGLT2" = "#E69F00", "SU" = "#CC79A7", "DPP4" = "#0072B2", "TZD" = "#D55E00")) +
    coord_cartesian(xlim = c(0.05, 0.3)) +
    ggtitle("Absolute risk") +
    theme_bw() +
    theme(
      legend.position = "bottom"
    ),
  
  
  interim_2 %>%
    mutate(
      obs = exp(obs),
      lci = exp(lci),
      uci = exp(uci)
    ) %>%
    ggplot() + 
    geom_vline(aes(xintercept = 1), colour = "grey") +
    geom_point(aes(y = group, x = obs)) + 
    geom_errorbar(aes(y = group, xmin = lci, xmax = uci, x = obs)) + 
    scale_y_discrete(
      limits= rev(interim_1$group),
      labels = rev(interim_1$label %>% unlist())
    ) +
    xlab("Odds Ratio") +
    scale_x_log10(breaks = c(0.25, 0.5, 1, 2, 4)) +
    coord_cartesian(xlim = c(0.2, 4.1)) +
    ggtitle(paste0("Relative benefit (below 1 = benefit on ", drug1, ")")) +
    theme_bw()
  
) + patchwork::plot_layout(axes = "collect") +
  patchwork::plot_annotation(
    title = paste0(drug1, " vs ", drug2)
  ) &
  theme(
    axis.title.y = element_blank(),
    axis.title.x = element_text(size = 13),
    axis.text = element_text(size = 13)
  )


##: Per drug combination: GLP1vsSU
drug1 = "GLP1"; drug2 = "SU"

title <- ATE.overall.no_adj[1,]
title[,] <- NA

interim_1 <- title %>% as.data.frame() %>%mutate(group = paste0("Predicted benefit on ", drug1, ":"), label = paste0("Predicted benefit on ", drug1, ":")) %>%
  rbind(
    ATE.overall.no_adj %>%
      filter(drug_1 == drug1 & drug_2 == drug2 & diff.pred < 0),
    title %>% as.data.frame() %>% mutate(group = paste0("Predicted benefit on ", drug2, ":"), label = paste0("Predicted benefit on ", drug2, ":")),
    ATE.overall.no_adj %>%
      filter(drug_1 == drug1 & drug_2 == drug2 & diff.pred > 0)
  )

title <- odds_ratios.overall.no_adj[1,]
title[,] <- NA


interim_2 <- title %>% as.data.frame() %>%mutate(group = paste0("Predicted benefit on ", drug1, ":"), label = paste0("Predicted benefit on ", drug1, ":")) %>%
  rbind(
    odds_ratios.overall.no_adj %>%
      filter(drug_1 == drug1 & drug_2 == drug2 & obs > 0),
    title %>% as.data.frame() %>% mutate(group = paste0("Predicted benefit on ", drug2, ":"), label = paste0("Predicted benefit on ", drug2, ":")),
    odds_ratios.overall.no_adj %>%
      filter(drug_1 == drug1 & drug_2 == drug2 & obs < 0)
  )

title <- absolute.overall.no_adj[1,]
title[,] <- NA

interim_3 <- title %>% as.data.frame() %>% mutate(group = paste0("Predicted benefit on ", drug1, ":"), label = paste0("Predicted benefit on ", drug1, ":")) %>%
  mutate(drug_type = drug1) %>%
  rbind(
    absolute.overall.no_adj %>%
      filter(drug_1 == drug1 & drug_2 == drug2 & diff.pred < 0),
    title %>% as.data.frame() %>% mutate(group = paste0("Predicted benefit on ", drug2, ":"), label = paste0("Predicted benefit on ", drug2, ":")) %>%
      mutate(drug_type = drug2),
    absolute.overall.no_adj %>%
      filter(drug_1 == drug1 & drug_2 == drug2 & diff.pred > 0)
  )


patchwork::wrap_plots(
  
  interim_1 %>%
    ggplot(aes(y = group, x = obs, xmin = lci, xmax = uci)) +
    geom_vline(aes(xintercept = 0), colour = "grey") +
    geom_point() +
    geom_errorbar() + 
    scale_y_discrete(
      limits= rev(interim_1$group),
      labels = rev(interim_1$label %>% unlist())
    ) +
    scale_x_continuous("Predicted Treatment Discontinuation Difference", labels = scales::percent, breaks = seq(-1, 1, 0.02)) +
    coord_cartesian(xlim = c(-0.15, 0.1)) +
    ggtitle(paste0("Absolute benefit (below 0 = benefit on ", drug1,")")) +
    theme_bw(),
  
  interim_3 %>%
    ggplot(aes(y = group, x = obs, xmin = lci, xmax = uci, colour = drug_type)) +
    geom_point(position = position_dodge(width = 0.5)) +
    geom_errorbar(position = position_dodge(width = 0.5)) +
    scale_y_discrete(
      limits= rev(interim_1$group),
      labels = rev(interim_1$label %>% unlist())
    ) +
    scale_x_continuous("Predicted Discontinuation", labels = scales::percent) +
    scale_colour_manual("Therapy", values = c("GLP1" = "#56B4E9", "SGLT2" = "#E69F00", "SU" = "#CC79A7", "DPP4" = "#0072B2", "TZD" = "#D55E00")) +
    coord_cartesian(xlim = c(0.05, 0.3)) +
    ggtitle("Absolute risk") +
    theme_bw() +
    theme(
      legend.position = "bottom"
    ),
  
  
  interim_2 %>%
    mutate(
      obs = exp(obs),
      lci = exp(lci),
      uci = exp(uci)
    ) %>%
    ggplot() + 
    geom_vline(aes(xintercept = 1), colour = "grey") +
    geom_point(aes(y = group, x = obs)) + 
    geom_errorbar(aes(y = group, xmin = lci, xmax = uci, x = obs)) + 
    scale_y_discrete(
      limits= rev(interim_1$group),
      labels = rev(interim_1$label %>% unlist())
    ) +
    xlab("Odds Ratio") +
    scale_x_log10(breaks = c(0.25, 0.5, 1, 2, 4)) +
    coord_cartesian(xlim = c(0.2, 4.1)) +
    ggtitle(paste0("Relative benefit (below 1 = benefit on ", drug1, ")")) +
    theme_bw()
  
) + patchwork::plot_layout(axes = "collect") +
  patchwork::plot_annotation(
    title = paste0(drug1, " vs ", drug2)
  ) &
  theme(
    axis.title.y = element_blank(),
    axis.title.x = element_text(size = 13),
    axis.text = element_text(size = 13)
  )


##: Per drug combination: GLP1vsTZD
drug1 = "GLP1"; drug2 = "TZD"

title <- ATE.overall.no_adj[1,]
title[,] <- NA

interim_1 <- title %>% as.data.frame() %>%mutate(group = paste0("Predicted benefit on ", drug1, ":"), label = paste0("Predicted benefit on ", drug1, ":")) %>%
  rbind(
    ATE.overall.no_adj %>%
      filter(drug_1 == drug1 & drug_2 == drug2 & diff.pred < 0),
    title %>% as.data.frame() %>% mutate(group = paste0("Predicted benefit on ", drug2, ":"), label = paste0("Predicted benefit on ", drug2, ":")),
    ATE.overall.no_adj %>%
      filter(drug_1 == drug1 & drug_2 == drug2 & diff.pred > 0)
  )

title <- odds_ratios.overall.no_adj[1,]
title[,] <- NA


interim_2 <- title %>% as.data.frame() %>%mutate(group = paste0("Predicted benefit on ", drug1, ":"), label = paste0("Predicted benefit on ", drug1, ":")) %>%
  rbind(
    odds_ratios.overall.no_adj %>%
      filter(drug_1 == drug1 & drug_2 == drug2 & obs > 0),
    title %>% as.data.frame() %>% mutate(group = paste0("Predicted benefit on ", drug2, ":"), label = paste0("Predicted benefit on ", drug2, ":")),
    odds_ratios.overall.no_adj %>%
      filter(drug_1 == drug1 & drug_2 == drug2 & obs < 0)
  )

title <- absolute.overall.no_adj[1,]
title[,] <- NA

interim_3 <- title %>% as.data.frame() %>% mutate(group = paste0("Predicted benefit on ", drug1, ":"), label = paste0("Predicted benefit on ", drug1, ":")) %>%
  mutate(drug_type = drug1) %>%
  rbind(
    absolute.overall.no_adj %>%
      filter(drug_1 == drug1 & drug_2 == drug2 & diff.pred < 0),
    title %>% as.data.frame() %>% mutate(group = paste0("Predicted benefit on ", drug2, ":"), label = paste0("Predicted benefit on ", drug2, ":")) %>%
      mutate(drug_type = drug2),
    absolute.overall.no_adj %>%
      filter(drug_1 == drug1 & drug_2 == drug2 & diff.pred > 0)
  )


patchwork::wrap_plots(
  
  interim_1 %>%
    ggplot(aes(y = group, x = obs, xmin = lci, xmax = uci)) +
    geom_vline(aes(xintercept = 0), colour = "grey") +
    geom_point() +
    geom_errorbar() + 
    scale_y_discrete(
      limits= rev(interim_1$group),
      labels = rev(interim_1$label %>% unlist())
    ) +
    scale_x_continuous("Predicted Treatment Discontinuation Difference", labels = scales::percent, breaks = seq(-1, 1, 0.02)) +
    coord_cartesian(xlim = c(-0.15, 0.1)) +
    ggtitle(paste0("Absolute benefit (below 0 = benefit on ", drug1,")")) +
    theme_bw(),
  
  interim_3 %>%
    ggplot(aes(y = group, x = obs, xmin = lci, xmax = uci, colour = drug_type)) +
    geom_point(position = position_dodge(width = 0.5)) +
    geom_errorbar(position = position_dodge(width = 0.5)) +
    scale_y_discrete(
      limits= rev(interim_1$group),
      labels = rev(interim_1$label %>% unlist())
    ) +
    scale_x_continuous("Predicted Discontinuation", labels = scales::percent) +
    scale_colour_manual("Therapy", values = c("GLP1" = "#56B4E9", "SGLT2" = "#E69F00", "SU" = "#CC79A7", "DPP4" = "#0072B2", "TZD" = "#D55E00")) +
    coord_cartesian(xlim = c(0.05, 0.3)) +
    ggtitle("Absolute risk") +
    theme_bw() +
    theme(
      legend.position = "bottom"
    ),
  
  
  interim_2 %>%
    mutate(
      obs = exp(obs),
      lci = exp(lci),
      uci = exp(uci)
    ) %>%
    ggplot() + 
    geom_vline(aes(xintercept = 1), colour = "grey") +
    geom_point(aes(y = group, x = obs)) + 
    geom_errorbar(aes(y = group, xmin = lci, xmax = uci, x = obs)) + 
    scale_y_discrete(
      limits= rev(interim_1$group),
      labels = rev(interim_1$label %>% unlist())
    ) +
    xlab("Odds Ratio") +
    scale_x_log10(breaks = c(0.25, 0.5, 1, 2, 4)) +
    coord_cartesian(xlim = c(0.2, 4.1)) +
    ggtitle(paste0("Relative benefit (below 1 = benefit on ", drug1, ")")) +
    theme_bw()
  
) + patchwork::plot_layout(axes = "collect") +
  patchwork::plot_annotation(
    title = paste0(drug1, " vs ", drug2)
  ) &
  theme(
    axis.title.y = element_blank(),
    axis.title.x = element_text(size = 13),
    axis.text = element_text(size = 13)
  )


##: Per drug combination: TZDvsSU
drug1 = "TZD"; drug2 = "SU"

title <- ATE.overall.no_adj[1,]
title[,] <- NA

interim_1 <- title %>% as.data.frame() %>%mutate(group = paste0("Predicted benefit on ", drug1, ":"), label = paste0("Predicted benefit on ", drug1, ":")) %>%
  rbind(
    ATE.overall.no_adj %>%
      filter(drug_1 == drug1 & drug_2 == drug2 & diff.pred < 0),
    title %>% as.data.frame() %>% mutate(group = paste0("Predicted benefit on ", drug2, ":"), label = paste0("Predicted benefit on ", drug2, ":")),
    ATE.overall.no_adj %>%
      filter(drug_1 == drug1 & drug_2 == drug2 & diff.pred > 0)
  )

title <- odds_ratios.overall.no_adj[1,]
title[,] <- NA


interim_2 <- title %>% as.data.frame() %>%mutate(group = paste0("Predicted benefit on ", drug1, ":"), label = paste0("Predicted benefit on ", drug1, ":")) %>%
  rbind(
    odds_ratios.overall.no_adj %>%
      filter(drug_1 == drug1 & drug_2 == drug2 & obs > 0),
    title %>% as.data.frame() %>% mutate(group = paste0("Predicted benefit on ", drug2, ":"), label = paste0("Predicted benefit on ", drug2, ":")),
    odds_ratios.overall.no_adj %>%
      filter(drug_1 == drug1 & drug_2 == drug2 & obs < 0)
  )

title <- absolute.overall.no_adj[1,]
title[,] <- NA

interim_3 <- title %>% as.data.frame() %>% mutate(group = paste0("Predicted benefit on ", drug1, ":"), label = paste0("Predicted benefit on ", drug1, ":")) %>%
  mutate(drug_type = drug1) %>%
  rbind(
    absolute.overall.no_adj %>%
      filter(drug_1 == drug1 & drug_2 == drug2 & diff.pred < 0),
    title %>% as.data.frame() %>% mutate(group = paste0("Predicted benefit on ", drug2, ":"), label = paste0("Predicted benefit on ", drug2, ":")) %>%
      mutate(drug_type = drug2),
    absolute.overall.no_adj %>%
      filter(drug_1 == drug1 & drug_2 == drug2 & diff.pred > 0)
  )


patchwork::wrap_plots(
  
  interim_1 %>%
    ggplot(aes(y = group, x = obs, xmin = lci, xmax = uci)) +
    geom_vline(aes(xintercept = 0), colour = "grey") +
    geom_point() +
    geom_errorbar() + 
    scale_y_discrete(
      limits= rev(interim_1$group),
      labels = rev(interim_1$label %>% unlist())
    ) +
    scale_x_continuous("Predicted Treatment Discontinuation Difference", labels = scales::percent, breaks = seq(-1, 1, 0.02)) +
    coord_cartesian(xlim = c(-0.15, 0.1)) +
    ggtitle(paste0("Absolute benefit (below 0 = benefit on ", drug1,")")) +
    theme_bw(),
  
  interim_3 %>%
    ggplot(aes(y = group, x = obs, xmin = lci, xmax = uci, colour = drug_type)) +
    geom_point(position = position_dodge(width = 0.5)) +
    geom_errorbar(position = position_dodge(width = 0.5)) +
    scale_y_discrete(
      limits= rev(interim_1$group),
      labels = rev(interim_1$label %>% unlist())
    ) +
    scale_x_continuous("Predicted Discontinuation", labels = scales::percent) +
    scale_colour_manual("Therapy", values = c("GLP1" = "#56B4E9", "SGLT2" = "#E69F00", "SU" = "#CC79A7", "DPP4" = "#0072B2", "TZD" = "#D55E00")) +
    coord_cartesian(xlim = c(0.05, 0.3)) +
    ggtitle("Absolute risk") +
    theme_bw() +
    theme(
      legend.position = "bottom"
    ),
  
  
  interim_2 %>%
    mutate(
      obs = exp(obs),
      lci = exp(lci),
      uci = exp(uci)
    ) %>%
    ggplot() + 
    geom_vline(aes(xintercept = 1), colour = "grey") +
    geom_point(aes(y = group, x = obs)) + 
    geom_errorbar(aes(y = group, xmin = lci, xmax = uci, x = obs)) + 
    scale_y_discrete(
      limits= rev(interim_1$group),
      labels = rev(interim_1$label %>% unlist())
    ) +
    xlab("Odds Ratio") +
    scale_x_log10(breaks = c(0.25, 0.5, 1, 2, 4)) +
    coord_cartesian(xlim = c(0.2, 4.1)) +
    ggtitle(paste0("Relative benefit (below 1 = benefit on ", drug1, ")")) +
    theme_bw()
  
) + patchwork::plot_layout(axes = "collect") +
  patchwork::plot_annotation(
    title = paste0(drug1, " vs ", drug2)
  ) &
  theme(
    axis.title.y = element_blank(),
    axis.title.x = element_text(size = 13),
    axis.text = element_text(size = 13)
  )


##: Per drug combination: TZDvsDPP4
drug1 = "TZD"; drug2 = "DPP4"

title <- ATE.overall.no_adj[1,]
title[,] <- NA

interim_1 <- title %>% as.data.frame() %>%mutate(group = paste0("Predicted benefit on ", drug1, ":"), label = paste0("Predicted benefit on ", drug1, ":")) %>%
  rbind(
    ATE.overall.no_adj %>%
      filter(drug_1 == drug1 & drug_2 == drug2 & diff.pred < 0),
    title %>% as.data.frame() %>% mutate(group = paste0("Predicted benefit on ", drug2, ":"), label = paste0("Predicted benefit on ", drug2, ":")),
    ATE.overall.no_adj %>%
      filter(drug_1 == drug1 & drug_2 == drug2 & diff.pred > 0)
  )

title <- odds_ratios.overall.no_adj[1,]
title[,] <- NA


interim_2 <- title %>% as.data.frame() %>%mutate(group = paste0("Predicted benefit on ", drug1, ":"), label = paste0("Predicted benefit on ", drug1, ":")) %>%
  rbind(
    odds_ratios.overall.no_adj %>%
      filter(drug_1 == drug1 & drug_2 == drug2 & obs > 0),
    title %>% as.data.frame() %>% mutate(group = paste0("Predicted benefit on ", drug2, ":"), label = paste0("Predicted benefit on ", drug2, ":")),
    odds_ratios.overall.no_adj %>%
      filter(drug_1 == drug1 & drug_2 == drug2 & obs < 0)
  )

title <- absolute.overall.no_adj[1,]
title[,] <- NA

interim_3 <- title %>% as.data.frame() %>% mutate(group = paste0("Predicted benefit on ", drug1, ":"), label = paste0("Predicted benefit on ", drug1, ":")) %>%
  mutate(drug_type = drug1) %>%
  rbind(
    absolute.overall.no_adj %>%
      filter(drug_1 == drug1 & drug_2 == drug2 & diff.pred < 0),
    title %>% as.data.frame() %>% mutate(group = paste0("Predicted benefit on ", drug2, ":"), label = paste0("Predicted benefit on ", drug2, ":")) %>%
      mutate(drug_type = drug2),
    absolute.overall.no_adj %>%
      filter(drug_1 == drug1 & drug_2 == drug2 & diff.pred > 0)
  )


patchwork::wrap_plots(
  
  interim_1 %>%
    ggplot(aes(y = group, x = obs, xmin = lci, xmax = uci)) +
    geom_vline(aes(xintercept = 0), colour = "grey") +
    geom_point() +
    geom_errorbar() + 
    scale_y_discrete(
      limits= rev(interim_1$group),
      labels = rev(interim_1$label %>% unlist())
    ) +
    scale_x_continuous("Predicted Treatment Discontinuation Difference", labels = scales::percent, breaks = seq(-1, 1, 0.02)) +
    coord_cartesian(xlim = c(-0.15, 0.1)) +
    ggtitle(paste0("Absolute benefit (below 0 = benefit on ", drug1,")")) +
    theme_bw(),
  
  interim_3 %>%
    ggplot(aes(y = group, x = obs, xmin = lci, xmax = uci, colour = drug_type)) +
    geom_point(position = position_dodge(width = 0.5)) +
    geom_errorbar(position = position_dodge(width = 0.5)) +
    scale_y_discrete(
      limits= rev(interim_1$group),
      labels = rev(interim_1$label %>% unlist())
    ) +
    scale_x_continuous("Predicted Discontinuation", labels = scales::percent) +
    scale_colour_manual("Therapy", values = c("GLP1" = "#56B4E9", "SGLT2" = "#E69F00", "SU" = "#CC79A7", "DPP4" = "#0072B2", "TZD" = "#D55E00")) +
    coord_cartesian(xlim = c(0.05, 0.3)) +
    ggtitle("Absolute risk") +
    theme_bw() +
    theme(
      legend.position = "bottom"
    ),
  
  
  interim_2 %>%
    mutate(
      obs = exp(obs),
      lci = exp(lci),
      uci = exp(uci)
    ) %>%
    ggplot() + 
    geom_vline(aes(xintercept = 1), colour = "grey") +
    geom_point(aes(y = group, x = obs)) + 
    geom_errorbar(aes(y = group, xmin = lci, xmax = uci, x = obs)) + 
    scale_y_discrete(
      limits= rev(interim_1$group),
      labels = rev(interim_1$label %>% unlist())
    ) +
    xlab("Odds Ratio") +
    scale_x_log10(breaks = c(0.25, 0.5, 1, 2, 4)) +
    coord_cartesian(xlim = c(0.2, 4.1)) +
    ggtitle(paste0("Relative benefit (below 1 = benefit on ", drug1, ")")) +
    theme_bw()
  
) + patchwork::plot_layout(axes = "collect") +
  patchwork::plot_annotation(
    title = paste0(drug1, " vs ", drug2)
  ) &
  theme(
    axis.title.y = element_blank(),
    axis.title.x = element_text(size = 13),
    axis.text = element_text(size = 13)
  )


##: Per drug combination: SUvsDPP4
drug1 = "SU"; drug2 = "DPP4"

title <- ATE.overall.no_adj[1,]
title[,] <- NA

interim_1 <- title %>% as.data.frame() %>%mutate(group = paste0("Predicted benefit on ", drug1, ":"), label = paste0("Predicted benefit on ", drug1, ":")) %>%
  rbind(
    ATE.overall.no_adj %>%
      filter(drug_1 == drug1 & drug_2 == drug2 & diff.pred < 0),
    title %>% as.data.frame() %>% mutate(group = paste0("Predicted benefit on ", drug2, ":"), label = paste0("Predicted benefit on ", drug2, ":")),
    ATE.overall.no_adj %>%
      filter(drug_1 == drug1 & drug_2 == drug2 & diff.pred > 0)
  )

title <- odds_ratios.overall.no_adj[1,]
title[,] <- NA


interim_2 <- title %>% as.data.frame() %>%mutate(group = paste0("Predicted benefit on ", drug1, ":"), label = paste0("Predicted benefit on ", drug1, ":")) %>%
  rbind(
    odds_ratios.overall.no_adj %>%
      filter(drug_1 == drug1 & drug_2 == drug2 & obs > 0),
    title %>% as.data.frame() %>% mutate(group = paste0("Predicted benefit on ", drug2, ":"), label = paste0("Predicted benefit on ", drug2, ":")),
    odds_ratios.overall.no_adj %>%
      filter(drug_1 == drug1 & drug_2 == drug2 & obs < 0)
  )

title <- absolute.overall.no_adj[1,]
title[,] <- NA

interim_3 <- title %>% as.data.frame() %>% mutate(group = paste0("Predicted benefit on ", drug1, ":"), label = paste0("Predicted benefit on ", drug1, ":")) %>%
  mutate(drug_type = drug1) %>%
  rbind(
    absolute.overall.no_adj %>%
      filter(drug_1 == drug1 & drug_2 == drug2 & diff.pred < 0),
    title %>% as.data.frame() %>% mutate(group = paste0("Predicted benefit on ", drug2, ":"), label = paste0("Predicted benefit on ", drug2, ":")) %>%
      mutate(drug_type = drug2),
    absolute.overall.no_adj %>%
      filter(drug_1 == drug1 & drug_2 == drug2 & diff.pred > 0)
  )



patchwork::wrap_plots(
  
  interim_1 %>%
    ggplot(aes(y = group, x = obs, xmin = lci, xmax = uci)) +
    geom_vline(aes(xintercept = 0), colour = "grey") +
    geom_point() +
    geom_errorbar() + 
    scale_y_discrete(
      limits= rev(interim_1$group),
      labels = rev(interim_1$label %>% unlist())
    ) +
    scale_x_continuous("Predicted Treatment Discontinuation Difference", labels = scales::percent, breaks = seq(-1, 1, 0.02)) +
    coord_cartesian(xlim = c(-0.15, 0.1)) +
    ggtitle(paste0("Absolute benefit (below 0 = benefit on ", drug1,")")) +
    theme_bw(),
  
  interim_3 %>%
    ggplot(aes(y = group, x = obs, xmin = lci, xmax = uci, colour = drug_type)) +
    geom_point(position = position_dodge(width = 0.5)) +
    geom_errorbar(position = position_dodge(width = 0.5)) +
    scale_y_discrete(
      limits= rev(interim_1$group),
      labels = rev(interim_1$label %>% unlist())
    ) +
    scale_x_continuous("Predicted Discontinuation", labels = scales::percent) +
    scale_colour_manual("Therapy", values = c("GLP1" = "#56B4E9", "SGLT2" = "#E69F00", "SU" = "#CC79A7", "DPP4" = "#0072B2", "TZD" = "#D55E00")) +
    coord_cartesian(xlim = c(0.05, 0.3)) +
    ggtitle("Absolute risk") +
    theme_bw() +
    theme(
      legend.position = "bottom"
    ),
  
  
  interim_2 %>%
    mutate(
      obs = exp(obs),
      lci = exp(lci),
      uci = exp(uci)
    ) %>%
    ggplot() + 
    geom_vline(aes(xintercept = 1), colour = "grey") +
    geom_point(aes(y = group, x = obs)) + 
    geom_errorbar(aes(y = group, xmin = lci, xmax = uci, x = obs)) + 
    scale_y_discrete(
      limits= rev(interim_1$group),
      labels = rev(interim_1$label %>% unlist())
    ) +
    xlab("Odds Ratio") +
    scale_x_log10(breaks = c(0.25, 0.5, 1, 2, 4)) +
    coord_cartesian(xlim = c(0.2, 4.1)) +
    ggtitle(paste0("Relative benefit (below 1 = benefit on ", drug1, ")")) +
    theme_bw()
  
) + patchwork::plot_layout(axes = "collect") +
  patchwork::plot_annotation(
    title = paste0(drug1, " vs ", drug2)
  ) &
  theme(
    axis.title.y = element_blank(),
    axis.title.x = element_text(size = 13),
    axis.text = element_text(size = 13)
  )



dev.off()




