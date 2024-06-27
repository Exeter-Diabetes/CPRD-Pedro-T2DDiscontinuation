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



###############################################################################
###############################################################################
############################### Read Data In ##################################
###############################################################################
###############################################################################

# load dataset
cprd_dataset.3m <- set_up_data(
  raw_data = "20240308_t2d_1stinstance",
  diagnosis = FALSE,
  therapies = c("GLP1", "DPP4", "SGLT2", "TZD", "SU"),
  dataset = "3m.disc.dataset.dev",
  follow_up = "6-months",
  full_prescribing_history = TRUE
) %>%
  drop_na(-stopdrug_6m_6mFU, -stopdrug_12m_6mFU)

# load dataset
cprd_dataset.6m <- set_up_data(
  raw_data = "20240308_t2d_1stinstance",
  diagnosis = FALSE,
  therapies = c("GLP1", "DPP4", "SGLT2", "TZD", "SU"),
  dataset = "6m.disc.dataset.dev",
  follow_up = "6-months",
  full_prescribing_history = TRUE
) %>%
  drop_na(-stopdrug_12m_6mFU)

# load dataset
cprd_dataset.12m <- set_up_data(
  raw_data = "20240308_t2d_1stinstance",
  diagnosis = FALSE,
  therapies = c("GLP1", "DPP4", "SGLT2", "TZD", "SU"),
  dataset = "12m.disc.dataset.dev",
  follow_up = "6-months",
  full_prescribing_history = TRUE
) %>%
  drop_na()

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


###############################################################################
###############################################################################
######################## Fit the model: whole cohort ##########################
###############################################################################
###############################################################################


# Run model
if (class(try(
  
  bartmachine_full_model_3m <- readRDS("results/Models/bartmachine/bartmachine_full_model_3m.rds")
  
  , silent = TRUE)) == "try-error") {
  
  bartmachine_full_model_3m <- bartMachine::build_bart_machine(X = cprd_dataset.3m %>%
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
                                                                     prehdl,
                                                                     stopdrug_3m_3mFU_MFN_hist),
                                                            y = cprd_dataset.3m %>% 
                                                              select(stopdrug_3m_6mFU) %>%
                                                              mutate(stopdrug_3m_6mFU = factor(stopdrug_3m_6mFU)) %>%
                                                              unlist(),
                                                            num_trees = 200,
                                                            use_missing_data = FALSE,
                                                            num_burn_in = 15000,
                                                            num_iterations_after_burn_in = 10000,
                                                            serialize = TRUE)
  
  saveRDS(bartmachine_full_model_3m, "results/Models/bartmachine/bartmachine_full_model_3m.rds")
  
}



plot_pooled_3m.overall <- calibration_plot(
  data = cprd_dataset.3m %>% mutate(pred = bartmachine_full_model_3m$p_hat) %>% mutate(stopdrug_3m_6mFU = as.numeric(stopdrug_3m_6mFU)-1) %>% mutate(grouping = "Pooled"), 
  obs = "stopdrug_3m_6mFU", 
  pred = "pred",
  group = "grouping",
  nTiles = 10)$calibration_plot +
  scale_x_continuous("Prediction", limits = c(0, 0.7), breaks = seq(0, 0.7, 0.1)) +
  scale_y_continuous("Observation", limits = c(0, 0.7), breaks = seq(0, 0.7, 0.1)) +
  theme_bw() +
  scale_colour_manual(values = c("Pooled" = "black", "SGLT2" = "#E69F00", "GLP1" = "#56B4E9", "SU" = "#CC79A7", "DPP4" = "#0072B2", "TZD" = "#D55E00", "MFN" = "grey"), breaks = rev(c("Pooled", "GLP1", "DPP4", "SGLT2", "TZD", "SU")), labels = rev(c("Pooled", "GLP-1RA", "DPP4i", "SGLT2i", "TZD", "SU")), name = "Therapy", guide = guide_legend(reverse = TRUE, nrow = 1)) +
  theme(
    legend.position = "bottom"
  )



plot_pooled_3m.drugs <- calibration_plot(
  data = cprd_dataset.3m %>% mutate(pred = bartmachine_full_model_3m$p_hat) %>% mutate(stopdrug_3m_6mFU = as.numeric(stopdrug_3m_6mFU)-1), 
  obs = "stopdrug_3m_6mFU", 
  pred = "pred", 
  group = "drugclass",
  nTiles = 10)$calibration_plot +
  scale_x_continuous("Prediction", limits = c(0, 0.7), breaks = seq(0, 0.7, 0.1)) +
  scale_y_continuous("Observation", limits = c(0, 0.7), breaks = seq(0, 0.7, 0.1)) +
  theme_bw() +
  scale_colour_manual(values = c("Pooled" = "black", "SGLT2" = "#E69F00", "GLP1" = "#56B4E9", "SU" = "#CC79A7", "DPP4" = "#0072B2", "TZD" = "#D55E00", "MFN" = "grey"), breaks = rev(c("Pooled", "GLP1", "DPP4", "SGLT2", "TZD", "SU")), labels = rev(c("Pooled", "GLP-1RA", "DPP4i", "SGLT2i", "TZD", "SU")), name = "Therapy", guide = guide_legend(reverse = TRUE, nrow = 1)) +
  theme(
    legend.position = "bottom"
  )


roc_values_3m_overall <- pROC::roc(response = cprd_dataset.3m %>%
                                     select(stopdrug_3m_6mFU) %>%
                                     mutate(stopdrug_3m_6mFU = factor(stopdrug_3m_6mFU)) %>%
                                     unlist(),
                                   predictor = bartmachine_full_model_3m$p_hat, ci = TRUE)$ci %>%
  as.data.frame() %>%
  unlist()


# roc_coords_3m <- pROC::ci(pROC::roc(
#   cprd_dataset.3m %>%
#     select(stopdrug_3m_6mFU) %>%
#     mutate(stopdrug_3m_6mFU = factor(stopdrug_3m_6mFU)) %>%
#     unlist(),
#   bartmachine_full_model_3m$p_hat
# ), of = "thresholds")
# 
# saveRDS(roc_coords_3m, "results/Models/bartmachine/roc_coords_3m.rds")

roc_coords_3m <- readRDS("results/Models/bartmachine/roc_coords_3m.rds")


plot_roc_3m <- roc_coords_3m %>%
  as.data.frame() %>%
  ggplot() +
  geom_abline(aes(intercept = 0, slope = 1), colour = "black", linetype = "dashed") +
  geom_path(aes(x = 1 - specificity.2.5., y = sensitivity.2.5.), alpha = 0.4) +
  geom_path(aes(x = 1 - specificity.50., y = sensitivity.50.)) +
  geom_path(aes(x = 1 - specificity.97.5., y = sensitivity.97.5.), alpha = 0.4) +
  annotate("label", x = 0.60, y = 0.125, size = 5, label = paste0("AUROC=", signif(roc_values_3m_overall[2], digits = 3), ")")) + 
  scale_x_continuous("1 - Specificity", limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
  scale_y_continuous("Sensitivity", limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
  theme_bw() +
  theme(
    legend.position = "bottom",
    axis.title = element_text(size = 18),
    axis.text = element_text(size = 16)
  )


pdf("results/figures/07.roc_3m_overall.pdf", width = 5, height = 5)

plot_roc_3m

dev.off()



# Run model
if (class(try(
  
  bartmachine_full_model_6m <- readRDS("results/Models/bartmachine/bartmachine_full_model_6m.rds")
  
  , silent = TRUE)) == "try-error") {
  
  bartmachine_full_model_6m <- bartMachine::build_bart_machine(X = cprd_dataset.6m %>%
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
                                                                        prehdl,
                                                                        stopdrug_3m_3mFU_MFN_hist),
                                                               y = cprd_dataset.6m %>% 
                                                                 select(stopdrug_6m_6mFU) %>%
                                                                 mutate(stopdrug_6m_6mFU = factor(stopdrug_6m_6mFU)) %>%
                                                                 unlist(),
                                                               num_trees = 200,
                                                               use_missing_data = FALSE,
                                                               num_burn_in = 15000,
                                                               num_iterations_after_burn_in = 10000,
                                                               serialize = TRUE)
  
  saveRDS(bartmachine_full_model_6m, "results/Models/bartmachine/bartmachine_full_model_6m.rds")
  
}



plot_pooled_6m.overall <- calibration_plot(
  data = cprd_dataset.6m %>% mutate(pred = bartmachine_full_model_6m$p_hat) %>% mutate(stopdrug_6m_6mFU = as.numeric(stopdrug_6m_6mFU)-1) %>% mutate(grouping = "Pooled"), 
  obs = "stopdrug_6m_6mFU", 
  pred = "pred",
  group = "grouping",
  nTiles = 10)$calibration_plot +
  scale_x_continuous("Prediction", limits = c(0, 0.7), breaks = seq(0, 0.7, 0.1)) +
  scale_y_continuous("Observation", limits = c(0, 0.7), breaks = seq(0, 0.7, 0.1)) +
  theme_bw() +
  scale_colour_manual(values = c("Pooled" = "black", "SGLT2" = "#E69F00", "GLP1" = "#56B4E9", "SU" = "#CC79A7", "DPP4" = "#0072B2", "TZD" = "#D55E00", "MFN" = "grey"), breaks = rev(c("Pooled", "GLP1", "DPP4", "SGLT2", "TZD", "SU")), labels = rev(c("Pooled", "GLP-1RA", "DPP4i", "SGLT2i", "TZD", "SU")), name = "Therapy", guide = guide_legend(reverse = TRUE, nrow = 1)) +
  theme(
    legend.position = "bottom"
  )



plot_pooled_6m.drugs <- calibration_plot(
  data = cprd_dataset.6m %>% mutate(pred = bartmachine_full_model_6m$p_hat) %>% mutate(stopdrug_6m_6mFU = as.numeric(stopdrug_6m_6mFU)-1), 
  obs = "stopdrug_6m_6mFU", 
  pred = "pred", 
  group = "drugclass",
  nTiles = 10)$calibration_plot +
  scale_x_continuous("Prediction", limits = c(0, 0.7), breaks = seq(0, 0.7, 0.1)) +
  scale_y_continuous("Observation", limits = c(0, 0.7), breaks = seq(0, 0.7, 0.1)) +
  theme_bw() +
  scale_colour_manual(values = c("Pooled" = "black", "SGLT2" = "#E69F00", "GLP1" = "#56B4E9", "SU" = "#CC79A7", "DPP4" = "#0072B2", "TZD" = "#D55E00", "MFN" = "grey"), breaks = rev(c("Pooled", "GLP1", "DPP4", "SGLT2", "TZD", "SU")), labels = rev(c("Pooled", "GLP-1RA", "DPP4i", "SGLT2i", "TZD", "SU")), name = "Therapy", guide = guide_legend(reverse = TRUE, nrow = 1)) +
  theme(
    legend.position = "bottom"
  )



# pROC::roc(response = cprd_dataset.6m %>%
#             select(stopdrug_6m_6mFU) %>%
#             mutate(stopdrug_6m_6mFU = factor(stopdrug_6m_6mFU)) %>%
#             unlist(),
#           predictor = bartmachine_full_model_6m$p_hat, ci = TRUE)


# roc_coords_6m <- pROC::ci(pROC::roc(
#   cprd_dataset.6m %>%
#     select(stopdrug_6m_6mFU) %>%
#     mutate(stopdrug_6m_6mFU = factor(stopdrug_6m_6mFU)) %>%
#     unlist(),
#   bartmachine_full_model_6m$p_hat
# ), of = "thresholds")
# 
# saveRDS(roc_coords_6m, "results/Models/bartmachine/roc_coords_6m.rds")

roc_coords_6m <- readRDS("results/Models/bartmachine/roc_coords_6m.rds")


# Run model
if (class(try(
  
  bartmachine_full_model_12m <- readRDS("results/Models/bartmachine/bartmachine_full_model_12m.rds")
  
  , silent = TRUE)) == "try-error") {
  
  bartmachine_full_model_12m <- bartMachine::build_bart_machine(X = cprd_dataset.12m %>%
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
                                                                        prehdl,
                                                                        stopdrug_3m_3mFU_MFN_hist),
                                                               y = cprd_dataset.12m %>% 
                                                                 select(stopdrug_12m_6mFU) %>%
                                                                 mutate(stopdrug_12m_6mFU = factor(stopdrug_12m_6mFU)) %>%
                                                                 unlist(),
                                                               num_trees = 200,
                                                               use_missing_data = FALSE,
                                                               num_burn_in = 15000,
                                                               num_iterations_after_burn_in = 10000,
                                                               serialize = TRUE)
  
  saveRDS(bartmachine_full_model_12m, "results/Models/bartmachine/bartmachine_full_model_12m.rds")
  
}

plot_pooled_12m.overall <- calibration_plot(
  data = cprd_dataset.12m %>% mutate(pred = bartmachine_full_model_12m$p_hat) %>% mutate(stopdrug_12m_6mFU = as.numeric(stopdrug_12m_6mFU)-1) %>% mutate(grouping = "Pooled"), 
  obs = "stopdrug_12m_6mFU", 
  pred = "pred",
  group = "grouping",
  nTiles = 10)$calibration_plot +
  scale_x_continuous("Prediction", limits = c(0, 0.7), breaks = seq(0, 0.7, 0.1)) +
  scale_y_continuous("Observation", limits = c(0, 0.7), breaks = seq(0, 0.7, 0.1)) +
  theme_bw() +
  scale_colour_manual(values = c("Pooled" = "black", "SGLT2" = "#E69F00", "GLP1" = "#56B4E9", "SU" = "#CC79A7", "DPP4" = "#0072B2", "TZD" = "#D55E00", "MFN" = "grey"), breaks = rev(c("Pooled", "GLP1", "DPP4", "SGLT2", "TZD", "SU")), labels = rev(c("Pooled", "GLP-1RA", "DPP4i", "SGLT2i", "TZD", "SU")), name = "Therapy", guide = guide_legend(reverse = TRUE, nrow = 1)) +
  theme(
    legend.position = "bottom"
  )



plot_pooled_12m.drugs <- calibration_plot(
  data = cprd_dataset.12m %>% mutate(pred = bartmachine_full_model_12m$p_hat) %>% mutate(stopdrug_12m_6mFU = as.numeric(stopdrug_12m_6mFU)-1), 
  obs = "stopdrug_12m_6mFU", 
  pred = "pred", 
  group = "drugclass", 
  nTiles = 10)$calibration_plot +
  scale_x_continuous("Prediction", limits = c(0, 0.7), breaks = seq(0, 0.7, 0.1)) +
  scale_y_continuous("Observation", limits = c(0, 0.7), breaks = seq(0, 0.7, 0.1)) +
  theme_bw() +
  scale_colour_manual(values = c("Pooled" = "black", "SGLT2" = "#E69F00", "GLP1" = "#56B4E9", "SU" = "#CC79A7", "DPP4" = "#0072B2", "TZD" = "#D55E00", "MFN" = "grey"), breaks = rev(c("Pooled", "GLP1", "DPP4", "SGLT2", "TZD", "SU")), labels = rev(c("Pooled", "GLP-1RA", "DPP4i", "SGLT2i", "TZD", "SU")), name = "Therapy", guide = guide_legend(reverse = TRUE, nrow = 1)) +
  theme(
    legend.position = "bottom"
  )


pdf("results/figures/07.calibration_overall_and_per_drug.pdf", width = 8, height = 11)

patchwork::wrap_plots(
  
  plot_pooled_3m.overall,
  
  plot_pooled_3m.drugs,
  
  plot_pooled_6m.overall,
  
  plot_pooled_6m.drugs,
  
  plot_pooled_12m.overall,
  
  plot_pooled_12m.drugs,
  
  nrow = 3, ncol = 2
  
) +
  patchwork::plot_layout(guides = 'collect') +
  patchwork::plot_annotation(tag_levels = list(c("A.1", "A.2", "B.1", "B.2", "C.1", "C.2"))) &
  theme(
    legend.position = "bottom",
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 14)
  )

dev.off()



# pROC::roc(response = cprd_dataset.12m %>%
#             select(stopdrug_12m_6mFU) %>%
#             mutate(stopdrug_12m_6mFU = factor(stopdrug_12m_6mFU)) %>%
#             unlist(),
#           predictor = bartmachine_full_model_12m$p_hat, ci = TRUE)



roc_coords_12m <- pROC::ci(pROC::roc(
  cprd_dataset.12m %>%
    select(stopdrug_12m_6mFU) %>%
    mutate(stopdrug_12m_6mFU = factor(stopdrug_12m_6mFU)) %>%
    unlist(),
  bartmachine_full_model_12m$p_hat
), of = "thresholds")

saveRDS(roc_coords_12m, "results/Models/bartmachine/roc_coords_12m.rds")

roc_coords_12m <- readRDS("results/Models/bartmachine/roc_coords_12m.rds")









###############################################################################
###############################################################################
######################## Fit the model: dev and val ###########################
###############################################################################
###############################################################################


# Run model
if (class(try(
  
  bartmachine_full_model <- readRDS("results/Models/bartmachine/bartmachine_full_model.rds")
  
  , silent = TRUE)) == "try-error") {
  
  bartmachine_full_model <- bartMachine::build_bart_machine(X = cprd_dataset.dev %>%
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
                                                       prehdl,
                                                       stopdrug_3m_3mFU_MFN_hist),
                                              y = cprd_dataset.dev %>% 
                                                select(stopdrug_3m_6mFU) %>%
                                                mutate(stopdrug_3m_6mFU = factor(stopdrug_3m_6mFU)) %>%
                                                unlist(),
                                              num_trees = 200,
                                              use_missing_data = FALSE,
                                              num_burn_in = 15000,
                                              num_iterations_after_burn_in = 10000,
                                              serialize = TRUE)
  
  saveRDS(bartmachine_full_model, "results/Models/bartmachine/bartmachine_full_model.rds")
  
}


# Check convergence
# plot_convergence_diagnostics(bartmachine_full_model, plots =  c("sigsqs", "mh_acceptance", "num_nodes", "tree_depths"))



# Predict Development dataset
if (class(try(
  
  bartmachine_pred_dev <- readRDS("results/Models/bartmachine/bartmachine_pred_dev.rds")
  
  , silent = TRUE)) == "try-error") {
  
  bartmachine_pred_dev <- predict(bartmachine_full_model,
                                  cprd_dataset.dev %>%
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
                                                                     prehdl,
                                                                     stopdrug_3m_3mFU_MFN_hist))
  
  saveRDS(bartmachine_pred_dev, "results/Models/bartmachine/bartmachine_pred_dev.rds")
  
}



# Predict Validation dataset
if (class(try(
  
  bartmachine_pred_val <- readRDS("results/Models/bartmachine/bartmachine_pred_val.rds")
  
  , silent = TRUE)) == "try-error") {
  
  bartmachine_pred_val <- predict(bartmachine_full_model,
                                  cprd_dataset.val %>%
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
                                           prehdl,
                                           stopdrug_3m_3mFU_MFN_hist))
  
  saveRDS(bartmachine_pred_val, "results/Models/bartmachine/bartmachine_pred_val.rds")
  
}


 
# pROC::roc(response = cprd_dataset.dev %>%
#             select(stopdrug_3m_6mFU) %>%
#             mutate(stopdrug_3m_6mFU = factor(stopdrug_3m_6mFU)) %>%
#             unlist(),
#           predictor = bartmachine_pred_dev, ci = TRUE)




################################################################################
################################################################################


################################################################################
################################################################################


## Try BART
# Run model
if (class(try(
  
  bart_full_model <- readRDS("results/Models/BART/bart_full_model.rds")
  
  , silent = TRUE)) == "try-error") {
  
  bart_full_model <- wbart(x.train = model.matrix(~., cprd_dataset.dev %>%
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
                                                           prehdl,
                                                           stopdrug_3m_3mFU_MFN_hist)) %>% 
                             as.data.frame() %>% 
                             select(-`(Intercept)`),
                           y.train = cprd_dataset.dev %>%
                             select(stopdrug_3m_6mFU) %>%
                             mutate(stopdrug_3m_6mFU = as.numeric(stopdrug_3m_6mFU) - 1) %>% 
                             unlist(),
                           nskip = 15000,
                           ndpost = 10000)
  
  saveRDS(bart_full_model, "results/Models/BART/bart_full_model.rds")
  
}


if (class(try(
  
  bart_pred_dev <- readRDS("results/Models/BART/bart_pred_dev.rds")
  
  , silent = TRUE)) == "try-error") {
  
  bart_pred_dev = predict(bart_full_model, model.matrix(~., cprd_dataset.dev %>%
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
                                                            prehdl,
                                                            stopdrug_3m_3mFU_MFN_hist)) %>% 
                       as.data.frame() %>% 
                       select(-`(Intercept)`))
  
  saveRDS(bart_pred_dev, "results/Models/BART/bart_pred_dev.rds")
  
}



if (class(try(
  
  bart_pred_val <- readRDS("results/Models/BART/bart_pred_val.rds")
  
  , silent = TRUE)) == "try-error") {
  
  bart_pred_val = predict(bart_full_model, model.matrix(~., cprd_dataset.val %>%
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
                                                            prehdl,
                                                            stopdrug_3m_3mFU_MFN_hist)) %>% 
                       as.data.frame() %>% 
                       select(-`(Intercept)`))
  
  saveRDS(bart_pred_val, "results/Models/BART/bart_pred_val.rds")
  
}




# pROC::roc(response = cprd_dataset.val %>%
#             select(stopdrug_3m_6mFU) %>%
#             mutate(stopdrug_3m_6mFU = factor(stopdrug_3m_6mFU)) %>%
#             unlist(),
#           predictor = colMeans(bart_pred_val), ci = TRUE)






