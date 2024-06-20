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
  therapies = c("MFN", "GLP1", "DPP4", "SGLT2", "TZD", "SU"),
  dataset = "3m.disc.dataset.dev",
  follow_up = "6-months",
  full_prescribing_history = TRUE
) %>%
  drop_na(-stopdrug_6m_6mFU, -stopdrug_12m_6mFU)

# load dataset
cprd_dataset.6m <- set_up_data(
  raw_data = "20240308_t2d_1stinstance",
  diagnosis = FALSE,
  therapies = c("MFN", "GLP1", "DPP4", "SGLT2", "TZD", "SU"),
  dataset = "6m.disc.dataset.dev",
  follow_up = "6-months",
  full_prescribing_history = TRUE
) %>%
  drop_na(-stopdrug_12m_6mFU)

# load dataset
cprd_dataset.12m <- set_up_data(
  raw_data = "20240308_t2d_1stinstance",
  diagnosis = FALSE,
  therapies = c("MFN", "GLP1", "DPP4", "SGLT2", "TZD", "SU"),
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


# Predict Development dataset
if (class(try(
  
  bartmachine_pred_3m <- readRDS("results/Models/bartmachine/bartmachine_pred_3m.rds")
  
  , silent = TRUE)) == "try-error") {
  
  bartmachine_pred_3m <- predict(bartmachine_full_model_3m,
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
                                           prehdl,
                                           stopdrug_3m_3mFU_MFN_hist))
  
  saveRDS(bartmachine_pred_3m, "results/Models/bartmachine/bartmachine_pred_3m.rds")
  
}


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


# Predict Development dataset
if (class(try(
  
  bartmachine_pred_6m <- readRDS("results/Models/bartmachine/bartmachine_pred_6m.rds")
  
  , silent = TRUE)) == "try-error") {
  
  bartmachine_pred_6m <- predict(bartmachine_full_model_6m,
                                 cprd_dataset.6m %>%
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
  
  saveRDS(bartmachine_pred_6m, "results/Models/bartmachine/bartmachine_pred_6m.rds")
  
}



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


# Predict Development dataset
if (class(try(
  
  bartmachine_pred_12m <- readRDS("results/Models/bartmachine/bartmachine_pred_12m.rds")
  
  , silent = TRUE)) == "try-error") {
  
  bartmachine_pred_12m <- predict(bartmachine_full_model_12m,
                                 cprd_dataset.12m %>%
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
  
  saveRDS(bartmachine_pred_12m, "results/Models/bartmachine/bartmachine_pred_12m.rds")
  
}












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






