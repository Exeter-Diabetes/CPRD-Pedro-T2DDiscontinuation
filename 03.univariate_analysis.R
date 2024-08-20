####################
## Description: 
##  - In this file we:
##    - Fit a linear model to each of the variables for pooled therapies and by therapy.
#################### 


# load working directory
setwd("Samples/T2D_Discontinuation")

# load functions
source("code/00.set_up_data_and_functions.R")

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
  raw_data = "20240814_t2d_1stinstance",
  diagnosis = FALSE,
  therapies = c("DPP4", "GLP1", "SGLT2", "SU", "TZD"),
  dataset = "3m.disc.dataset",
  follow_up = "6-months",
  full_prescribing_history = TRUE
) %>%
  drop_na(-prehdl, -stopdrug_6m_6mFU, -stopdrug_12m_6mFU)
  # drop_na(-prehdl, stopdrug_3m_6mFU)


###############################################################################
###############################################################################

## Generating the data needs for the plots
# discontinuation by drug simplified (MFN vs others)

discontinuation_per_drug_combined <- cprd_dataset %>%
  select(stopdrug_3m_6mFU, drugclass) %>%
  mutate(drugclass = factor(drugclass, levels = c("GLP1", "DPP4", "SGLT2", "TZD", "SU"), labels = c("other", "other", "other", "other", "other"))) %>%
  rename("analysis_var" = "stopdrug_3m_6mFU") %>%
  table() %>%
  as.data.frame() %>%
  group_by(drugclass) %>%
  mutate(total_3m = sum(Freq), perc_3m = Freq/total_3m) %>%
  rename("Freq_3m" = "Freq") %>%
  filter(analysis_var == 1) %>%
  left_join(
    cprd_dataset %>%
      drop_na(stopdrug_6m_6mFU) %>%
      select(stopdrug_6m_6mFU, drugclass) %>%
      mutate(drugclass = factor(drugclass, levels = c("GLP1", "DPP4", "SGLT2", "TZD", "SU"), labels = c("other", "other", "other", "other", "other"))) %>%
      rename("analysis_var" = "stopdrug_6m_6mFU") %>%
      table() %>%
      as.data.frame() %>%
      group_by(drugclass) %>%
      mutate(total_6m = sum(Freq), perc_6m = Freq/total_6m) %>%
      rename("Freq_6m" = "Freq") %>%
      filter(analysis_var == 1),
    by = c("drugclass", "analysis_var")
  ) %>%
  left_join(
    cprd_dataset %>%
      drop_na(stopdrug_12m_6mFU) %>%
      select(stopdrug_12m_6mFU, drugclass) %>%
      mutate(drugclass = factor(drugclass, levels = c("GLP1", "DPP4", "SGLT2", "TZD", "SU"), labels = c("other", "other", "other", "other", "other"))) %>%
      rename("analysis_var" = "stopdrug_12m_6mFU") %>%
      table() %>%
      as.data.frame() %>%
      group_by(drugclass) %>%
      mutate(total_12m = sum(Freq), perc_12m = Freq/total_12m) %>%
      rename("Freq_12m" = "Freq") %>%
      filter(analysis_var == 1),
    by = c("drugclass", "analysis_var")
  ) %>%
  select(-Freq_3m, -total_3m, -Freq_6m, -total_6m, -Freq_12m, -total_12m) %>%
  mutate(variable = "overall", group = NA) %>%
  rbind(
    cprd_dataset %>%
      select(stopdrug_3m_3mFU_MFN_hist, stopdrug_3m_6mFU, drugclass) %>%
      mutate(drugclass = factor(drugclass, levels = c("GLP1", "DPP4", "SGLT2", "TZD", "SU"), labels = c("other", "other", "other", "other", "other"))) %>%
      rename("analysis_var" = "stopdrug_3m_6mFU") %>%
      table() %>%
      as.data.frame() %>%
      filter(Freq != 0) %>%
      group_by(stopdrug_3m_3mFU_MFN_hist, drugclass) %>%
      mutate(total_3m = sum(Freq), perc_3m = Freq/total_3m) %>%
      rename("Freq_3m" = "Freq") %>%
      filter(analysis_var == 1) %>%
      left_join(
        cprd_dataset %>%
          drop_na(stopdrug_6m_6mFU) %>%
          select(stopdrug_3m_3mFU_MFN_hist, stopdrug_6m_6mFU, drugclass) %>%
          mutate(drugclass = factor(drugclass, levels = c("GLP1", "DPP4", "SGLT2", "TZD", "SU"), labels = c("other", "other", "other", "other", "other"))) %>%
          rename("analysis_var" = "stopdrug_6m_6mFU") %>%
          table() %>%
          as.data.frame() %>%
          filter(Freq != 0) %>%
          group_by(stopdrug_3m_3mFU_MFN_hist, drugclass) %>%
          mutate(total_6m = sum(Freq), perc_6m = Freq/total_6m) %>%
          rename("Freq_6m" = "Freq") %>%
          filter(analysis_var == 1),
        by = c("analysis_var", "stopdrug_3m_3mFU_MFN_hist", "drugclass")
      ) %>%
      left_join(
        cprd_dataset %>%
          drop_na(stopdrug_12m_6mFU) %>%
          select(stopdrug_3m_3mFU_MFN_hist, stopdrug_12m_6mFU, drugclass) %>%
          mutate(drugclass = factor(drugclass, levels = c("GLP1", "DPP4", "SGLT2", "TZD", "SU"), labels = c("other", "other", "other", "other", "other"))) %>%
          rename("analysis_var" = "stopdrug_12m_6mFU") %>%
          table() %>%
          as.data.frame() %>%
          filter(Freq != 0) %>%
          group_by(stopdrug_3m_3mFU_MFN_hist, drugclass) %>%
          mutate(total_12m = sum(Freq), perc_12m = Freq/total_12m) %>%
          rename("Freq_12m" = "Freq") %>%
          filter(analysis_var == 1),
        by = c("analysis_var", "stopdrug_3m_3mFU_MFN_hist", "drugclass")
      ) %>%
      select(-Freq_3m, -total_3m, -Freq_6m, -total_6m, -Freq_12m, -total_12m) %>%
      mutate(variable = "stopdrug_MFN") %>%
      rename("group" = "stopdrug_3m_3mFU_MFN_hist"),
    cprd_dataset %>%
      select(predrug_statins, stopdrug_3m_6mFU, drugclass) %>%
      mutate(drugclass = factor(drugclass, levels = c("GLP1", "DPP4", "SGLT2", "TZD", "SU"), labels = c("other", "other", "other", "other", "other"))) %>%
      rename("analysis_var" = "stopdrug_3m_6mFU") %>%
      table() %>%
      as.data.frame() %>%
      filter(Freq != 0) %>%
      group_by(predrug_statins, drugclass) %>%
      mutate(total_3m = sum(Freq), perc_3m = Freq/total_3m) %>%
      rename("Freq_3m" = "Freq") %>%
      filter(analysis_var == 1) %>%
      left_join(
        cprd_dataset %>%
          drop_na(stopdrug_6m_6mFU) %>%
          select(predrug_statins, stopdrug_6m_6mFU, drugclass) %>%
          mutate(drugclass = factor(drugclass, levels = c("GLP1", "DPP4", "SGLT2", "TZD", "SU"), labels = c("other", "other", "other", "other", "other"))) %>%
          rename("analysis_var" = "stopdrug_6m_6mFU") %>%
          table() %>%
          as.data.frame() %>%
          filter(Freq != 0) %>%
          group_by(predrug_statins, drugclass) %>%
          mutate(total_6m = sum(Freq), perc_6m = Freq/total_6m) %>%
          rename("Freq_6m" = "Freq") %>%
          filter(analysis_var == 1),
        by = c("analysis_var", "predrug_statins", "drugclass")
      ) %>%
      left_join(
        cprd_dataset %>%
          drop_na(stopdrug_12m_6mFU) %>%
          select(predrug_statins, stopdrug_12m_6mFU, drugclass) %>%
          mutate(drugclass = factor(drugclass, levels = c("GLP1", "DPP4", "SGLT2", "TZD", "SU"), labels = c("other", "other", "other", "other", "other"))) %>%
          rename("analysis_var" = "stopdrug_12m_6mFU") %>%
          table() %>%
          as.data.frame() %>%
          filter(Freq != 0) %>%
          group_by(predrug_statins, drugclass) %>%
          mutate(total_12m = sum(Freq), perc_12m = Freq/total_12m) %>%
          rename("Freq_12m" = "Freq") %>%
          filter(analysis_var == 1),
        by = c("analysis_var", "predrug_statins", "drugclass")
      ) %>%
      select(-Freq_3m, -total_3m, -Freq_6m, -total_6m, -Freq_12m, -total_12m) %>%
      mutate(variable = "statins") %>%
      rename("group" = "predrug_statins"),
    cprd_dataset %>%
      select(predrug_bloodmed, stopdrug_3m_6mFU, drugclass) %>%
      mutate(drugclass = factor(drugclass, levels = c("GLP1", "DPP4", "SGLT2", "TZD", "SU"), labels = c("other", "other", "other", "other", "other"))) %>%
      rename("analysis_var" = "stopdrug_3m_6mFU") %>%
      table() %>%
      as.data.frame() %>%
      filter(Freq != 0) %>%
      group_by(predrug_bloodmed, drugclass) %>%
      mutate(total_3m = sum(Freq), perc_3m = Freq/total_3m) %>%
      rename("Freq_3m" = "Freq") %>%
      filter(analysis_var == 1) %>%
      left_join(
        cprd_dataset %>%
          drop_na(stopdrug_6m_6mFU) %>%
          select(predrug_bloodmed, stopdrug_6m_6mFU, drugclass) %>%
          mutate(drugclass = factor(drugclass, levels = c("GLP1", "DPP4", "SGLT2", "TZD", "SU"), labels = c("other", "other", "other", "other", "other"))) %>%
          rename("analysis_var" = "stopdrug_6m_6mFU") %>%
          table() %>%
          as.data.frame() %>%
          filter(Freq != 0) %>%
          group_by(predrug_bloodmed, drugclass) %>%
          mutate(total_6m = sum(Freq), perc_6m = Freq/total_6m) %>%
          rename("Freq_6m" = "Freq") %>%
          filter(analysis_var == 1),
        by = c("analysis_var", "predrug_bloodmed", "drugclass")
      ) %>%
      left_join(
        cprd_dataset %>%
          drop_na(stopdrug_12m_6mFU) %>%
          select(predrug_bloodmed, stopdrug_12m_6mFU, drugclass) %>%
          mutate(drugclass = factor(drugclass, levels = c("GLP1", "DPP4", "SGLT2", "TZD", "SU"), labels = c("other", "other", "other", "other", "other"))) %>%
          rename("analysis_var" = "stopdrug_12m_6mFU") %>%
          table() %>%
          as.data.frame() %>%
          filter(Freq != 0) %>%
          group_by(predrug_bloodmed, drugclass) %>%
          mutate(total_12m = sum(Freq), perc_12m = Freq/total_12m) %>%
          rename("Freq_12m" = "Freq") %>%
          filter(analysis_var == 1),
        by = c("analysis_var", "predrug_bloodmed", "drugclass")
      ) %>%
      select(-Freq_3m, -total_3m, -Freq_6m, -total_6m, -Freq_12m, -total_12m) %>%
      mutate(variable = "bloodmed") %>%
      rename("group" = "predrug_bloodmed"),
    cprd_dataset %>%
      select(dstartdate_age, stopdrug_3m_6mFU, drugclass) %>%
      mutate(drugclass = factor(drugclass, levels = c("GLP1", "DPP4", "SGLT2", "TZD", "SU"), labels = c("other", "other", "other", "other", "other"))) %>%
      mutate(group = cut(dstartdate_age, breaks = c(18, 50, 60, 70, 150))) %>%
      rename("analysis_var" = "stopdrug_3m_6mFU") %>%
      select(-dstartdate_age) %>%
      table() %>%
      as.data.frame() %>%
      group_by(group, drugclass) %>%
      mutate(total_3m = sum(Freq), perc_3m = Freq/total_3m) %>%
      rename("Freq_3m" = "Freq") %>%
      filter(analysis_var == 1) %>%
      left_join(
        cprd_dataset %>%
          drop_na(stopdrug_6m_6mFU) %>%
          select(dstartdate_age, stopdrug_6m_6mFU, drugclass) %>%
          mutate(drugclass = factor(drugclass, levels = c("GLP1", "DPP4", "SGLT2", "TZD", "SU"), labels = c("other", "other", "other", "other", "other"))) %>%
          mutate(group = cut(dstartdate_age, breaks = c(18, 50, 60, 70, 150))) %>%
          rename("analysis_var" = "stopdrug_6m_6mFU") %>%
          select(-dstartdate_age) %>%
          table() %>%
          as.data.frame() %>%
          group_by(group, drugclass) %>%
          mutate(total_6m = sum(Freq), perc_6m = Freq/total_6m) %>%
          rename("Freq_6m" = "Freq") %>%
          filter(analysis_var == 1),
        by = c("drugclass", "group", "analysis_var")
      ) %>%
      left_join(
        cprd_dataset %>%
          drop_na(stopdrug_12m_6mFU) %>%
          select(dstartdate_age, stopdrug_12m_6mFU, drugclass) %>%
          mutate(drugclass = factor(drugclass, levels = c("GLP1", "DPP4", "SGLT2", "TZD", "SU"), labels = c("other", "other", "other", "other", "other"))) %>%
          mutate(group = cut(dstartdate_age, breaks = c(18, 50, 60, 70, 150))) %>%
          rename("analysis_var" = "stopdrug_12m_6mFU") %>%
          select(-dstartdate_age) %>%
          table() %>%
          as.data.frame() %>%
          group_by(group, drugclass) %>%
          mutate(total_12m = sum(Freq), perc_12m = Freq/total_12m) %>%
          rename("Freq_12m" = "Freq") %>%
          filter(analysis_var == 1),
        by = c("drugclass", "group", "analysis_var")
      ) %>%
      select(-Freq_3m, -total_3m, -Freq_6m, -total_6m, -Freq_12m, -total_12m) %>%
      mutate(variable = "dstartdate_age"),
    cprd_dataset %>%
      select(gender, stopdrug_3m_6mFU, drugclass) %>%
      mutate(drugclass = factor(drugclass, levels = c("GLP1", "DPP4", "SGLT2", "TZD", "SU"), labels = c("other", "other", "other", "other", "other"))) %>%
      rename("analysis_var" = "stopdrug_3m_6mFU") %>%
      table() %>%
      as.data.frame() %>%
      group_by(gender, drugclass) %>%
      mutate(total_3m = sum(Freq), perc_3m = Freq/total_3m) %>%
      rename("Freq_3m" = "Freq") %>%
      filter(analysis_var == 1) %>%
      left_join(
        cprd_dataset %>%
          drop_na(stopdrug_6m_6mFU) %>%
          select(gender, stopdrug_6m_6mFU, drugclass) %>%
          mutate(drugclass = factor(drugclass, levels = c("GLP1", "DPP4", "SGLT2", "TZD", "SU"), labels = c("other", "other", "other", "other", "other"))) %>%
          rename("analysis_var" = "stopdrug_6m_6mFU") %>%
          table() %>%
          as.data.frame() %>%
          group_by(gender, drugclass) %>%
          mutate(total_6m = sum(Freq), perc_6m = Freq/total_6m) %>%
          rename("Freq_6m" = "Freq") %>%
          filter(analysis_var == 1),
        by = c("analysis_var", "gender", "drugclass")
      ) %>%
      left_join(
        cprd_dataset %>%
          drop_na(stopdrug_12m_6mFU) %>%
          select(gender, stopdrug_12m_6mFU, drugclass) %>%
          mutate(drugclass = factor(drugclass, levels = c("GLP1", "DPP4", "SGLT2", "TZD", "SU"), labels = c("other", "other", "other", "other", "other"))) %>%
          rename("analysis_var" = "stopdrug_12m_6mFU") %>%
          table() %>%
          as.data.frame() %>%
          group_by(gender, drugclass) %>%
          mutate(total_12m = sum(Freq), perc_12m = Freq/total_12m) %>%
          rename("Freq_12m" = "Freq") %>%
          filter(analysis_var == 1),
        by = c("analysis_var", "gender", "drugclass")
      ) %>%
      select(-Freq_3m, -total_3m, -Freq_6m, -total_6m, -Freq_12m, -total_12m) %>%
      mutate(variable = "gender") %>%
      rename("group" = "gender"),
    cprd_dataset %>%
      select(imd2015_10, stopdrug_3m_6mFU, drugclass) %>%
      mutate(drugclass = factor(drugclass, levels = c("GLP1", "DPP4", "SGLT2", "TZD", "SU"), labels = c("other", "other", "other", "other", "other"))) %>%
      rename("analysis_var" = "stopdrug_3m_6mFU") %>%
      table() %>%
      as.data.frame() %>%
      group_by(imd2015_10, drugclass) %>%
      mutate(total_3m = sum(Freq), perc_3m = Freq/total_3m) %>%
      rename("Freq_3m" = "Freq") %>%
      filter(analysis_var == 1) %>%
      left_join(
        cprd_dataset %>%
          drop_na(stopdrug_6m_6mFU) %>%
          select(imd2015_10, stopdrug_6m_6mFU, drugclass) %>%
          mutate(drugclass = factor(drugclass, levels = c("GLP1", "DPP4", "SGLT2", "TZD", "SU"), labels = c("other", "other", "other", "other", "other"))) %>%
          rename("analysis_var" = "stopdrug_6m_6mFU") %>%
          table() %>%
          as.data.frame() %>%
          group_by(imd2015_10, drugclass) %>%
          mutate(total_6m = sum(Freq), perc_6m = Freq/total_6m) %>%
          rename("Freq_6m" = "Freq") %>%
          filter(analysis_var == 1),
        by = c("analysis_var", "imd2015_10", "drugclass")
      ) %>%
      left_join(
        cprd_dataset %>%
          drop_na(stopdrug_12m_6mFU) %>%
          select(imd2015_10, stopdrug_12m_6mFU, drugclass) %>%
          mutate(drugclass = factor(drugclass, levels = c("GLP1", "DPP4", "SGLT2", "TZD", "SU"), labels = c("other", "other", "other", "other", "other"))) %>%
          rename("analysis_var" = "stopdrug_12m_6mFU") %>%
          table() %>%
          as.data.frame() %>%
          group_by(imd2015_10, drugclass) %>%
          mutate(total_12m = sum(Freq), perc_12m = Freq/total_12m) %>%
          rename("Freq_12m" = "Freq") %>%
          filter(analysis_var == 1),
        by = c("analysis_var", "imd2015_10", "drugclass")
      ) %>%
      select(-Freq_3m, -total_3m, -Freq_6m, -total_6m, -Freq_12m, -total_12m) %>%
      mutate(variable = "imd2015") %>%
      rename("group" = "imd2015_10"),
    cprd_dataset %>%
      select(prebmi, stopdrug_3m_6mFU, drugclass) %>%
      mutate(drugclass = factor(drugclass, levels = c("GLP1", "DPP4", "SGLT2", "TZD", "SU"), labels = c("other", "other", "other", "other", "other"))) %>%
      mutate(group = cut(prebmi, breaks = c(1, 25, 30, 35, 500))) %>%
      rename("analysis_var" = "stopdrug_3m_6mFU") %>%
      select(-prebmi) %>%
      table() %>%
      as.data.frame() %>%
      group_by(group, drugclass) %>%
      mutate(total_3m = sum(Freq), perc_3m = Freq/total_3m) %>%
      rename("Freq_3m" = "Freq") %>%
      filter(analysis_var == 1) %>%
      left_join(
        cprd_dataset %>%
          drop_na(stopdrug_6m_6mFU) %>%
          select(prebmi, stopdrug_6m_6mFU, drugclass) %>%
          mutate(drugclass = factor(drugclass, levels = c("GLP1", "DPP4", "SGLT2", "TZD", "SU"), labels = c("other", "other", "other", "other", "other"))) %>%
          mutate(group = cut(prebmi, breaks = c(1, 25, 30, 35, 500))) %>%
          rename("analysis_var" = "stopdrug_6m_6mFU") %>%
          select(-prebmi) %>%
          table() %>%
          as.data.frame() %>%
          group_by(group, drugclass) %>%
          mutate(total_6m = sum(Freq), perc_6m = Freq/total_6m) %>%
          rename("Freq_6m" = "Freq") %>%
          filter(analysis_var == 1),
        by = c("drugclass", "group", "analysis_var")
      ) %>%
      left_join(
        cprd_dataset %>%
          drop_na(stopdrug_12m_6mFU) %>%
          select(prebmi, stopdrug_12m_6mFU, drugclass) %>%
          mutate(drugclass = factor(drugclass, levels = c("GLP1", "DPP4", "SGLT2", "TZD", "SU"), labels = c("other", "other", "other", "other", "other"))) %>%
          mutate(group = cut(prebmi, breaks = c(1, 25, 30, 35, 500))) %>%
          rename("analysis_var" = "stopdrug_12m_6mFU") %>%
          select(-prebmi) %>%
          table() %>%
          as.data.frame() %>%
          group_by(group, drugclass) %>%
          mutate(total_12m = sum(Freq), perc_12m = Freq/total_12m) %>%
          rename("Freq_12m" = "Freq") %>%
          filter(analysis_var == 1),
        by = c("drugclass", "group", "analysis_var")
      ) %>%
      select(-Freq_3m, -total_3m, -Freq_6m, -total_6m, -Freq_12m, -total_12m) %>%
      mutate(variable = "bmi"),
    cprd_dataset %>%
      select(dstartdate_dm_dur, stopdrug_3m_6mFU, drugclass) %>%
      mutate(drugclass = factor(drugclass, levels = c("GLP1", "DPP4", "SGLT2", "TZD", "SU"), labels = c("other", "other", "other", "other", "other"))) %>%
      mutate(group = cut(dstartdate_dm_dur, breaks = c(-1, 1, 3, 6, 10, 150))) %>%
      rename("analysis_var" = "stopdrug_3m_6mFU") %>%
      select(-dstartdate_dm_dur) %>%
      table() %>%
      as.data.frame() %>%
      group_by(group, drugclass) %>%
      mutate(total_3m = sum(Freq), perc_3m = Freq/total_3m) %>%
      rename("Freq_3m" = "Freq") %>%
      filter(analysis_var == 1) %>%
      left_join(
        cprd_dataset %>%
          drop_na(stopdrug_6m_6mFU) %>%
          select(dstartdate_dm_dur, stopdrug_6m_6mFU, drugclass) %>%
          mutate(drugclass = factor(drugclass, levels = c("GLP1", "DPP4", "SGLT2", "TZD", "SU"), labels = c("other", "other", "other", "other", "other"))) %>%
          mutate(group = cut(dstartdate_dm_dur, breaks = c(-1, 1, 3, 6, 10, 150))) %>%
          rename("analysis_var" = "stopdrug_6m_6mFU") %>%
          select(-dstartdate_dm_dur) %>%
          table() %>%
          as.data.frame() %>%
          group_by(group, drugclass) %>%
          mutate(total_6m = sum(Freq), perc_6m = Freq/total_6m) %>%
          rename("Freq_6m" = "Freq") %>%
          filter(analysis_var == 1),
        by = c("drugclass", "group", "analysis_var")
      ) %>%
      left_join(
        cprd_dataset %>%
          drop_na(stopdrug_12m_6mFU) %>%
          select(dstartdate_dm_dur, stopdrug_12m_6mFU, drugclass) %>%
          mutate(drugclass = factor(drugclass, levels = c("GLP1", "DPP4", "SGLT2", "TZD", "SU"), labels = c("other", "other", "other", "other", "other"))) %>%
          mutate(group = cut(dstartdate_dm_dur, breaks = c(-1, 1, 3, 6, 10, 150))) %>%
          rename("analysis_var" = "stopdrug_12m_6mFU") %>%
          select(-dstartdate_dm_dur) %>%
          table() %>%
          as.data.frame() %>%
          group_by(group, drugclass) %>%
          mutate(total_12m = sum(Freq), perc_12m = Freq/total_12m) %>%
          rename("Freq_12m" = "Freq") %>%
          filter(analysis_var == 1),
        by = c("drugclass", "group", "analysis_var")
      ) %>%
      select(-Freq_3m, -total_3m, -Freq_6m, -total_6m, -Freq_12m, -total_12m) %>%
      mutate(variable = "dstartdate_dm_dur"),
    cprd_dataset %>%
      select(prehba1c, stopdrug_3m_6mFU, drugclass) %>%
      mutate(drugclass = factor(drugclass, levels = c("GLP1", "DPP4", "SGLT2", "TZD", "SU"), labels = c("other", "other", "other", "other", "other"))) %>%
      mutate(group = cut(prehba1c, breaks = c(52, 64, 75, 86, 250))) %>%
      rename("analysis_var" = "stopdrug_3m_6mFU") %>%
      select(-prehba1c) %>%
      table() %>%
      as.data.frame() %>%
      group_by(group, drugclass) %>%
      mutate(total_3m = sum(Freq), perc_3m = Freq/total_3m) %>%
      rename("Freq_3m" = "Freq") %>%
      filter(analysis_var == 1) %>%
      left_join(
        cprd_dataset %>%
          drop_na(stopdrug_6m_6mFU) %>%
          select(prehba1c, stopdrug_6m_6mFU, drugclass) %>%
          mutate(drugclass = factor(drugclass, levels = c("GLP1", "DPP4", "SGLT2", "TZD", "SU"), labels = c("other", "other", "other", "other", "other"))) %>%
          mutate(group = cut(prehba1c, breaks = c(52, 64, 75, 86, 250))) %>%
          rename("analysis_var" = "stopdrug_6m_6mFU") %>%
          select(-prehba1c) %>%
          table() %>%
          as.data.frame() %>%
          group_by(group, drugclass) %>%
          mutate(total_6m = sum(Freq), perc_6m = Freq/total_6m) %>%
          rename("Freq_6m" = "Freq") %>%
          filter(analysis_var == 1),
        by = c("drugclass", "group", "analysis_var")
      ) %>%
      left_join(
        cprd_dataset %>%
          drop_na(stopdrug_12m_6mFU) %>%
          select(prehba1c, stopdrug_12m_6mFU, drugclass) %>%
          mutate(drugclass = factor(drugclass, levels = c("GLP1", "DPP4", "SGLT2", "TZD", "SU"), labels = c("other", "other", "other", "other", "other"))) %>%
          mutate(group = cut(prehba1c, breaks = c(52, 64, 75, 86, 250))) %>%
          rename("analysis_var" = "stopdrug_12m_6mFU") %>%
          select(-prehba1c) %>%
          table() %>%
          as.data.frame() %>%
          group_by(group, drugclass) %>%
          mutate(total_12m = sum(Freq), perc_12m = Freq/total_12m) %>%
          rename("Freq_12m" = "Freq") %>%
          filter(analysis_var == 1),
        by = c("drugclass", "group", "analysis_var")
      ) %>%
      select(-Freq_3m, -total_3m, -Freq_6m, -total_6m, -Freq_12m, -total_12m) %>%
      mutate(variable = "hba1c"),
    cprd_dataset %>%
      select(drugline, stopdrug_3m_6mFU, drugclass) %>%
      mutate(drugline = factor(drugline, levels = c("2", "3", "4", "5+"), labels = c("2", "3", "4+", "4+"))) %>%
      mutate(drugclass = factor(drugclass, levels = c("GLP1", "DPP4", "SGLT2", "TZD", "SU"), labels = c("other", "other", "other", "other", "other"))) %>%
      rename("analysis_var" = "stopdrug_3m_6mFU") %>%
      table() %>%
      as.data.frame() %>%
      group_by(drugline, drugclass) %>%
      mutate(total_3m = sum(Freq), perc_3m = Freq/total_3m) %>%
      rename("Freq_3m" = "Freq") %>%
      filter(analysis_var == 1) %>%
      left_join(
        cprd_dataset %>%
          drop_na(stopdrug_6m_6mFU) %>%
          select(drugline, stopdrug_6m_6mFU, drugclass) %>%
          mutate(drugline = factor(drugline, levels = c("2", "3", "4", "5+"), labels = c("2", "3", "4+", "4+"))) %>%
          mutate(drugclass = factor(drugclass, levels = c("GLP1", "DPP4", "SGLT2", "TZD", "SU"), labels = c("other", "other", "other", "other", "other"))) %>%
          rename("analysis_var" = "stopdrug_6m_6mFU") %>%
          table() %>%
          as.data.frame() %>%
          group_by(drugline, drugclass) %>%
          mutate(total_6m = sum(Freq), perc_6m = Freq/total_6m) %>%
          rename("Freq_6m" = "Freq") %>%
          filter(analysis_var == 1),
        by = c("analysis_var", "drugline", "drugclass")
      ) %>%
      left_join(
        cprd_dataset %>%
          drop_na(stopdrug_12m_6mFU) %>%
          select(drugline, stopdrug_12m_6mFU, drugclass) %>%
          mutate(drugline = factor(drugline, levels = c("2", "3", "4", "5+"), labels = c("2", "3", "4+", "4+"))) %>%
          mutate(drugclass = factor(drugclass, levels = c("GLP1", "DPP4", "SGLT2", "TZD", "SU"), labels = c("other", "other", "other", "other", "other"))) %>%
          rename("analysis_var" = "stopdrug_12m_6mFU") %>%
          table() %>%
          as.data.frame() %>%
          group_by(drugline, drugclass) %>%
          mutate(total_12m = sum(Freq), perc_12m = Freq/total_12m) %>%
          rename("Freq_12m" = "Freq") %>%
          filter(analysis_var == 1),
        by = c("analysis_var", "drugline", "drugclass")
      ) %>%
      select(-Freq_3m, -total_3m, -Freq_6m, -total_6m, -Freq_12m, -total_12m) %>%
      mutate(variable = "drugline") %>%
      rename("group" = "drugline"),
    cprd_dataset %>%
      select(predrug_frailty_proxy, stopdrug_3m_6mFU, drugclass) %>%
      mutate(drugclass = factor(drugclass, levels = c("GLP1", "DPP4", "SGLT2", "TZD", "SU"), labels = c("other", "other", "other", "other", "other"))) %>%
      rename("analysis_var" = "stopdrug_3m_6mFU") %>%
      table() %>%
      as.data.frame() %>%
      group_by(predrug_frailty_proxy, drugclass) %>%
      mutate(total_3m = sum(Freq), perc_3m = Freq/total_3m) %>%
      rename("Freq_3m" = "Freq") %>%
      filter(analysis_var == 1) %>%
      left_join(
        cprd_dataset %>%
          drop_na(stopdrug_6m_6mFU) %>%
          select(predrug_frailty_proxy, stopdrug_6m_6mFU, drugclass) %>%
          mutate(drugclass = factor(drugclass, levels = c("GLP1", "DPP4", "SGLT2", "TZD", "SU"), labels = c("other", "other", "other", "other", "other"))) %>%
          rename("analysis_var" = "stopdrug_6m_6mFU") %>%
          table() %>%
          as.data.frame() %>%
          group_by(predrug_frailty_proxy, drugclass) %>%
          mutate(total_6m = sum(Freq), perc_6m = Freq/total_6m) %>%
          rename("Freq_6m" = "Freq") %>%
          filter(analysis_var == 1),
        by = c("analysis_var", "predrug_frailty_proxy", "drugclass")
      ) %>%
      left_join(
        cprd_dataset %>%
          drop_na(stopdrug_12m_6mFU) %>%
          select(predrug_frailty_proxy, stopdrug_12m_6mFU, drugclass) %>%
          mutate(drugclass = factor(drugclass, levels = c("GLP1", "DPP4", "SGLT2", "TZD", "SU"), labels = c("other", "other", "other", "other", "other"))) %>%
          rename("analysis_var" = "stopdrug_12m_6mFU") %>%
          table() %>%
          as.data.frame() %>%
          group_by(predrug_frailty_proxy, drugclass) %>%
          mutate(total_12m = sum(Freq), perc_12m = Freq/total_12m) %>%
          rename("Freq_12m" = "Freq") %>%
          filter(analysis_var == 1),
        by = c("analysis_var", "predrug_frailty_proxy", "drugclass")
      ) %>%
      select(-Freq_3m, -total_3m, -Freq_6m, -total_6m, -Freq_12m, -total_12m) %>%
      mutate(variable = "frailty") %>%
      rename("group" = "predrug_frailty_proxy"),
    cprd_dataset %>%
      select(preegfr, stopdrug_3m_6mFU, drugclass) %>%
      mutate(drugclass = factor(drugclass, levels = c("GLP1", "DPP4", "SGLT2", "TZD", "SU"), labels = c("other", "other", "other", "other", "other"))) %>%
      mutate(group = cut(preegfr, breaks = c(0, 60, 90, 500))) %>%
      rename("analysis_var" = "stopdrug_3m_6mFU") %>%
      select(-preegfr) %>%
      table() %>%
      as.data.frame() %>%
      group_by(group, drugclass) %>%
      mutate(total_3m = sum(Freq), perc_3m = Freq/total_3m) %>%
      rename("Freq_3m" = "Freq") %>%
      filter(analysis_var == 1) %>%
      left_join(
        cprd_dataset %>%
          drop_na(stopdrug_6m_6mFU) %>%
          select(preegfr, stopdrug_6m_6mFU, drugclass) %>%
          mutate(drugclass = factor(drugclass, levels = c("GLP1", "DPP4", "SGLT2", "TZD", "SU"), labels = c("other", "other", "other", "other", "other"))) %>%
          mutate(group = cut(preegfr, breaks = c(0, 60, 90, 500))) %>%
          rename("analysis_var" = "stopdrug_6m_6mFU") %>%
          select(-preegfr) %>%
          table() %>%
          as.data.frame() %>%
          group_by(group, drugclass) %>%
          mutate(total_6m = sum(Freq), perc_6m = Freq/total_6m) %>%
          rename("Freq_6m" = "Freq") %>%
          filter(analysis_var == 1),
        by = c("drugclass", "group", "analysis_var")
      ) %>%
      left_join(
        cprd_dataset %>%
          drop_na(stopdrug_12m_6mFU) %>%
          select(preegfr, stopdrug_12m_6mFU, drugclass) %>%
          mutate(drugclass = factor(drugclass, levels = c("GLP1", "DPP4", "SGLT2", "TZD", "SU"), labels = c("other", "other", "other", "other", "other"))) %>%
          mutate(group = cut(preegfr, breaks = c(0, 60, 90, 500))) %>%
          rename("analysis_var" = "stopdrug_12m_6mFU") %>%
          select(-preegfr) %>%
          table() %>%
          as.data.frame() %>%
          group_by(group, drugclass) %>%
          mutate(total_12m = sum(Freq), perc_12m = Freq/total_12m) %>%
          rename("Freq_12m" = "Freq") %>%
          filter(analysis_var == 1),
        by = c("drugclass", "group", "analysis_var")
      ) %>%
      select(-Freq_3m, -total_3m, -Freq_6m, -total_6m, -Freq_12m, -total_12m) %>%
      mutate(variable = "preegfr"),
    cprd_dataset %>%
      select(smoking_cat, stopdrug_3m_6mFU, drugclass) %>%
      mutate(drugclass = factor(drugclass, levels = c("GLP1", "DPP4", "SGLT2", "TZD", "SU"), labels = c("other", "other", "other", "other", "other"))) %>%
      rename("analysis_var" = "stopdrug_3m_6mFU") %>%
      table() %>%
      as.data.frame() %>%
      group_by(smoking_cat, drugclass) %>%
      mutate(total_3m = sum(Freq), perc_3m = Freq/total_3m) %>%
      rename("Freq_3m" = "Freq") %>%
      filter(analysis_var == 1) %>%
      left_join(
        cprd_dataset %>%
          drop_na(stopdrug_6m_6mFU) %>%
          select(smoking_cat, stopdrug_6m_6mFU, drugclass) %>%
          mutate(drugclass = factor(drugclass, levels = c("GLP1", "DPP4", "SGLT2", "TZD", "SU"), labels = c("other", "other", "other", "other", "other"))) %>%
          rename("analysis_var" = "stopdrug_6m_6mFU") %>%
          table() %>%
          as.data.frame() %>%
          group_by(smoking_cat, drugclass) %>%
          mutate(total_6m = sum(Freq), perc_6m = Freq/total_6m) %>%
          rename("Freq_6m" = "Freq") %>%
          filter(analysis_var == 1),
        by = c("analysis_var", "smoking_cat", "drugclass")
      ) %>%
      left_join(
        cprd_dataset %>%
          drop_na(stopdrug_12m_6mFU) %>%
          select(smoking_cat, stopdrug_12m_6mFU, drugclass) %>%
          mutate(drugclass = factor(drugclass, levels = c("GLP1", "DPP4", "SGLT2", "TZD", "SU"), labels = c("other", "other", "other", "other", "other"))) %>%
          rename("analysis_var" = "stopdrug_12m_6mFU") %>%
          table() %>%
          as.data.frame() %>%
          group_by(smoking_cat, drugclass) %>%
          mutate(total_12m = sum(Freq), perc_12m = Freq/total_12m) %>%
          rename("Freq_12m" = "Freq") %>%
          filter(analysis_var == 1),
        by = c("analysis_var", "smoking_cat", "drugclass")
      ) %>%
      select(-Freq_3m, -total_3m, -Freq_6m, -total_6m, -Freq_12m, -total_12m) %>%
      mutate(variable = "smoking") %>%
      rename("group" = "smoking_cat"),
    cprd_dataset %>%
      select(numdrugs, stopdrug_3m_6mFU, drugclass) %>%
      mutate(drugclass = factor(drugclass, levels = c("GLP1", "DPP4", "SGLT2", "TZD", "SU"), labels = c("other", "other", "other", "other", "other"))) %>%
      rename("analysis_var" = "stopdrug_3m_6mFU") %>%
      table() %>%
      as.data.frame() %>%
      group_by(numdrugs, drugclass) %>%
      mutate(total_3m = sum(Freq), perc_3m = Freq/total_3m) %>%
      rename("Freq_3m" = "Freq") %>%
      filter(analysis_var == 1) %>%
      left_join(
        cprd_dataset %>%
          drop_na(stopdrug_6m_6mFU) %>%
          select(numdrugs, stopdrug_6m_6mFU, drugclass) %>%
          mutate(drugclass = factor(drugclass, levels = c("GLP1", "DPP4", "SGLT2", "TZD", "SU"), labels = c("other", "other", "other", "other", "other"))) %>%
          rename("analysis_var" = "stopdrug_6m_6mFU") %>%
          table() %>%
          as.data.frame() %>%
          group_by(numdrugs, drugclass) %>%
          mutate(total_6m = sum(Freq), perc_6m = Freq/total_6m) %>%
          rename("Freq_6m" = "Freq") %>%
          filter(analysis_var == 1),
        by = c("analysis_var", "numdrugs", "drugclass")
      ) %>%
      left_join(
        cprd_dataset %>%
          drop_na(stopdrug_12m_6mFU) %>%
          select(numdrugs, stopdrug_12m_6mFU, drugclass) %>%
          mutate(drugclass = factor(drugclass, levels = c("GLP1", "DPP4", "SGLT2", "TZD", "SU"), labels = c("other", "other", "other", "other", "other"))) %>%
          rename("analysis_var" = "stopdrug_12m_6mFU") %>%
          table() %>%
          as.data.frame() %>%
          group_by(numdrugs, drugclass) %>%
          mutate(total_12m = sum(Freq), perc_12m = Freq/total_12m) %>%
          rename("Freq_12m" = "Freq") %>%
          filter(analysis_var == 1),
        by = c("analysis_var", "numdrugs", "drugclass")
      ) %>%
      select(-Freq_3m, -total_3m, -Freq_6m, -total_6m, -Freq_12m, -total_12m) %>%
      mutate(variable = "numdrugs") %>%
      rename("group" = "numdrugs"),
    cprd_dataset %>%
      select(ethnicity_5cat, stopdrug_3m_6mFU, drugclass) %>%
      mutate(ethnicity_5cat = factor(ethnicity_5cat, levels = c("4", "3", "2", "1", "0"), labels = c("Mixed & Other", "Mixed & Other", "Black", "South Asian", "White"))) %>%
      mutate(drugclass = factor(drugclass, levels = c("GLP1", "DPP4", "SGLT2", "TZD", "SU"), labels = c("other", "other", "other", "other", "other"))) %>%
      rename("analysis_var" = "stopdrug_3m_6mFU") %>%
      table() %>%
      as.data.frame() %>%
      group_by(ethnicity_5cat, drugclass) %>%
      mutate(total_3m = sum(Freq), perc_3m = Freq/total_3m) %>%
      rename("Freq_3m" = "Freq") %>%
      filter(analysis_var == 1) %>%
      left_join(
        cprd_dataset %>%
          drop_na(stopdrug_6m_6mFU) %>%
          select(ethnicity_5cat, stopdrug_6m_6mFU, drugclass) %>%
          mutate(ethnicity_5cat = factor(ethnicity_5cat, levels = c("4", "3", "2", "1", "0"), labels = c("Mixed & Other", "Mixed & Other", "Black", "South Asian", "White"))) %>%
          mutate(drugclass = factor(drugclass, levels = c("GLP1", "DPP4", "SGLT2", "TZD", "SU"), labels = c("other", "other", "other", "other", "other"))) %>%
          rename("analysis_var" = "stopdrug_6m_6mFU") %>%
          table() %>%
          as.data.frame() %>%
          group_by(ethnicity_5cat, drugclass) %>%
          mutate(total_6m = sum(Freq), perc_6m = Freq/total_6m) %>%
          rename("Freq_6m" = "Freq") %>%
          filter(analysis_var == 1),
        by = c("analysis_var", "ethnicity_5cat", "drugclass")
      ) %>%
      left_join(
        cprd_dataset %>%
          drop_na(stopdrug_12m_6mFU) %>%
          select(ethnicity_5cat, stopdrug_12m_6mFU, drugclass) %>%
          mutate(ethnicity_5cat = factor(ethnicity_5cat, levels = c("4", "3", "2", "1", "0"), labels = c("Mixed & Other", "Mixed & Other", "Black", "South Asian", "White"))) %>%
          mutate(drugclass = factor(drugclass, levels = c("GLP1", "DPP4", "SGLT2", "TZD", "SU"), labels = c("other", "other", "other", "other", "other"))) %>%
          rename("analysis_var" = "stopdrug_12m_6mFU") %>%
          table() %>%
          as.data.frame() %>%
          group_by(ethnicity_5cat, drugclass) %>%
          mutate(total_12m = sum(Freq), perc_12m = Freq/total_12m) %>%
          rename("Freq_12m" = "Freq") %>%
          filter(analysis_var == 1),
        by = c("analysis_var", "ethnicity_5cat", "drugclass")
      ) %>%
      select(-Freq_3m, -total_3m, -Freq_6m, -total_6m, -Freq_12m, -total_12m) %>%
      mutate(variable = "ethnicity") %>%
      rename("group" = "ethnicity_5cat")
  )




discontinuation_per_drug_not_combined <- cprd_dataset %>%
  select(stopdrug_3m_6mFU, drugclass) %>%
  mutate(drugclass = factor(drugclass, levels = c("GLP1", "DPP4", "SGLT2", "TZD", "SU"))) %>%
  rename("analysis_var" = "stopdrug_3m_6mFU") %>%
  table() %>%
  as.data.frame() %>%
  group_by(drugclass) %>%
  mutate(total_3m = sum(Freq), perc_3m = Freq/total_3m) %>%
  rename("Freq_3m" = "Freq") %>%
  filter(analysis_var == 1) %>%
  left_join(
    cprd_dataset %>%
      drop_na(stopdrug_6m_6mFU) %>%
      select(stopdrug_6m_6mFU, drugclass) %>%
      mutate(drugclass = factor(drugclass, levels = c("GLP1", "DPP4", "SGLT2", "TZD", "SU"))) %>%
      rename("analysis_var" = "stopdrug_6m_6mFU") %>%
      table() %>%
      as.data.frame() %>%
      group_by(drugclass) %>%
      mutate(total_6m = sum(Freq), perc_6m = Freq/total_6m) %>%
      rename("Freq_6m" = "Freq") %>%
      filter(analysis_var == 1),
    by = c("drugclass", "analysis_var")
  ) %>%
  left_join(
    cprd_dataset %>%
      drop_na(stopdrug_12m_6mFU) %>%
      select(stopdrug_12m_6mFU, drugclass) %>%
      mutate(drugclass = factor(drugclass, levels = c("GLP1", "DPP4", "SGLT2", "TZD", "SU"))) %>%
      rename("analysis_var" = "stopdrug_12m_6mFU") %>%
      table() %>%
      as.data.frame() %>%
      group_by(drugclass) %>%
      mutate(total_12m = sum(Freq), perc_12m = Freq/total_12m) %>%
      rename("Freq_12m" = "Freq") %>%
      filter(analysis_var == 1),
    by = c("drugclass", "analysis_var")
  ) %>%
  select(-Freq_3m, -total_3m, -Freq_6m, -total_6m, -Freq_12m, -total_12m) %>%
  mutate(variable = "overall", group = NA) %>%
  rbind(
    cprd_dataset %>%
      select(stopdrug_3m_3mFU_MFN_hist, stopdrug_3m_6mFU, drugclass) %>%
      mutate(drugclass = factor(drugclass, levels = c("GLP1", "DPP4", "SGLT2", "TZD", "SU"))) %>%
      rename("analysis_var" = "stopdrug_3m_6mFU") %>%
      table() %>%
      as.data.frame() %>%
      filter(Freq != 0) %>%
      group_by(stopdrug_3m_3mFU_MFN_hist, drugclass) %>%
      mutate(total_3m = sum(Freq), perc_3m = Freq/total_3m) %>%
      rename("Freq_3m" = "Freq") %>%
      filter(analysis_var == 1) %>%
      left_join(
        cprd_dataset %>%
          drop_na(stopdrug_6m_6mFU) %>%
          select(stopdrug_3m_3mFU_MFN_hist, stopdrug_6m_6mFU, drugclass) %>%
          mutate(drugclass = factor(drugclass, levels = c("GLP1", "DPP4", "SGLT2", "TZD", "SU"))) %>%
          rename("analysis_var" = "stopdrug_6m_6mFU") %>%
          table() %>%
          as.data.frame() %>%
          filter(Freq != 0) %>%
          group_by(stopdrug_3m_3mFU_MFN_hist, drugclass) %>%
          mutate(total_6m = sum(Freq), perc_6m = Freq/total_6m) %>%
          rename("Freq_6m" = "Freq") %>%
          filter(analysis_var == 1),
        by = c("analysis_var", "stopdrug_3m_3mFU_MFN_hist", "drugclass")
      ) %>%
      left_join(
        cprd_dataset %>%
          drop_na(stopdrug_12m_6mFU) %>%
          select(stopdrug_3m_3mFU_MFN_hist, stopdrug_12m_6mFU, drugclass) %>%
          mutate(drugclass = factor(drugclass, levels = c("GLP1", "DPP4", "SGLT2", "TZD", "SU"))) %>%
          rename("analysis_var" = "stopdrug_12m_6mFU") %>%
          table() %>%
          as.data.frame() %>%
          filter(Freq != 0) %>%
          group_by(stopdrug_3m_3mFU_MFN_hist, drugclass) %>%
          mutate(total_12m = sum(Freq), perc_12m = Freq/total_12m) %>%
          rename("Freq_12m" = "Freq") %>%
          filter(analysis_var == 1),
        by = c("analysis_var", "stopdrug_3m_3mFU_MFN_hist", "drugclass")
      ) %>%
      select(-Freq_3m, -total_3m, -Freq_6m, -total_6m, -Freq_12m, -total_12m) %>%
      mutate(variable = "stopdrug_MFN") %>%
      rename("group" = "stopdrug_3m_3mFU_MFN_hist"),
    cprd_dataset %>%
      select(numdrugs, stopdrug_3m_6mFU, drugclass) %>%
      mutate(drugclass = factor(drugclass, levels = c("GLP1", "DPP4", "SGLT2", "TZD", "SU"))) %>%
      rename("analysis_var" = "stopdrug_3m_6mFU") %>%
      table() %>%
      as.data.frame() %>%
      group_by(numdrugs, drugclass) %>%
      mutate(total_3m = sum(Freq), perc_3m = Freq/total_3m) %>%
      rename("Freq_3m" = "Freq") %>%
      filter(analysis_var == 1) %>%
      left_join(
        cprd_dataset %>%
          drop_na(stopdrug_6m_6mFU) %>%
          select(numdrugs, stopdrug_6m_6mFU, drugclass) %>%
          mutate(drugclass = factor(drugclass, levels = c("GLP1", "DPP4", "SGLT2", "TZD", "SU"))) %>%
          rename("analysis_var" = "stopdrug_6m_6mFU") %>%
          table() %>%
          as.data.frame() %>%
          group_by(numdrugs, drugclass) %>%
          mutate(total_6m = sum(Freq), perc_6m = Freq/total_6m) %>%
          rename("Freq_6m" = "Freq") %>%
          filter(analysis_var == 1),
        by = c("analysis_var", "numdrugs", "drugclass")
      ) %>%
      left_join(
        cprd_dataset %>%
          drop_na(stopdrug_12m_6mFU) %>%
          select(numdrugs, stopdrug_12m_6mFU, drugclass) %>%
          mutate(drugclass = factor(drugclass, levels = c("GLP1", "DPP4", "SGLT2", "TZD", "SU"))) %>%
          rename("analysis_var" = "stopdrug_12m_6mFU") %>%
          table() %>%
          as.data.frame() %>%
          group_by(numdrugs, drugclass) %>%
          mutate(total_12m = sum(Freq), perc_12m = Freq/total_12m) %>%
          rename("Freq_12m" = "Freq") %>%
          filter(analysis_var == 1),
        by = c("analysis_var", "numdrugs", "drugclass")
      ) %>%
      select(-Freq_3m, -total_3m, -Freq_6m, -total_6m, -Freq_12m, -total_12m) %>%
      mutate(variable = "numdrugs") %>%
      rename("group" = "numdrugs"),
    cprd_dataset %>%
      select(predrug_statins, stopdrug_3m_6mFU, drugclass) %>%
      mutate(drugclass = factor(drugclass, levels = c("GLP1", "DPP4", "SGLT2", "TZD", "SU"))) %>%
      rename("analysis_var" = "stopdrug_3m_6mFU") %>%
      table() %>%
      as.data.frame() %>%
      filter(Freq != 0) %>%
      group_by(predrug_statins, drugclass) %>%
      mutate(total_3m = sum(Freq), perc_3m = Freq/total_3m) %>%
      rename("Freq_3m" = "Freq") %>%
      filter(analysis_var == 1) %>%
      left_join(
        cprd_dataset %>%
          drop_na(stopdrug_6m_6mFU) %>%
          select(predrug_statins, stopdrug_6m_6mFU, drugclass) %>%
          mutate(drugclass = factor(drugclass, levels = c("GLP1", "DPP4", "SGLT2", "TZD", "SU"))) %>%
          rename("analysis_var" = "stopdrug_6m_6mFU") %>%
          table() %>%
          as.data.frame() %>%
          filter(Freq != 0) %>%
          group_by(predrug_statins, drugclass) %>%
          mutate(total_6m = sum(Freq), perc_6m = Freq/total_6m) %>%
          rename("Freq_6m" = "Freq") %>%
          filter(analysis_var == 1),
        by = c("analysis_var", "predrug_statins", "drugclass")
      ) %>%
      left_join(
        cprd_dataset %>%
          drop_na(stopdrug_12m_6mFU) %>%
          select(predrug_statins, stopdrug_12m_6mFU, drugclass) %>%
          mutate(drugclass = factor(drugclass, levels = c("GLP1", "DPP4", "SGLT2", "TZD", "SU"))) %>%
          rename("analysis_var" = "stopdrug_12m_6mFU") %>%
          table() %>%
          as.data.frame() %>%
          filter(Freq != 0) %>%
          group_by(predrug_statins, drugclass) %>%
          mutate(total_12m = sum(Freq), perc_12m = Freq/total_12m) %>%
          rename("Freq_12m" = "Freq") %>%
          filter(analysis_var == 1),
        by = c("analysis_var", "predrug_statins", "drugclass")
      ) %>%
      select(-Freq_3m, -total_3m, -Freq_6m, -total_6m, -Freq_12m, -total_12m) %>%
      mutate(variable = "statins") %>%
      rename("group" = "predrug_statins"),
    cprd_dataset %>%
      select(predrug_bloodmed, stopdrug_3m_6mFU, drugclass) %>%
      mutate(drugclass = factor(drugclass, levels = c("GLP1", "DPP4", "SGLT2", "TZD", "SU"))) %>%
      rename("analysis_var" = "stopdrug_3m_6mFU") %>%
      table() %>%
      as.data.frame() %>%
      filter(Freq != 0) %>%
      group_by(predrug_bloodmed, drugclass) %>%
      mutate(total_3m = sum(Freq), perc_3m = Freq/total_3m) %>%
      rename("Freq_3m" = "Freq") %>%
      filter(analysis_var == 1) %>%
      left_join(
        cprd_dataset %>%
          drop_na(stopdrug_6m_6mFU) %>%
          select(predrug_bloodmed, stopdrug_6m_6mFU, drugclass) %>%
          mutate(drugclass = factor(drugclass, levels = c("GLP1", "DPP4", "SGLT2", "TZD", "SU"))) %>%
          rename("analysis_var" = "stopdrug_6m_6mFU") %>%
          table() %>%
          as.data.frame() %>%
          filter(Freq != 0) %>%
          group_by(predrug_bloodmed, drugclass) %>%
          mutate(total_6m = sum(Freq), perc_6m = Freq/total_6m) %>%
          rename("Freq_6m" = "Freq") %>%
          filter(analysis_var == 1),
        by = c("analysis_var", "predrug_bloodmed", "drugclass")
      ) %>%
      left_join(
        cprd_dataset %>%
          drop_na(stopdrug_12m_6mFU) %>%
          select(predrug_bloodmed, stopdrug_12m_6mFU, drugclass) %>%
          mutate(drugclass = factor(drugclass, levels = c("GLP1", "DPP4", "SGLT2", "TZD", "SU"))) %>%
          rename("analysis_var" = "stopdrug_12m_6mFU") %>%
          table() %>%
          as.data.frame() %>%
          filter(Freq != 0) %>%
          group_by(predrug_bloodmed, drugclass) %>%
          mutate(total_12m = sum(Freq), perc_12m = Freq/total_12m) %>%
          rename("Freq_12m" = "Freq") %>%
          filter(analysis_var == 1),
        by = c("analysis_var", "predrug_bloodmed", "drugclass")
      ) %>%
      select(-Freq_3m, -total_3m, -Freq_6m, -total_6m, -Freq_12m, -total_12m) %>%
      mutate(variable = "bloodmed") %>%
      rename("group" = "predrug_bloodmed"),
    cprd_dataset %>%
      select(dstartdate_age, stopdrug_3m_6mFU, drugclass) %>%
      mutate(drugclass = factor(drugclass, levels = c("GLP1", "DPP4", "SGLT2", "TZD", "SU"))) %>%
      mutate(group = cut(dstartdate_age, breaks = c(18, 50, 60, 70, 150))) %>%
      rename("analysis_var" = "stopdrug_3m_6mFU") %>%
      select(-dstartdate_age) %>%
      table() %>%
      as.data.frame() %>%
      group_by(group, drugclass) %>%
      mutate(total_3m = sum(Freq), perc_3m = Freq/total_3m) %>%
      rename("Freq_3m" = "Freq") %>%
      filter(analysis_var == 1) %>%
      left_join(
        cprd_dataset %>%
          drop_na(stopdrug_6m_6mFU) %>%
          select(dstartdate_age, stopdrug_6m_6mFU, drugclass) %>%
          mutate(drugclass = factor(drugclass, levels = c("GLP1", "DPP4", "SGLT2", "TZD", "SU"))) %>%
          mutate(group = cut(dstartdate_age, breaks = c(18, 50, 60, 70, 150))) %>%
          rename("analysis_var" = "stopdrug_6m_6mFU") %>%
          select(-dstartdate_age) %>%
          table() %>%
          as.data.frame() %>%
          group_by(group, drugclass) %>%
          mutate(total_6m = sum(Freq), perc_6m = Freq/total_6m) %>%
          rename("Freq_6m" = "Freq") %>%
          filter(analysis_var == 1),
        by = c("drugclass", "group", "analysis_var")
      ) %>%
      left_join(
        cprd_dataset %>%
          drop_na(stopdrug_12m_6mFU) %>%
          select(dstartdate_age, stopdrug_12m_6mFU, drugclass) %>%
          mutate(drugclass = factor(drugclass, levels = c("GLP1", "DPP4", "SGLT2", "TZD", "SU"))) %>%
          mutate(group = cut(dstartdate_age, breaks = c(18, 50, 60, 70, 150))) %>%
          rename("analysis_var" = "stopdrug_12m_6mFU") %>%
          select(-dstartdate_age) %>%
          table() %>%
          as.data.frame() %>%
          group_by(group, drugclass) %>%
          mutate(total_12m = sum(Freq), perc_12m = Freq/total_12m) %>%
          rename("Freq_12m" = "Freq") %>%
          filter(analysis_var == 1),
        by = c("drugclass", "group", "analysis_var")
      ) %>%
      select(-Freq_3m, -total_3m, -Freq_6m, -total_6m, -Freq_12m, -total_12m) %>%
      mutate(variable = "dstartdate_age"),
    cprd_dataset %>%
      select(preegfr, stopdrug_3m_6mFU, drugclass) %>%
      mutate(drugclass = factor(drugclass, levels = c("GLP1", "DPP4", "SGLT2", "TZD", "SU"))) %>%
      mutate(group = cut(preegfr, breaks = c(0, 60, 90, 500))) %>%
      rename("analysis_var" = "stopdrug_3m_6mFU") %>%
      select(-preegfr) %>%
      table() %>%
      as.data.frame() %>%
      group_by(group, drugclass) %>%
      mutate(total_3m = sum(Freq), perc_3m = Freq/total_3m) %>%
      rename("Freq_3m" = "Freq") %>%
      filter(analysis_var == 1) %>%
      left_join(
        cprd_dataset %>%
          drop_na(stopdrug_6m_6mFU) %>%
          select(preegfr, stopdrug_6m_6mFU, drugclass) %>%
          mutate(drugclass = factor(drugclass, levels = c("GLP1", "DPP4", "SGLT2", "TZD", "SU"))) %>%
          mutate(group = cut(preegfr, breaks = c(0, 60, 90, 500))) %>%
          rename("analysis_var" = "stopdrug_6m_6mFU") %>%
          select(-preegfr) %>%
          table() %>%
          as.data.frame() %>%
          group_by(group, drugclass) %>%
          mutate(total_6m = sum(Freq), perc_6m = Freq/total_6m) %>%
          rename("Freq_6m" = "Freq") %>%
          filter(analysis_var == 1),
        by = c("drugclass", "group", "analysis_var")
      ) %>%
      left_join(
        cprd_dataset %>%
          drop_na(stopdrug_12m_6mFU) %>%
          select(preegfr, stopdrug_12m_6mFU, drugclass) %>%
          mutate(drugclass = factor(drugclass, levels = c("GLP1", "DPP4", "SGLT2", "TZD", "SU"))) %>%
          mutate(group = cut(preegfr, breaks = c(0, 60, 90, 500))) %>%
          rename("analysis_var" = "stopdrug_12m_6mFU") %>%
          select(-preegfr) %>%
          table() %>%
          as.data.frame() %>%
          group_by(group, drugclass) %>%
          mutate(total_12m = sum(Freq), perc_12m = Freq/total_12m) %>%
          rename("Freq_12m" = "Freq") %>%
          filter(analysis_var == 1),
        by = c("drugclass", "group", "analysis_var")
      ) %>%
      select(-Freq_3m, -total_3m, -Freq_6m, -total_6m, -Freq_12m, -total_12m) %>%
      mutate(variable = "preegfr"),
    cprd_dataset %>%
      select(gender, stopdrug_3m_6mFU, drugclass) %>%
      mutate(drugclass = factor(drugclass, levels = c("GLP1", "DPP4", "SGLT2", "TZD", "SU"))) %>%
      rename("analysis_var" = "stopdrug_3m_6mFU") %>%
      table() %>%
      as.data.frame() %>%
      group_by(gender, drugclass) %>%
      mutate(total_3m = sum(Freq), perc_3m = Freq/total_3m) %>%
      rename("Freq_3m" = "Freq") %>%
      filter(analysis_var == 1) %>%
      left_join(
        cprd_dataset %>%
          drop_na(stopdrug_6m_6mFU) %>%
          select(gender, stopdrug_6m_6mFU, drugclass) %>%
          mutate(drugclass = factor(drugclass, levels = c("GLP1", "DPP4", "SGLT2", "TZD", "SU"))) %>%
          rename("analysis_var" = "stopdrug_6m_6mFU") %>%
          table() %>%
          as.data.frame() %>%
          group_by(gender, drugclass) %>%
          mutate(total_6m = sum(Freq), perc_6m = Freq/total_6m) %>%
          rename("Freq_6m" = "Freq") %>%
          filter(analysis_var == 1),
        by = c("analysis_var", "gender", "drugclass")
      ) %>%
      left_join(
        cprd_dataset %>%
          drop_na(stopdrug_12m_6mFU) %>%
          select(gender, stopdrug_12m_6mFU, drugclass) %>%
          mutate(drugclass = factor(drugclass, levels = c("GLP1", "DPP4", "SGLT2", "TZD", "SU"))) %>%
          rename("analysis_var" = "stopdrug_12m_6mFU") %>%
          table() %>%
          as.data.frame() %>%
          group_by(gender, drugclass) %>%
          mutate(total_12m = sum(Freq), perc_12m = Freq/total_12m) %>%
          rename("Freq_12m" = "Freq") %>%
          filter(analysis_var == 1),
        by = c("analysis_var", "gender", "drugclass")
      ) %>%
      select(-Freq_3m, -total_3m, -Freq_6m, -total_6m, -Freq_12m, -total_12m) %>%
      mutate(variable = "gender") %>%
      rename("group" = "gender"),
    cprd_dataset %>%
      select(imd2015_10, stopdrug_3m_6mFU, drugclass) %>%
      mutate(drugclass = factor(drugclass, levels = c("GLP1", "DPP4", "SGLT2", "TZD", "SU"))) %>%
      rename("analysis_var" = "stopdrug_3m_6mFU") %>%
      table() %>%
      as.data.frame() %>%
      group_by(imd2015_10, drugclass) %>%
      mutate(total_3m = sum(Freq), perc_3m = Freq/total_3m) %>%
      rename("Freq_3m" = "Freq") %>%
      filter(analysis_var == 1) %>%
      left_join(
        cprd_dataset %>%
          drop_na(stopdrug_6m_6mFU) %>%
          select(imd2015_10, stopdrug_6m_6mFU, drugclass) %>%
          mutate(drugclass = factor(drugclass, levels = c("GLP1", "DPP4", "SGLT2", "TZD", "SU"))) %>%
          rename("analysis_var" = "stopdrug_6m_6mFU") %>%
          table() %>%
          as.data.frame() %>%
          group_by(imd2015_10, drugclass) %>%
          mutate(total_6m = sum(Freq), perc_6m = Freq/total_6m) %>%
          rename("Freq_6m" = "Freq") %>%
          filter(analysis_var == 1),
        by = c("analysis_var", "imd2015_10", "drugclass")
      ) %>%
      left_join(
        cprd_dataset %>%
          drop_na(stopdrug_12m_6mFU) %>%
          select(imd2015_10, stopdrug_12m_6mFU, drugclass) %>%
          mutate(drugclass = factor(drugclass, levels = c("GLP1", "DPP4", "SGLT2", "TZD", "SU"))) %>%
          rename("analysis_var" = "stopdrug_12m_6mFU") %>%
          table() %>%
          as.data.frame() %>%
          group_by(imd2015_10, drugclass) %>%
          mutate(total_12m = sum(Freq), perc_12m = Freq/total_12m) %>%
          rename("Freq_12m" = "Freq") %>%
          filter(analysis_var == 1),
        by = c("analysis_var", "imd2015_10", "drugclass")
      ) %>%
      select(-Freq_3m, -total_3m, -Freq_6m, -total_6m, -Freq_12m, -total_12m) %>%
      mutate(variable = "imd2015") %>%
      rename("group" = "imd2015_10"),
    cprd_dataset %>%
      select(smoking_cat, stopdrug_3m_6mFU, drugclass) %>%
      mutate(drugclass = factor(drugclass, levels = c("GLP1", "DPP4", "SGLT2", "TZD", "SU"))) %>%
      rename("analysis_var" = "stopdrug_3m_6mFU") %>%
      table() %>%
      as.data.frame() %>%
      group_by(smoking_cat, drugclass) %>%
      mutate(total_3m = sum(Freq), perc_3m = Freq/total_3m) %>%
      rename("Freq_3m" = "Freq") %>%
      filter(analysis_var == 1) %>%
      left_join(
        cprd_dataset %>%
          drop_na(stopdrug_6m_6mFU) %>%
          select(smoking_cat, stopdrug_6m_6mFU, drugclass) %>%
          mutate(drugclass = factor(drugclass, levels = c("GLP1", "DPP4", "SGLT2", "TZD", "SU"))) %>%
          rename("analysis_var" = "stopdrug_6m_6mFU") %>%
          table() %>%
          as.data.frame() %>%
          group_by(smoking_cat, drugclass) %>%
          mutate(total_6m = sum(Freq), perc_6m = Freq/total_6m) %>%
          rename("Freq_6m" = "Freq") %>%
          filter(analysis_var == 1),
        by = c("analysis_var", "smoking_cat", "drugclass")
      ) %>%
      left_join(
        cprd_dataset %>%
          drop_na(stopdrug_12m_6mFU) %>%
          select(smoking_cat, stopdrug_12m_6mFU, drugclass) %>%
          mutate(drugclass = factor(drugclass, levels = c("GLP1", "DPP4", "SGLT2", "TZD", "SU"))) %>%
          rename("analysis_var" = "stopdrug_12m_6mFU") %>%
          table() %>%
          as.data.frame() %>%
          group_by(smoking_cat, drugclass) %>%
          mutate(total_12m = sum(Freq), perc_12m = Freq/total_12m) %>%
          rename("Freq_12m" = "Freq") %>%
          filter(analysis_var == 1),
        by = c("analysis_var", "smoking_cat", "drugclass")
      ) %>%
      select(-Freq_3m, -total_3m, -Freq_6m, -total_6m, -Freq_12m, -total_12m) %>%
      mutate(variable = "smoking") %>%
      rename("group" = "smoking_cat"),
    cprd_dataset %>%
      select(prebmi, stopdrug_3m_6mFU, drugclass) %>%
      mutate(drugclass = factor(drugclass, levels = c("GLP1", "DPP4", "SGLT2", "TZD", "SU"))) %>%
      mutate(group = cut(prebmi, breaks = c(1, 25, 30, 35, 500))) %>%
      rename("analysis_var" = "stopdrug_3m_6mFU") %>%
      select(-prebmi) %>%
      table() %>%
      as.data.frame() %>%
      group_by(group, drugclass) %>%
      mutate(total_3m = sum(Freq), perc_3m = Freq/total_3m) %>%
      rename("Freq_3m" = "Freq") %>%
      filter(analysis_var == 1) %>%
      left_join(
        cprd_dataset %>%
          drop_na(stopdrug_6m_6mFU) %>%
          select(prebmi, stopdrug_6m_6mFU, drugclass) %>%
          mutate(drugclass = factor(drugclass, levels = c("GLP1", "DPP4", "SGLT2", "TZD", "SU"))) %>%
          mutate(group = cut(prebmi, breaks = c(1, 25, 30, 35, 500))) %>%
          rename("analysis_var" = "stopdrug_6m_6mFU") %>%
          select(-prebmi) %>%
          table() %>%
          as.data.frame() %>%
          group_by(group, drugclass) %>%
          mutate(total_6m = sum(Freq), perc_6m = Freq/total_6m) %>%
          rename("Freq_6m" = "Freq") %>%
          filter(analysis_var == 1),
        by = c("drugclass", "group", "analysis_var")
      ) %>%
      left_join(
        cprd_dataset %>%
          drop_na(stopdrug_12m_6mFU) %>%
          select(prebmi, stopdrug_12m_6mFU, drugclass) %>%
          mutate(drugclass = factor(drugclass, levels = c("GLP1", "DPP4", "SGLT2", "TZD", "SU"))) %>%
          mutate(group = cut(prebmi, breaks = c(1, 25, 30, 35, 500))) %>%
          rename("analysis_var" = "stopdrug_12m_6mFU") %>%
          select(-prebmi) %>%
          table() %>%
          as.data.frame() %>%
          group_by(group, drugclass) %>%
          mutate(total_12m = sum(Freq), perc_12m = Freq/total_12m) %>%
          rename("Freq_12m" = "Freq") %>%
          filter(analysis_var == 1),
        by = c("drugclass", "group", "analysis_var")
      ) %>%
      select(-Freq_3m, -total_3m, -Freq_6m, -total_6m, -Freq_12m, -total_12m) %>%
      mutate(variable = "bmi"),
    cprd_dataset %>%
      select(dstartdate_dm_dur, stopdrug_3m_6mFU, drugclass) %>%
      mutate(drugclass = factor(drugclass, levels = c("GLP1", "DPP4", "SGLT2", "TZD", "SU"))) %>%
      mutate(group = cut(dstartdate_dm_dur, breaks = c(-1, 1, 3, 6, 10, 150))) %>%
      rename("analysis_var" = "stopdrug_3m_6mFU") %>%
      select(-dstartdate_dm_dur) %>%
      table() %>%
      as.data.frame() %>%
      group_by(group, drugclass) %>%
      mutate(total_3m = sum(Freq), perc_3m = Freq/total_3m) %>%
      rename("Freq_3m" = "Freq") %>%
      filter(analysis_var == 1) %>%
      left_join(
        cprd_dataset %>%
          drop_na(stopdrug_6m_6mFU) %>%
          select(dstartdate_dm_dur, stopdrug_6m_6mFU, drugclass) %>%
          mutate(drugclass = factor(drugclass, levels = c("GLP1", "DPP4", "SGLT2", "TZD", "SU"))) %>%
          mutate(group = cut(dstartdate_dm_dur, breaks = c(-1, 1, 3, 6, 10, 150))) %>%
          rename("analysis_var" = "stopdrug_6m_6mFU") %>%
          select(-dstartdate_dm_dur) %>%
          table() %>%
          as.data.frame() %>%
          group_by(group, drugclass) %>%
          mutate(total_6m = sum(Freq), perc_6m = Freq/total_6m) %>%
          rename("Freq_6m" = "Freq") %>%
          filter(analysis_var == 1),
        by = c("drugclass", "group", "analysis_var")
      ) %>%
      left_join(
        cprd_dataset %>%
          drop_na(stopdrug_12m_6mFU) %>%
          select(dstartdate_dm_dur, stopdrug_12m_6mFU, drugclass) %>%
          mutate(drugclass = factor(drugclass, levels = c("GLP1", "DPP4", "SGLT2", "TZD", "SU"))) %>%
          mutate(group = cut(dstartdate_dm_dur, breaks = c(-1, 1, 3, 6, 10, 150))) %>%
          rename("analysis_var" = "stopdrug_12m_6mFU") %>%
          select(-dstartdate_dm_dur) %>%
          table() %>%
          as.data.frame() %>%
          group_by(group, drugclass) %>%
          mutate(total_12m = sum(Freq), perc_12m = Freq/total_12m) %>%
          rename("Freq_12m" = "Freq") %>%
          filter(analysis_var == 1),
        by = c("drugclass", "group", "analysis_var")
      ) %>%
      select(-Freq_3m, -total_3m, -Freq_6m, -total_6m, -Freq_12m, -total_12m) %>%
      mutate(variable = "dstartdate_dm_dur"),
    cprd_dataset %>%
      select(prehba1c, stopdrug_3m_6mFU, drugclass) %>%
      mutate(drugclass = factor(drugclass, levels = c("GLP1", "DPP4", "SGLT2", "TZD", "SU"))) %>%
      mutate(group = cut(prehba1c, breaks = c(52, 64, 75, 86, 250))) %>%
      rename("analysis_var" = "stopdrug_3m_6mFU") %>%
      select(-prehba1c) %>%
      table() %>%
      as.data.frame() %>%
      group_by(group, drugclass) %>%
      mutate(total_3m = sum(Freq), perc_3m = Freq/total_3m) %>%
      rename("Freq_3m" = "Freq") %>%
      filter(analysis_var == 1) %>%
      left_join(
        cprd_dataset %>%
          drop_na(stopdrug_6m_6mFU) %>%
          select(prehba1c, stopdrug_6m_6mFU, drugclass) %>%
          mutate(drugclass = factor(drugclass, levels = c("GLP1", "DPP4", "SGLT2", "TZD", "SU"))) %>%
          mutate(group = cut(prehba1c, breaks = c(52, 64, 75, 86, 250))) %>%
          rename("analysis_var" = "stopdrug_6m_6mFU") %>%
          select(-prehba1c) %>%
          table() %>%
          as.data.frame() %>%
          group_by(group, drugclass) %>%
          mutate(total_6m = sum(Freq), perc_6m = Freq/total_6m) %>%
          rename("Freq_6m" = "Freq") %>%
          filter(analysis_var == 1),
        by = c("drugclass", "group", "analysis_var")
      ) %>%
      left_join(
        cprd_dataset %>%
          drop_na(stopdrug_12m_6mFU) %>%
          select(prehba1c, stopdrug_12m_6mFU, drugclass) %>%
          mutate(drugclass = factor(drugclass, levels = c("GLP1", "DPP4", "SGLT2", "TZD", "SU"))) %>%
          mutate(group = cut(prehba1c, breaks = c(52, 64, 75, 86, 250))) %>%
          rename("analysis_var" = "stopdrug_12m_6mFU") %>%
          select(-prehba1c) %>%
          table() %>%
          as.data.frame() %>%
          group_by(group, drugclass) %>%
          mutate(total_12m = sum(Freq), perc_12m = Freq/total_12m) %>%
          rename("Freq_12m" = "Freq") %>%
          filter(analysis_var == 1),
        by = c("drugclass", "group", "analysis_var")
      ) %>%
      select(-Freq_3m, -total_3m, -Freq_6m, -total_6m, -Freq_12m, -total_12m) %>%
      mutate(variable = "hba1c"),
    cprd_dataset %>%
      select(drugline, stopdrug_3m_6mFU, drugclass) %>%
      mutate(drugline = factor(drugline, levels = c("2", "3", "4", "5+"), labels = c("2", "3", "4+", "4+"))) %>%
      mutate(drugclass = factor(drugclass, levels = c("GLP1", "DPP4", "SGLT2", "TZD", "SU"))) %>%
      rename("analysis_var" = "stopdrug_3m_6mFU") %>%
      table() %>%
      as.data.frame() %>%
      group_by(drugline, drugclass) %>%
      mutate(total_3m = sum(Freq), perc_3m = Freq/total_3m) %>%
      rename("Freq_3m" = "Freq") %>%
      filter(analysis_var == 1) %>%
      left_join(
        cprd_dataset %>%
          drop_na(stopdrug_6m_6mFU) %>%
          select(drugline, stopdrug_6m_6mFU, drugclass) %>%
          mutate(drugline = factor(drugline, levels = c("2", "3", "4", "5+"), labels = c("2", "3", "4+", "4+"))) %>%
          mutate(drugclass = factor(drugclass, levels = c("GLP1", "DPP4", "SGLT2", "TZD", "SU"))) %>%
          rename("analysis_var" = "stopdrug_6m_6mFU") %>%
          table() %>%
          as.data.frame() %>%
          group_by(drugline, drugclass) %>%
          mutate(total_6m = sum(Freq), perc_6m = Freq/total_6m) %>%
          rename("Freq_6m" = "Freq") %>%
          filter(analysis_var == 1),
        by = c("analysis_var", "drugline", "drugclass")
      ) %>%
      left_join(
        cprd_dataset %>%
          drop_na(stopdrug_12m_6mFU) %>%
          select(drugline, stopdrug_12m_6mFU, drugclass) %>%
          mutate(drugline = factor(drugline, levels = c("2", "3", "4", "5+"), labels = c("2", "3", "4+", "4+"))) %>%
          mutate(drugclass = factor(drugclass, levels = c("GLP1", "DPP4", "SGLT2", "TZD", "SU"))) %>%
          rename("analysis_var" = "stopdrug_12m_6mFU") %>%
          table() %>%
          as.data.frame() %>%
          group_by(drugline, drugclass) %>%
          mutate(total_12m = sum(Freq), perc_12m = Freq/total_12m) %>%
          rename("Freq_12m" = "Freq") %>%
          filter(analysis_var == 1),
        by = c("analysis_var", "drugline", "drugclass")
      ) %>%
      select(-Freq_3m, -total_3m, -Freq_6m, -total_6m, -Freq_12m, -total_12m) %>%
      mutate(variable = "drugline") %>%
      rename("group" = "drugline"),
    cprd_dataset %>%
      select(predrug_frailty_proxy, stopdrug_3m_6mFU, drugclass) %>%
      mutate(drugclass = factor(drugclass, levels = c("GLP1", "DPP4", "SGLT2", "TZD", "SU"))) %>%
      rename("analysis_var" = "stopdrug_3m_6mFU") %>%
      table() %>%
      as.data.frame() %>%
      group_by(predrug_frailty_proxy, drugclass) %>%
      mutate(total_3m = sum(Freq), perc_3m = Freq/total_3m) %>%
      rename("Freq_3m" = "Freq") %>%
      filter(analysis_var == 1) %>%
      left_join(
        cprd_dataset %>%
          drop_na(stopdrug_6m_6mFU) %>%
          select(predrug_frailty_proxy, stopdrug_6m_6mFU, drugclass) %>%
          mutate(drugclass = factor(drugclass, levels = c("GLP1", "DPP4", "SGLT2", "TZD", "SU"))) %>%
          rename("analysis_var" = "stopdrug_6m_6mFU") %>%
          table() %>%
          as.data.frame() %>%
          group_by(predrug_frailty_proxy, drugclass) %>%
          mutate(total_6m = sum(Freq), perc_6m = Freq/total_6m) %>%
          rename("Freq_6m" = "Freq") %>%
          filter(analysis_var == 1),
        by = c("analysis_var", "predrug_frailty_proxy", "drugclass")
      ) %>%
      left_join(
        cprd_dataset %>%
          drop_na(stopdrug_12m_6mFU) %>%
          select(predrug_frailty_proxy, stopdrug_12m_6mFU, drugclass) %>%
          mutate(drugclass = factor(drugclass, levels = c("GLP1", "DPP4", "SGLT2", "TZD", "SU"))) %>%
          rename("analysis_var" = "stopdrug_12m_6mFU") %>%
          table() %>%
          as.data.frame() %>%
          group_by(predrug_frailty_proxy, drugclass) %>%
          mutate(total_12m = sum(Freq), perc_12m = Freq/total_12m) %>%
          rename("Freq_12m" = "Freq") %>%
          filter(analysis_var == 1),
        by = c("analysis_var", "predrug_frailty_proxy", "drugclass")
      ) %>%
      select(-Freq_3m, -total_3m, -Freq_6m, -total_6m, -Freq_12m, -total_12m) %>%
      mutate(variable = "frailty") %>%
      rename("group" = "predrug_frailty_proxy"),
    cprd_dataset %>%
      select(ethnicity_5cat, stopdrug_3m_6mFU, drugclass) %>%
      mutate(ethnicity_5cat = factor(ethnicity_5cat, levels = c("4", "3", "2", "1", "0"), labels = c("Mixed & Other", "Mixed & Other", "Black", "South Asian", "White"))) %>%
      mutate(drugclass = factor(drugclass, levels = c("GLP1", "DPP4", "SGLT2", "TZD", "SU"))) %>%
      rename("analysis_var" = "stopdrug_3m_6mFU") %>%
      table() %>%
      as.data.frame() %>%
      group_by(ethnicity_5cat, drugclass) %>%
      mutate(total_3m = sum(Freq), perc_3m = Freq/total_3m) %>%
      rename("Freq_3m" = "Freq") %>%
      filter(analysis_var == 1) %>%
      left_join(
        cprd_dataset %>%
          drop_na(stopdrug_6m_6mFU) %>%
          select(ethnicity_5cat, stopdrug_6m_6mFU, drugclass) %>%
          mutate(ethnicity_5cat = factor(ethnicity_5cat, levels = c("4", "3", "2", "1", "0"), labels = c("Mixed & Other", "Mixed & Other", "Black", "South Asian", "White"))) %>%
          mutate(drugclass = factor(drugclass, levels = c("GLP1", "DPP4", "SGLT2", "TZD", "SU"))) %>%
          rename("analysis_var" = "stopdrug_6m_6mFU") %>%
          table() %>%
          as.data.frame() %>%
          group_by(ethnicity_5cat, drugclass) %>%
          mutate(total_6m = sum(Freq), perc_6m = Freq/total_6m) %>%
          rename("Freq_6m" = "Freq") %>%
          filter(analysis_var == 1),
        by = c("analysis_var", "ethnicity_5cat", "drugclass")
      ) %>%
      left_join(
        cprd_dataset %>%
          drop_na(stopdrug_12m_6mFU) %>%
          select(ethnicity_5cat, stopdrug_12m_6mFU, drugclass) %>%
          mutate(ethnicity_5cat = factor(ethnicity_5cat, levels = c("4", "3", "2", "1", "0"), labels = c("Mixed & Other", "Mixed & Other", "Black", "South Asian", "White"))) %>%
          mutate(drugclass = factor(drugclass, levels = c("GLP1", "DPP4", "SGLT2", "TZD", "SU"))) %>%
          rename("analysis_var" = "stopdrug_12m_6mFU") %>%
          table() %>%
          as.data.frame() %>%
          group_by(ethnicity_5cat, drugclass) %>%
          mutate(total_12m = sum(Freq), perc_12m = Freq/total_12m) %>%
          rename("Freq_12m" = "Freq") %>%
          filter(analysis_var == 1),
        by = c("analysis_var", "ethnicity_5cat", "drugclass")
      ) %>%
      select(-Freq_3m, -total_3m, -Freq_6m, -total_6m, -Freq_12m, -total_12m) %>%
      mutate(variable = "ethnicity") %>%
      rename("group" = "ethnicity_5cat")
  )





###############################################################################
###############################################################################

# plot for the therapy discontinuation

pdf("results/figures/03.discontinuation_per_drug_combined.pdf", width = 6, height = 5)

cprd_dataset %>%
  group_by(drugclass) %>%
  mutate(total_num_drug = n()) %>%
  group_by(drugclass, stopdrug_3m_6mFU) %>%
  mutate(total_num_drug_disc_3m = n()) %>%
  ungroup() %>%
  mutate(perc_3m = total_num_drug_disc_3m/total_num_drug) %>%
  select(drugclass, perc_3m, stopdrug_3m_6mFU) %>%
  unique() %>%
  filter(stopdrug_3m_6mFU == 1) %>%
  select(-stopdrug_3m_6mFU) %>%
  left_join(
    cprd_dataset %>%
      drop_na(stopdrug_6m_6mFU) %>%
      group_by(drugclass) %>%
      mutate(total_num_drug = n()) %>%
      group_by(drugclass, stopdrug_6m_6mFU) %>%
      mutate(total_num_drug_disc_6m = n()) %>%
      ungroup() %>%
      mutate(perc_6m = total_num_drug_disc_6m/total_num_drug) %>%
      select(drugclass, perc_6m, stopdrug_6m_6mFU) %>%
      unique() %>%
      filter(stopdrug_6m_6mFU == 1) %>%
      select(-stopdrug_6m_6mFU),
    by = c("drugclass")
  ) %>%
  left_join(
    cprd_dataset %>%
      drop_na(stopdrug_12m_6mFU) %>%
      group_by(drugclass) %>%
      mutate(total_num_drug = n()) %>%
      group_by(drugclass, stopdrug_12m_6mFU) %>%
      mutate(total_num_drug_disc_12m = n()) %>%
      ungroup() %>%
      mutate(perc_12m = total_num_drug_disc_12m/total_num_drug) %>%
      select(drugclass, perc_12m, stopdrug_12m_6mFU) %>%
      unique() %>%
      filter(stopdrug_12m_6mFU == 1) %>%
      select(-stopdrug_12m_6mFU),
    by = c("drugclass")
  ) %>%
  mutate(
    perc_12m = perc_12m-perc_6m,
    perc_6m = perc_6m-perc_3m
  ) %>%
  gather("key", "value", -drugclass) %>%
  mutate(
    drugclass = factor(drugclass, levels = c("GLP1", "DPP4", "SGLT2", "TZD", "SU")),
    key = factor(key, levels = c("perc_12m", "perc_6m", "perc_3m"))
  ) %>%
  ggplot(aes(x = value, y = drugclass, fill = drugclass, alpha = key)) +
  geom_bar(stat = "identity", colour = "black") +
  scale_fill_manual(values = c("Pooled" = "black", "SGLT2" = "#E69F00", "GLP1" = "#56B4E9", "SU" = "#CC79A7", "DPP4" = "#0072B2", "TZD" = "#D55E00", "MFN" = "grey"), 
                    breaks = rev(c("Pooled", "GLP1", "DPP4", "SGLT2", "TZD", "SU")), labels = rev(c("Pooled", "GLP-1RA", "DPP4i", "SGLT2i", "TZD", "SU")),
                    name = "Therapy", guide = guide_legend(reverse = TRUE)) +
  scale_y_discrete(limits = rev, breaks = c("GLP1", "DPP4", "SGLT2", "TZD", "SU"), labels = c("GLP-1RA", "DPP4i", "SGLT2i", "TZD", "SU")) +
  scale_alpha_discrete(limits = rev, range = c(1, 0.4), breaks = c("perc_12m", "perc_6m", "perc_3m"), labels = c("12-month", "6-month", "3-month"), name = "Discontinuation", guide = guide_legend(reverse = TRUE, row = 1)) +
  scale_x_continuous("Proportion of discontinuation", labels = scales::percent, limits = c(0, 0.4), breaks = seq(0, 0.4, 0.1)) +
  guides(fill = "none") +
  theme_bw() +
  theme(
    legend.position = "bottom",
    legend.box = "vertical",
    plot.title.position = "plot",
    axis.title.y = element_blank(),
    axis.text.x = element_text(size = 14),
    axis.text.y = element_text(size = 14),
    axis.title.x = element_text(size = 15)
  )



# cprd_dataset %>%
#   filter(dstartdate > "2017-01-01") %>%
#   group_by(drugclass) %>%
#   mutate(total_num_drug = n()) %>%
#   group_by(drugclass, stopdrug_3m_6mFU) %>%
#   mutate(total_num_drug_disc_3m = n()) %>%
#   ungroup() %>%
#   mutate(perc_3m = total_num_drug_disc_3m/total_num_drug) %>%
#   select(drugclass, perc_3m, stopdrug_3m_6mFU) %>%
#   unique() %>%
#   filter(stopdrug_3m_6mFU == 1) %>%
#   select(-stopdrug_3m_6mFU) %>%
#   left_join(
#     cprd_dataset %>%
#       filter(dstartdate > "2017-01-01") %>%
#       drop_na(stopdrug_6m_6mFU) %>%
#       group_by(drugclass) %>%
#       mutate(total_num_drug = n()) %>%
#       group_by(drugclass, stopdrug_6m_6mFU) %>%
#       mutate(total_num_drug_disc_6m = n()) %>%
#       ungroup() %>%
#       mutate(perc_6m = total_num_drug_disc_6m/total_num_drug) %>%
#       select(drugclass, perc_6m, stopdrug_6m_6mFU) %>%
#       unique() %>%
#       filter(stopdrug_6m_6mFU == 1) %>%
#       select(-stopdrug_6m_6mFU),
#     by = c("drugclass")
#   ) %>%
#   left_join(
#     cprd_dataset %>%
#       filter(dstartdate > "2017-01-01") %>%
#       drop_na(stopdrug_12m_6mFU) %>%
#       group_by(drugclass) %>%
#       mutate(total_num_drug = n()) %>%
#       group_by(drugclass, stopdrug_12m_6mFU) %>%
#       mutate(total_num_drug_disc_12m = n()) %>%
#       ungroup() %>%
#       mutate(perc_12m = total_num_drug_disc_12m/total_num_drug) %>%
#       select(drugclass, perc_12m, stopdrug_12m_6mFU) %>%
#       unique() %>%
#       filter(stopdrug_12m_6mFU == 1) %>%
#       select(-stopdrug_12m_6mFU),
#     by = c("drugclass")
#   ) %>%
#   mutate(
#     perc_12m = perc_12m-perc_6m,
#     perc_6m = perc_6m-perc_3m
#   ) %>%
#   gather("key", "value", -drugclass) %>%
#   mutate(
#     drugclass = factor(drugclass, levels = c("GLP1", "DPP4", "SGLT2", "TZD", "SU")),
#     key = factor(key, levels = c("perc_12m", "perc_6m", "perc_3m"))
#   ) %>%
#   ggplot(aes(x = value, y = drugclass, fill = drugclass, alpha = key)) +
#   geom_bar(stat = "identity", colour = "black") +
#   scale_fill_manual(values = c("Pooled" = "black", "SGLT2" = "#E69F00", "GLP1" = "#56B4E9", "SU" = "#CC79A7", "DPP4" = "#0072B2", "TZD" = "#D55E00", "MFN" = "grey"), 
#                     breaks = rev(c("Pooled", "GLP1", "DPP4", "SGLT2", "TZD", "SU")), labels = rev(c("Pooled", "GLP-1RA", "DPP4i", "SGLT2i", "TZD", "SU")),
#                     name = "Therapy", guide = guide_legend(reverse = TRUE)) +
#   scale_y_discrete(limits = rev, breaks = c("GLP1", "DPP4", "SGLT2", "TZD", "SU"), labels = c("GLP-1RA", "DPP4i", "SGLT2i", "TZD", "SU")) +
#   scale_alpha_discrete(limits = rev, range = c(1, 0.4), breaks = c("perc_12m", "perc_6m", "perc_3m"), labels = c("12-month", "6-month", "3-month"), name = "Discontinuation", guide = guide_legend(reverse = TRUE, row = 1)) +
#   scale_x_continuous("Proportion of discontinuation", labels = scales::percent, limits = c(0, 0.35), breaks = seq(0, 0.35, 0.05)) +
#   ggtitle("Therapies (2017-onwards)") +
#   guides(fill = "none") +
#   theme_bw() +
#   theme(
#     legend.position = "bottom",
#     legend.box = "vertical",
#     plot.title.position = "plot",
#     axis.title.y = element_blank(),
#     axis.text.x = element_text(size = 14),
#     axis.text.y = element_text(size = 14),
#     axis.title.x = element_text(size = 15)
#   )


dev.off()

###############################################################################
###############################################################################


pdf("results/figures/03.discontinuation_by_variable.pdf", width = 7, height = 15)

patchwork::wrap_plots(
  
  # Overall
  discontinuation_per_drug_combined %>%
    mutate(
      perc_12m = perc_12m-perc_6m,
      perc_6m = perc_6m-perc_3m
    ) %>%
    gather("key", "value", -drugclass, -analysis_var, -group, -variable) %>%
    mutate(
      key = factor(key, levels = c("perc_12m", "perc_6m", "perc_3m"))
    ) %>%
    filter(variable == "overall") %>%
    ggplot(aes(y = analysis_var, x = value, fill = drugclass, alpha = key)) +
    geom_bar(stat = "identity", colour = "black") +
    scale_fill_manual(values = c("other" = "black", "SGLT2" = "#E69F00", "GLP1" = "#56B4E9", "SU" = "#CC79A7", "DPP4" = "#0072B2", "TZD" = "#D55E00", "MFN" = "grey"), breaks = c("other"), labels = c("Other glucose-lowering therapies"), name = " ", guide = guide_legend(reverse = FALSE)) +
    scale_y_discrete(breaks = c(1), labels = c("Overall")) +
    scale_alpha_discrete(limits = rev, range = c(1, 0.4), breaks = c("perc_12m", "perc_6m", "perc_3m"), labels = c("12-month", "6-month", "3-month"), name = "Discontinuation", guide = guide_legend(reverse = TRUE, row = 1)) +
    scale_x_continuous("Proportion of discontinuation", labels = scales::percent, limits = c(0, 0.4)) +
    facet_wrap(~drugclass, nrow = 1, labeller = as_labeller(c("other" = "Other glucose-lowering therapies"))) +
    theme_bw() +
    guides(fill = "none") +
    theme(
      legend.position = "bottom",
      axis.text.x = element_blank(),
      axis.title.x = element_blank()
    ),
  
  # Age
  discontinuation_per_drug_combined %>%
    mutate(
      perc_12m = perc_12m-perc_6m,
      perc_6m = perc_6m-perc_3m
    ) %>%
    gather("key", "value", -drugclass, -analysis_var, -group, -variable) %>%
    mutate(
      key = factor(key, levels = c("perc_12m", "perc_6m", "perc_3m"))
    ) %>%
    filter(variable == "dstartdate_age") %>%
    ggplot(aes(y = group, x = value, fill = drugclass, alpha = key)) +
    geom_bar(stat = "identity", colour = "black") +
    scale_fill_manual(values = c("other" = "black", "SGLT2" = "#E69F00", "GLP1" = "#56B4E9", "SU" = "#CC79A7", "DPP4" = "#0072B2", "TZD" = "#D55E00", "MFN" = "grey"), breaks = c("other"), labels = c("Other glucose-lowering therapies"), name = " ", guide = guide_legend(reverse = FALSE)) +
    scale_y_discrete(limits = rev, labels = c("70+", "60-69", "50-59", "<50")) +
    scale_alpha_discrete(limits = rev, range = c(1, 0.4), breaks = c("perc_12m", "perc_6m", "perc_3m"), labels = c("12-month", "6-month", "3-month"), name = "Discontinuation", guide = guide_legend(reverse = TRUE, row = 1)) +
    scale_x_continuous("Proportion of discontinuation", labels = scales::percent, limits = c(0, 0.4)) +
    ggtitle("Age group, years") +
    facet_wrap(~drugclass, nrow = 1, labeller = as_labeller(c("other" = "Other glucose-lowering therapies"))) +
    theme_bw() +
    guides(fill = "none") +
    theme(
      plot.title.position = "plot",
      legend.position = "bottom",
      axis.text.x = element_blank(),
      axis.title.x = element_blank(),
      strip.text.x = element_blank()
    ),
  
  # Diabetes duration
  discontinuation_per_drug_combined %>%
    mutate(
      perc_12m = perc_12m-perc_6m,
      perc_6m = perc_6m-perc_3m
    ) %>%
    gather("key", "value", -drugclass, -analysis_var, -group, -variable) %>%
    mutate(
      key = factor(key, levels = c("perc_12m", "perc_6m", "perc_3m"))
    ) %>%
    filter(variable == "dstartdate_dm_dur") %>%
    ggplot(aes(y = group, x = value, fill = drugclass, alpha = key)) +
    geom_bar(stat = "identity", colour = "black") +
    scale_fill_manual(values = c("other" = "black", "SGLT2" = "#E69F00", "GLP1" = "#56B4E9", "SU" = "#CC79A7", "DPP4" = "#0072B2", "TZD" = "#D55E00", "MFN" = "grey"), breaks = c("other"), labels = c("Other glucose-lowering therapies"), name = " ", guide = guide_legend(reverse = FALSE)) +
    scale_y_discrete(limits = rev, labels = c("10+", "6-9", "3-5", "1-2", "<1")) +
    scale_alpha_discrete(limits = rev, range = c(1, 0.4), breaks = c("perc_12m", "perc_6m", "perc_3m"), labels = c("12-month", "6-month", "3-month"), name = "Discontinuation", guide = guide_legend(reverse = TRUE, row = 1)) +
    scale_x_continuous("Proportion of discontinuation", labels = scales::percent, limits = c(0, 0.4)) +
    ggtitle("Duration of diabetes, years") +
    facet_wrap(~drugclass, nrow = 1, labeller = as_labeller(c("other" = "Other glucose-lowering therapies"))) +
    theme_bw() +
    guides(fill = "none") +
    theme(
      plot.title.position = "plot",
      legend.position = "bottom",
      axis.text.x = element_blank(),
      axis.title.x = element_blank(),
      strip.text.x = element_blank()
    ),
  
  # Sex
  discontinuation_per_drug_combined %>%
    mutate(
      perc_12m = perc_12m-perc_6m,
      perc_6m = perc_6m-perc_3m
    ) %>%
    gather("key", "value", -drugclass, -analysis_var, -group, -variable) %>%
    mutate(
      key = factor(key, levels = c("perc_12m", "perc_6m", "perc_3m"))
    ) %>%
    filter(variable == "gender") %>%
    ggplot(aes(y = group, x = value, fill = drugclass, alpha = key)) +
    geom_bar(stat = "identity", colour = "black") +
    scale_fill_manual(values = c("other" = "black", "SGLT2" = "#E69F00", "GLP1" = "#56B4E9", "SU" = "#CC79A7", "DPP4" = "#0072B2", "TZD" = "#D55E00", "MFN" = "grey"), breaks = c("other"), labels = c("Other glucose-lowering therapies"), name = " ", guide = guide_legend(reverse = FALSE)) +
    scale_y_discrete(limits = rev, labels = c("Female", "Male")) +
    scale_alpha_discrete(limits = rev, range = c(1, 0.4), breaks = c("perc_12m", "perc_6m", "perc_3m"), labels = c("12-month", "6-month", "3-month"), name = "Discontinuation", guide = guide_legend(reverse = TRUE, row = 1)) +
    scale_x_continuous("Proportion of discontinuation", labels = scales::percent, limits = c(0, 0.4)) +
    ggtitle("Sex") +
    facet_wrap(~drugclass, nrow = 1, labeller = as_labeller(c("other" = "Other glucose-lowering therapies"))) +
    theme_bw() +
    guides(fill = "none") +
    theme(
      plot.title.position = "plot",
      legend.position = "bottom",
      axis.text.x = element_blank(),
      axis.title.x = element_blank(),
      strip.text.x = element_blank()
    ),
  
  # Smoking
  discontinuation_per_drug_combined %>%
    mutate(
      perc_12m = perc_12m-perc_6m,
      perc_6m = perc_6m-perc_3m
    ) %>%
    gather("key", "value", -drugclass, -analysis_var, -group, -variable) %>%
    mutate(
      key = factor(key, levels = c("perc_12m", "perc_6m", "perc_3m"))
    ) %>%
    filter(variable == "smoking") %>%
    ggplot(aes(y = group, x = value, fill = drugclass, alpha = key)) +
    geom_bar(stat = "identity", colour = "black") +
    scale_fill_manual(values = c("other" = "black", "SGLT2" = "#E69F00", "GLP1" = "#56B4E9", "SU" = "#CC79A7", "DPP4" = "#0072B2", "TZD" = "#D55E00", "MFN" = "grey"), breaks = c("other"), labels = c("Other glucose-lowering therapies"), name = " ", guide = guide_legend(reverse = FALSE)) +
    scale_alpha_discrete(limits = rev, range = c(1, 0.4), breaks = c("perc_12m", "perc_6m", "perc_3m"), labels = c("12-month", "6-month", "3-month"), name = "Discontinuation", guide = guide_legend(reverse = TRUE, row = 1)) +
    scale_x_continuous("Proportion of discontinuation", labels = scales::percent, limits = c(0, 0.4)) +
    ggtitle("Smoking") +
    facet_wrap(~drugclass, nrow = 1, labeller = as_labeller(c("other" = "Other glucose-lowering therapies"))) +
    theme_bw() +
    guides(fill = "none") +
    theme(
      plot.title.position = "plot",
      legend.position = "bottom",
      axis.text.x = element_blank(),
      axis.title.x = element_blank(),
      strip.text.x = element_blank()
    ),
  
  # Ethnicity
  discontinuation_per_drug_combined %>%
    mutate(
      perc_12m = perc_12m-perc_6m,
      perc_6m = perc_6m-perc_3m
    ) %>%
    gather("key", "value", -drugclass, -analysis_var, -group, -variable) %>%
    mutate(
      key = factor(key, levels = c("perc_12m", "perc_6m", "perc_3m"))
    ) %>%
    filter(variable == "ethnicity") %>%
    ggplot(aes(y = group, x = value, fill = drugclass, alpha = key)) +
    geom_bar(stat = "identity", colour = "black") +
    scale_fill_manual(values = c("other" = "black", "SGLT2" = "#E69F00", "GLP1" = "#56B4E9", "SU" = "#CC79A7", "DPP4" = "#0072B2", "TZD" = "#D55E00", "MFN" = "grey"), breaks = c("other"), labels = c("Other glucose-lowering therapies"), name = " ", guide = guide_legend(reverse = FALSE)) +
    scale_alpha_discrete(limits = rev, range = c(1, 0.4), breaks = c("perc_12m", "perc_6m", "perc_3m"), labels = c("12-month", "6-month", "3-month"), name = "Discontinuation", guide = guide_legend(reverse = TRUE, row = 1)) +
    scale_x_continuous("Proportion of discontinuation", labels = scales::percent, limits = c(0, 0.4)) +
    ggtitle("Ethnicity") +
    facet_wrap(~drugclass, nrow = 1, labeller = as_labeller(c("other" = "Other glucose-lowering therapies"))) +
    theme_bw() +
    guides(fill = "none") +
    theme(
      plot.title.position = "plot",
      legend.position = "bottom",
      axis.text.x = element_blank(),
      axis.title.x = element_blank(),
      strip.text.x = element_blank()
    ),
  
  # Quintile
  discontinuation_per_drug_combined %>%
    mutate(
      perc_12m = perc_12m-perc_6m,
      perc_6m = perc_6m-perc_3m
    ) %>%
    gather("key", "value", -drugclass, -analysis_var, -group, -variable) %>%
    mutate(
      key = factor(key, levels = c("perc_12m", "perc_6m", "perc_3m"))
    ) %>%
    filter(variable == "imd2015") %>%
    ggplot(aes(y = group, x = value, fill = drugclass, alpha = key)) +
    geom_bar(stat = "identity", colour = "black") +
    scale_fill_manual(values = c("other" = "black", "SGLT2" = "#E69F00", "GLP1" = "#56B4E9", "SU" = "#CC79A7", "DPP4" = "#0072B2", "TZD" = "#D55E00", "MFN" = "grey"), breaks = c("other"), labels = c("Other glucose-lowering therapies"), name = " ", guide = guide_legend(reverse = FALSE)) +
    scale_y_discrete(limits = rev, labels = c("5 (most deprived)", "4", "3", "2", "1 (least deprived)")) +
    scale_alpha_discrete(limits = rev, range = c(1, 0.4), breaks = c("perc_12m", "perc_6m", "perc_3m"), labels = c("12-month", "6-month", "3-month"), name = "Discontinuation", guide = guide_legend(reverse = TRUE, row = 1)) +
    scale_x_continuous("Proportion of discontinuation", labels = scales::percent, limits = c(0, 0.4)) +
    ggtitle("Index of multiple deprivation") +
    facet_wrap(~drugclass, nrow = 1, labeller = as_labeller(c("other" = "Other glucose-lowering therapies"))) +
    theme_bw() +
    guides(fill = "none") +
    theme(
      plot.title.position = "plot",
      legend.position = "bottom",
      axis.text.x = element_blank(),
      axis.title.x = element_blank(),
      strip.text.x = element_blank()
    ),
  
  # BMI
  discontinuation_per_drug_combined %>%
    mutate(
      perc_12m = perc_12m-perc_6m,
      perc_6m = perc_6m-perc_3m
    ) %>%
    gather("key", "value", -drugclass, -analysis_var, -group, -variable) %>%
    mutate(
      key = factor(key, levels = c("perc_12m", "perc_6m", "perc_3m"))
    ) %>%
    filter(variable == "bmi") %>%
    ggplot(aes(y = group, x = value, fill = drugclass, alpha = key)) +
    geom_bar(stat = "identity", colour = "black") +
    scale_fill_manual(values = c("other" = "black", "SGLT2" = "#E69F00", "GLP1" = "#56B4E9", "SU" = "#CC79A7", "DPP4" = "#0072B2", "TZD" = "#D55E00", "MFN" = "grey"), breaks = c("other"), labels = c("Other glucose-lowering therapies"), name = " ", guide = guide_legend(reverse = FALSE)) +
    scale_y_discrete(limits = rev, labels = c("35+", "30-34.9", "25-29.9", "<25")) +
    scale_alpha_discrete(limits = rev, range = c(1, 0.4), breaks = c("perc_12m", "perc_6m", "perc_3m"), labels = c("12-month", "6-month", "3-month"), name = "Discontinuation", guide = guide_legend(reverse = TRUE, row = 1)) +
    scale_x_continuous("Proportion of discontinuation", labels = scales::percent, limits = c(0, 0.4)) +
    ggtitle("BMI, kg/m2") +
    facet_wrap(~drugclass, nrow = 1, labeller = as_labeller(c("other" = "Other glucose-lowering therapies"))) +
    theme_bw() +
    guides(fill = "none") +
    theme(
      plot.title.position = "plot",
      legend.position = "bottom",
      axis.text.x = element_blank(),
      axis.title.x = element_blank(),
      strip.text.x = element_blank()
    ),
  
  # HbA1c
  discontinuation_per_drug_combined %>%
    mutate(
      perc_12m = perc_12m-perc_6m,
      perc_6m = perc_6m-perc_3m
    ) %>%
    gather("key", "value", -drugclass, -analysis_var, -group, -variable) %>%
    mutate(
      key = factor(key, levels = c("perc_12m", "perc_6m", "perc_3m"))
    ) %>%
    filter(variable == "hba1c") %>%
    ggplot(aes(y = group, x = value, fill = drugclass, alpha = key)) +
    geom_bar(stat = "identity", colour = "black") +
    scale_fill_manual(values = c("other" = "black", "SGLT2" = "#E69F00", "GLP1" = "#56B4E9", "SU" = "#CC79A7", "DPP4" = "#0072B2", "TZD" = "#D55E00", "MFN" = "grey"), breaks = c("other"), labels = c("Other glucose-lowering therapies"), name = " ", guide = guide_legend(reverse = FALSE)) +
    scale_y_discrete(limits = rev, labels = c("86+", "75-86", "64-75", "53-64")) +
    scale_alpha_discrete(limits = rev, range = c(1, 0.4), breaks = c("perc_12m", "perc_6m", "perc_3m"), labels = c("12-month", "6-month", "3-month"), name = "Discontinuation", guide = guide_legend(reverse = TRUE, row = 1)) +
    scale_x_continuous("Proportion of discontinuation", labels = scales::percent, limits = c(0, 0.4)) +
    ggtitle("HbA1c, mmol/mol") +
    facet_wrap(~drugclass, nrow = 1, labeller = as_labeller(c("other" = "Other glucose-lowering therapies"))) +
    theme_bw() +
    guides(fill = "none") +
    theme(
      plot.title.position = "plot",
      legend.position = "bottom",
      axis.text.x = element_blank(),
      axis.title.x = element_blank(),
      strip.text.x = element_blank()
    ),
  
  # eGFR
  discontinuation_per_drug_combined %>%
    mutate(
      perc_12m = perc_12m-perc_6m,
      perc_6m = perc_6m-perc_3m
    ) %>%
    gather("key", "value", -drugclass, -analysis_var, -group, -variable) %>%
    mutate(
      key = factor(key, levels = c("perc_12m", "perc_6m", "perc_3m"))
    ) %>%
    filter(variable == "preegfr") %>%
    ggplot(aes(y = group, x = value, fill = drugclass, alpha = key)) +
    geom_bar(stat = "identity", colour = "black") +
    scale_fill_manual(values = c("other" = "black", "SGLT2" = "#E69F00", "GLP1" = "#56B4E9", "SU" = "#CC79A7", "DPP4" = "#0072B2", "TZD" = "#D55E00", "MFN" = "grey"), breaks = c("other"), labels = c("Other glucose-lowering therapies"), name = " ", guide = guide_legend(reverse = FALSE)) +
    scale_y_discrete(limits = rev, labels = c("90+", "60-90", "<60")) +
    scale_alpha_discrete(limits = rev, range = c(1, 0.4), breaks = c("perc_12m", "perc_6m", "perc_3m"), labels = c("12-month", "6-month", "3-month"), name = "Discontinuation", guide = guide_legend(reverse = TRUE, row = 1)) +
    scale_x_continuous("Proportion of discontinuation", labels = scales::percent, limits = c(0, 0.4)) +
    ggtitle("eGFR, ml/min per 1.73 m2") +
    facet_wrap(~drugclass, nrow = 1, labeller = as_labeller(c("other" = "Other glucose-lowering therapies"))) +
    theme_bw() +
    guides(fill = "none") +
    theme(
      plot.title.position = "plot",
      legend.position = "bottom",
      axis.text.x = element_blank(),
      axis.title.x = element_blank(),
      strip.text.x = element_blank()
    ),
  
  # Statins
  discontinuation_per_drug_combined %>%
    mutate(
      perc_12m = perc_12m-perc_6m,
      perc_6m = perc_6m-perc_3m
    ) %>%
    gather("key", "value", -drugclass, -analysis_var, -group, -variable) %>%
    mutate(
      key = factor(key, levels = c("perc_12m", "perc_6m", "perc_3m"))
    ) %>%
    filter(variable == "statins") %>%
    ggplot(aes(y = group, x = value, fill = drugclass, alpha = key)) +
    geom_bar(stat = "identity", colour = "black") +
    scale_fill_manual(values = c("other" = "black", "SGLT2" = "#E69F00", "GLP1" = "#56B4E9", "SU" = "#CC79A7", "DPP4" = "#0072B2", "TZD" = "#D55E00", "MFN" = "grey"), breaks = c("other"), labels = c("Other glucose-lowering therapies"), name = " ", guide = guide_legend(reverse = FALSE)) +
    scale_y_discrete(limits = rev, breaks = c(0, 1), labels = c("No", "Yes")) +
    scale_alpha_discrete(limits = rev, range = c(1, 0.4), breaks = c("perc_12m", "perc_6m", "perc_3m"), labels = c("12-month", "6-month", "3-month"), name = "Discontinuation", guide = guide_legend(reverse = TRUE, row = 1)) +
    scale_x_continuous("Proportion of discontinuation", labels = scales::percent, limits = c(0, 0.4)) +
    ggtitle("Statins") +
    facet_wrap(~drugclass, nrow = 1, labeller = as_labeller(c("other" = "Other glucose-lowering therapies"))) +
    theme_bw() +
    guides(fill = "none") +
    theme(
      plot.title.position = "plot",
      legend.position = "bottom",
      axis.text.x = element_blank(),
      axis.title.x = element_blank(),
      strip.text.x = element_blank()
    ),
  
  # Blood pressure medication
  discontinuation_per_drug_combined %>%
    mutate(
      perc_12m = perc_12m-perc_6m,
      perc_6m = perc_6m-perc_3m
    ) %>%
    gather("key", "value", -drugclass, -analysis_var, -group, -variable) %>%
    mutate(
      key = factor(key, levels = c("perc_12m", "perc_6m", "perc_3m"))
    ) %>%
    filter(variable == "bloodmed") %>%
    ggplot(aes(y = group, x = value, fill = drugclass, alpha = key)) +
    geom_bar(stat = "identity", colour = "black") +
    scale_fill_manual(values = c("other" = "black", "SGLT2" = "#E69F00", "GLP1" = "#56B4E9", "SU" = "#CC79A7", "DPP4" = "#0072B2", "TZD" = "#D55E00", "MFN" = "grey"), breaks = c("other"), labels = c("Other glucose-lowering therapies"), name = " ", guide = guide_legend(reverse = FALSE)) +
    scale_y_discrete(limits = rev, breaks = c(0, 1), labels = c("No", "Yes")) +
    scale_alpha_discrete(limits = rev, range = c(1, 0.4), breaks = c("perc_12m", "perc_6m", "perc_3m"), labels = c("12-month", "6-month", "3-month"), name = "Discontinuation", guide = guide_legend(reverse = TRUE, row = 1)) +
    scale_x_continuous("Proportion of discontinuation", labels = scales::percent, limits = c(0, 0.4)) +
    ggtitle("Blood pressure medication") +
    facet_wrap(~drugclass, nrow = 1, labeller = as_labeller(c("other" = "Other glucose-lowering therapies"))) +
    theme_bw() +
    guides(fill = "none") +
    theme(
      plot.title.position = "plot",
      legend.position = "bottom",
      axis.text.x = element_blank(),
      axis.title.x = element_blank(),
      strip.text.x = element_blank()
    ),
  
  # Frailty
  discontinuation_per_drug_combined %>%
    mutate(
      perc_12m = perc_12m-perc_6m,
      perc_6m = perc_6m-perc_3m
    ) %>%
    gather("key", "value", -drugclass, -analysis_var, -group, -variable) %>%
    mutate(
      key = factor(key, levels = c("perc_12m", "perc_6m", "perc_3m"))
    ) %>%
    filter(variable == "frailty") %>%
    ggplot(aes(y = group, x = value, fill = drugclass, alpha = key)) +
    geom_bar(stat = "identity", colour = "black") +
    scale_fill_manual(values = c("other" = "black", "SGLT2" = "#E69F00", "GLP1" = "#56B4E9", "SU" = "#CC79A7", "DPP4" = "#0072B2", "TZD" = "#D55E00", "MFN" = "grey"), breaks = c("other"), labels = c("Other glucose-lowering therapies"), name = " ", guide = guide_legend(reverse = FALSE)) +
    scale_y_discrete(limits = rev, breaks = c(0, 1), labels = c("No", "Yes")) +
    scale_alpha_discrete(limits = rev, range = c(1, 0.4), breaks = c("perc_12m", "perc_6m", "perc_3m"), labels = c("12-month", "6-month", "3-month"), name = "Discontinuation", guide = guide_legend(reverse = TRUE, row = 1)) +
    scale_x_continuous("Proportion of discontinuation", labels = scales::percent, limits = c(0, 0.4)) +
    ggtitle("Frailty") +
    facet_wrap(~drugclass, nrow = 1, labeller = as_labeller(c("other" = "Other glucose-lowering therapies"))) +
    theme_bw() +
    guides(fill = "none") +
    theme(
      plot.title.position = "plot",
      legend.position = "bottom",
      axis.text.x = element_blank(),
      axis.title.x = element_blank(),
      strip.text.x = element_blank()
    ),
  
  # stopping metformin
  discontinuation_per_drug_combined %>%
    mutate(
      perc_12m = perc_12m-perc_6m,
      perc_6m = perc_6m-perc_3m
    ) %>%
    gather("key", "value", -drugclass, -analysis_var, -group, -variable) %>%
    mutate(
      key = factor(key, levels = c("perc_12m", "perc_6m", "perc_3m"))
    ) %>%
    filter(variable == "stopdrug_MFN") %>%
    ggplot(aes(y = group, x = value, fill = drugclass, alpha = key)) +
    geom_bar(stat = "identity", colour = "black") +
    scale_fill_manual(values = c("other" = "black", "SGLT2" = "#E69F00", "GLP1" = "#56B4E9", "SU" = "#CC79A7", "DPP4" = "#0072B2", "TZD" = "#D55E00", "MFN" = "grey"), breaks = c("other"), labels = c("Other glucose-lowering therapies"), name = " ", guide = guide_legend(reverse = FALSE)) +
    scale_y_discrete(limits = rev, breaks = c(0, 1), labels = c("No", "Yes")) +
    scale_alpha_discrete(limits = rev, range = c(1, 0.4), breaks = c("perc_12m", "perc_6m", "perc_3m"), labels = c("12-month", "6-month", "3-month"), name = "Discontinuation", guide = guide_legend(reverse = TRUE, row = 1)) +
    scale_x_continuous("Proportion of discontinuation", labels = scales::percent, limits = c(0, 0.4)) +
    ggtitle("History of Metformin discontinuation within 3-months") +
    facet_wrap(~drugclass, nrow = 1, labeller = as_labeller(c("other" = "Other glucose-lowering therapies"))) +
    theme_bw() +
    guides(fill = "none") +
    theme(
      plot.title.position = "plot",
      legend.position = "bottom",
      axis.text.x = element_blank(),
      axis.title.x = element_blank(),
      strip.text.x = element_blank()
    ),
  
  # Numdrugs
  discontinuation_per_drug_combined %>%
    mutate(
      perc_12m = perc_12m-perc_6m,
      perc_6m = perc_6m-perc_3m
    ) %>%
    gather("key", "value", -drugclass, -analysis_var, -group, -variable) %>%
    mutate(
      key = factor(key, levels = c("perc_12m", "perc_6m", "perc_3m"))
    ) %>%
    filter(variable == "numdrugs") %>%
    ggplot(aes(y = group, x = value, fill = drugclass, alpha = key)) +
    geom_bar(stat = "identity", colour = "black") +
    scale_fill_manual(values = c("other" = "black", "SGLT2" = "#E69F00", "GLP1" = "#56B4E9", "SU" = "#CC79A7", "DPP4" = "#0072B2", "TZD" = "#D55E00", "MFN" = "grey"), breaks = c("other"), labels = c("Other glucose-lowering therapies"), name = " ", guide = guide_legend(reverse = FALSE)) +
    scale_y_discrete(limits = rev, labels = c("2+", "1", "0")) +
    scale_alpha_discrete(limits = rev, range = c(1, 0.4), breaks = c("perc_12m", "perc_6m", "perc_3m"), labels = c("12-month", "6-month", "3-month"), name = "Discontinuation", guide = guide_legend(reverse = TRUE, row = 1)) +
    scale_x_continuous("Proportion of discontinuation", labels = scales::percent, limits = c(0, 0.4)) +
    ggtitle("Number of other current glucose-lowering therapies") +
    facet_wrap(~drugclass, nrow = 1, labeller = as_labeller(c("other" = "Other glucose-lowering therapies"))) +
    theme_bw() +
    guides(fill = "none") +
    theme(
      plot.title.position = "plot",
      legend.position = "bottom",
      axis.text.x = element_blank(),
      axis.title.x = element_blank(),
      strip.text.x = element_blank()
    ),
  
  # Drugline
  discontinuation_per_drug_combined %>%
    mutate(
      perc_12m = perc_12m-perc_6m,
      perc_6m = perc_6m-perc_3m
    ) %>%
    gather("key", "value", -drugclass, -analysis_var, -group, -variable) %>%
    mutate(
      key = factor(key, levels = c("perc_12m", "perc_6m", "perc_3m"))
    ) %>%
    filter(variable == "drugline") %>%
    ggplot(aes(y = group, x = value, fill = drugclass, alpha = key)) +
    geom_bar(stat = "identity", colour = "black") +
    scale_fill_manual(values = c("other" = "black", "SGLT2" = "#E69F00", "GLP1" = "#56B4E9", "SU" = "#CC79A7", "DPP4" = "#0072B2", "TZD" = "#D55E00", "MFN" = "grey"), breaks = c("other"), labels = c("Other glucose-lowering therapies"), name = " ", guide = guide_legend(reverse = FALSE)) +
    scale_y_discrete(limits = rev) +
    scale_alpha_discrete(limits = rev, range = c(1, 0.4), breaks = c("perc_12m", "perc_6m", "perc_3m"), labels = c("12-month", "6-month", "3-month"), name = "Discontinuation", guide = guide_legend(reverse = TRUE, row = 1)) +
    scale_x_continuous("Proportion of discontinuation", labels = scales::percent, limits = c(0, 0.4)) +
    ggtitle("Number of therapy classes ever prescribed") +
    facet_wrap(~drugclass, nrow = 1, labeller = as_labeller(c("other" = "Other glucose-lowering therapies"))) +
    theme_bw() +
    guides(fill = "none") +
    theme(
      plot.title.position = "plot",
      legend.position = "bottom",
      axis.text.x = element_text(size = 10),
      axis.title.x = element_text(size = 11),
      strip.text.x = element_blank()
    )
  
  , ncol = 1
  
) + 
  plot_layout(
    guides = "collect",
    height = c(1, 4, 5, 2, 3, 4, 5, 4, 4, 3, 2, 2, 2, 2, 3, 3),
    # axis_titles = "collect" # not working, manually done in the code
  ) +
  plot_layout(
    axis_titles = "collect"
  ) &
  theme(
    legend.position = "bottom",
    plot.title.position = "plot",
    legend.box = 'vertical',
    axis.title.y = element_blank(),
    axis.text.y = element_text(size = 10),
    panel.spacing = unit(0.3, "cm", data = NULL)
  )


dev.off()



pdf("results/figures/03.discontinuation_by_variable_by_drug.pdf", width = 14, height = 15)

patchwork::wrap_plots(
  
  # Age
  discontinuation_per_drug_not_combined %>%
    mutate(
      perc_12m = perc_12m-perc_6m,
      perc_6m = perc_6m-perc_3m
    ) %>%
    gather("key", "value", -drugclass, -analysis_var, -group, -variable) %>%
    mutate(
      key = factor(key, levels = c("perc_12m", "perc_6m", "perc_3m"))
    ) %>%
    filter(variable == "dstartdate_age") %>%
    ggplot(aes(y = group, x = value, fill = drugclass, alpha = key)) +
    geom_bar(stat = "identity", colour = "black") +
    scale_fill_manual(values = c("other" = "black", "SGLT2" = "#E69F00", "GLP1" = "#56B4E9", "SU" = "#CC79A7", "DPP4" = "#0072B2", "TZD" = "#D55E00", "MFN" = "grey"), breaks = c("other"), labels = c("Other glucose-lowering therapies"), name = " ", guide = guide_legend(reverse = FALSE)) +
    scale_y_discrete(limits = rev, labels = c("70+", "60-69", "50-59", "<50")) +
    scale_alpha_discrete(limits = rev, range = c(1, 0.4), breaks = c("perc_12m", "perc_6m", "perc_3m"), labels = c("12-month", "6-month", "3-month"), name = "Discontinuation", guide = guide_legend(reverse = TRUE, row = 1)) +
    scale_x_continuous("Proportion of discontinuation", labels = scales::percent, limits = c(0, 0.6)) +
    ggtitle("Age group, years") +
    facet_wrap(~drugclass, nrow = 1, labeller = as_labeller(c("GLP1" = "GLP-1RA", "DPP4" = "DPP4i", "SGLT2" = "SGLT2i", "TZD" = "TZD", "SU" = "SU"))) +
    theme_bw() +
    guides(fill = "none") +
    theme(
      plot.title.position = "plot",
      legend.position = "bottom",
      axis.text.x = element_blank(),
      axis.title.x = element_blank()
    ),
  
  # Diabetes duration
  discontinuation_per_drug_not_combined %>%
    mutate(
      perc_12m = perc_12m-perc_6m,
      perc_6m = perc_6m-perc_3m
    ) %>%
    gather("key", "value", -drugclass, -analysis_var, -group, -variable) %>%
    mutate(
      key = factor(key, levels = c("perc_12m", "perc_6m", "perc_3m"))
    ) %>%
    filter(variable == "dstartdate_dm_dur") %>%
    ggplot(aes(y = group, x = value, fill = drugclass, alpha = key)) +
    geom_bar(stat = "identity", colour = "black") +
    scale_fill_manual(values = c("other" = "black", "SGLT2" = "#E69F00", "GLP1" = "#56B4E9", "SU" = "#CC79A7", "DPP4" = "#0072B2", "TZD" = "#D55E00", "MFN" = "grey"), breaks = c("other"), labels = c("Other glucose-lowering therapies"), name = " ", guide = guide_legend(reverse = FALSE)) +
    scale_y_discrete(limits = rev, labels = c("10+", "6-9", "3-5", "1-2", "<1")) +
    scale_alpha_discrete(limits = rev, range = c(1, 0.4), breaks = c("perc_12m", "perc_6m", "perc_3m"), labels = c("12-month", "6-month", "3-month"), name = "Discontinuation", guide = guide_legend(reverse = TRUE, row = 1)) +
    scale_x_continuous("Proportion of discontinuation", labels = scales::percent, limits = c(0, 0.6)) +
    ggtitle("Duration of diabetes, years") +
    facet_wrap(~drugclass, nrow = 1, labeller = as_labeller(c("GLP1" = "GLP-1RA", "DPP4" = "DPP4i", "SGLT2" = "SGLT2i", "TZD" = "TZD", "SU" = "SU"))) +
    theme_bw() +
    guides(fill = "none") +
    theme(
      plot.title.position = "plot",
      legend.position = "bottom",
      axis.text.x = element_blank(),
      axis.title.x = element_blank(),
      strip.text.x = element_blank()
    ),
  
  # Sex
  discontinuation_per_drug_not_combined %>%
    mutate(
      perc_12m = perc_12m-perc_6m,
      perc_6m = perc_6m-perc_3m
    ) %>%
    gather("key", "value", -drugclass, -analysis_var, -group, -variable) %>%
    mutate(
      key = factor(key, levels = c("perc_12m", "perc_6m", "perc_3m"))
    ) %>%
    filter(variable == "gender") %>%
    ggplot(aes(y = group, x = value, fill = drugclass, alpha = key)) +
    geom_bar(stat = "identity", colour = "black") +
    scale_fill_manual(values = c("other" = "black", "SGLT2" = "#E69F00", "GLP1" = "#56B4E9", "SU" = "#CC79A7", "DPP4" = "#0072B2", "TZD" = "#D55E00", "MFN" = "grey"), breaks = c("other"), labels = c("Other glucose-lowering therapies"), name = " ", guide = guide_legend(reverse = FALSE)) +
    scale_y_discrete(limits = rev, labels = c("Female", "Male")) +
    scale_alpha_discrete(limits = rev, range = c(1, 0.4), breaks = c("perc_12m", "perc_6m", "perc_3m"), labels = c("12-month", "6-month", "3-month"), name = "Discontinuation", guide = guide_legend(reverse = TRUE, row = 1)) +
    scale_x_continuous("Proportion of discontinuation", labels = scales::percent, limits = c(0, 0.6)) +
    ggtitle("Sex") +
    facet_wrap(~drugclass, nrow = 1, labeller = as_labeller(c("GLP1" = "GLP-1RA", "DPP4" = "DPP4i", "SGLT2" = "SGLT2i", "TZD" = "TZD", "SU" = "SU"))) +
    theme_bw() +
    guides(fill = "none") +
    theme(
      plot.title.position = "plot",
      legend.position = "bottom",
      axis.text.x = element_blank(),
      axis.title.x = element_blank(),
      strip.text.x = element_blank()
    ),
  
  # Smoking
  discontinuation_per_drug_not_combined %>%
    mutate(
      perc_12m = perc_12m-perc_6m,
      perc_6m = perc_6m-perc_3m
    ) %>%
    gather("key", "value", -drugclass, -analysis_var, -group, -variable) %>%
    mutate(
      key = factor(key, levels = c("perc_12m", "perc_6m", "perc_3m"))
    ) %>%
    filter(variable == "smoking") %>%
    ggplot(aes(y = group, x = value, fill = drugclass, alpha = key)) +
    geom_bar(stat = "identity", colour = "black") +
    scale_fill_manual(values = c("other" = "black", "SGLT2" = "#E69F00", "GLP1" = "#56B4E9", "SU" = "#CC79A7", "DPP4" = "#0072B2", "TZD" = "#D55E00", "MFN" = "grey"), breaks = c("other"), labels = c("Other glucose-lowering therapies"), name = " ", guide = guide_legend(reverse = FALSE)) +
    scale_alpha_discrete(limits = rev, range = c(1, 0.4), breaks = c("perc_12m", "perc_6m", "perc_3m"), labels = c("12-month", "6-month", "3-month"), name = "Discontinuation", guide = guide_legend(reverse = TRUE, row = 1)) +
    scale_x_continuous("Proportion of discontinuation", labels = scales::percent, limits = c(0, 0.6)) +
    ggtitle("Smoking") +
    facet_wrap(~drugclass, nrow = 1, labeller = as_labeller(c("GLP1" = "GLP-1RA", "DPP4" = "DPP4i", "SGLT2" = "SGLT2i", "TZD" = "TZD", "SU" = "SU"))) +
    theme_bw() +
    guides(fill = "none") +
    theme(
      plot.title.position = "plot",
      legend.position = "bottom",
      axis.text.x = element_blank(),
      axis.title.x = element_blank(),
      strip.text.x = element_blank()
    ),
  
  # Ethnicity
  discontinuation_per_drug_not_combined %>%
    mutate(
      perc_12m = perc_12m-perc_6m,
      perc_6m = perc_6m-perc_3m
    ) %>%
    gather("key", "value", -drugclass, -analysis_var, -group, -variable) %>%
    mutate(
      key = factor(key, levels = c("perc_12m", "perc_6m", "perc_3m"))
    ) %>%
    filter(variable == "ethnicity") %>%
    ggplot(aes(y = group, x = value, fill = drugclass, alpha = key)) +
    geom_bar(stat = "identity", colour = "black") +
    scale_fill_manual(values = c("other" = "black", "SGLT2" = "#E69F00", "GLP1" = "#56B4E9", "SU" = "#CC79A7", "DPP4" = "#0072B2", "TZD" = "#D55E00", "MFN" = "grey"), breaks = c("other"), labels = c("Other glucose-lowering therapies"), name = " ", guide = guide_legend(reverse = FALSE)) +
    scale_alpha_discrete(limits = rev, range = c(1, 0.4), breaks = c("perc_12m", "perc_6m", "perc_3m"), labels = c("12-month", "6-month", "3-month"), name = "Discontinuation", guide = guide_legend(reverse = TRUE, row = 1)) +
    scale_x_continuous("Proportion of discontinuation", labels = scales::percent, limits = c(0, 0.6)) +
    ggtitle("Ethnicity") +
    facet_wrap(~drugclass, nrow = 1, labeller = as_labeller(c("GLP1" = "GLP-1RA", "DPP4" = "DPP4i", "SGLT2" = "SGLT2i", "TZD" = "TZD", "SU" = "SU"))) +
    theme_bw() +
    guides(fill = "none") +
    theme(
      plot.title.position = "plot",
      legend.position = "bottom",
      axis.text.x = element_blank(),
      axis.title.x = element_blank(),
      strip.text.x = element_blank()
    ),
  
  # Quintile
  discontinuation_per_drug_not_combined %>%
    mutate(
      perc_12m = perc_12m-perc_6m,
      perc_6m = perc_6m-perc_3m
    ) %>%
    gather("key", "value", -drugclass, -analysis_var, -group, -variable) %>%
    mutate(
      key = factor(key, levels = c("perc_12m", "perc_6m", "perc_3m"))
    ) %>%
    filter(variable == "imd2015") %>%
    ggplot(aes(y = group, x = value, fill = drugclass, alpha = key)) +
    geom_bar(stat = "identity", colour = "black") +
    scale_fill_manual(values = c("other" = "black", "SGLT2" = "#E69F00", "GLP1" = "#56B4E9", "SU" = "#CC79A7", "DPP4" = "#0072B2", "TZD" = "#D55E00", "MFN" = "grey"), breaks = c("other"), labels = c("Other glucose-lowering therapies"), name = " ", guide = guide_legend(reverse = FALSE)) +
    scale_y_discrete(limits = rev, labels = c("5 (most deprived)", "4", "3", "2", "1 (least deprived)")) +
    scale_alpha_discrete(limits = rev, range = c(1, 0.4), breaks = c("perc_12m", "perc_6m", "perc_3m"), labels = c("12-month", "6-month", "3-month"), name = "Discontinuation", guide = guide_legend(reverse = TRUE, row = 1)) +
    scale_x_continuous("Proportion of discontinuation", labels = scales::percent, limits = c(0, 0.6)) +
    ggtitle("Index of multiple deprivation") +
    facet_wrap(~drugclass, nrow = 1, labeller = as_labeller(c("GLP1" = "GLP-1RA", "DPP4" = "DPP4i", "SGLT2" = "SGLT2i", "TZD" = "TZD", "SU" = "SU"))) +
    theme_bw() +
    guides(fill = "none") +
    theme(
      plot.title.position = "plot",
      legend.position = "bottom",
      axis.text.x = element_blank(),
      axis.title.x = element_blank(),
      strip.text.x = element_blank()
    ),
  
  # BMI
  discontinuation_per_drug_not_combined %>%
    mutate(
      perc_12m = perc_12m-perc_6m,
      perc_6m = perc_6m-perc_3m
    ) %>%
    gather("key", "value", -drugclass, -analysis_var, -group, -variable) %>%
    mutate(
      key = factor(key, levels = c("perc_12m", "perc_6m", "perc_3m"))
    ) %>%
    filter(variable == "bmi") %>%
    ggplot(aes(y = group, x = value, fill = drugclass, alpha = key)) +
    geom_bar(stat = "identity", colour = "black") +
    scale_fill_manual(values = c("other" = "black", "SGLT2" = "#E69F00", "GLP1" = "#56B4E9", "SU" = "#CC79A7", "DPP4" = "#0072B2", "TZD" = "#D55E00", "MFN" = "grey"), breaks = c("other"), labels = c("Other glucose-lowering therapies"), name = " ", guide = guide_legend(reverse = FALSE)) +
    scale_y_discrete(limits = rev, labels = c("35+", "30-34.9", "25-29.9", "<25")) +
    scale_alpha_discrete(limits = rev, range = c(1, 0.4), breaks = c("perc_12m", "perc_6m", "perc_3m"), labels = c("12-month", "6-month", "3-month"), name = "Discontinuation", guide = guide_legend(reverse = TRUE, row = 1)) +
    scale_x_continuous("Proportion of discontinuation", labels = scales::percent, limits = c(0, 0.6)) +
    ggtitle("BMI, kg/m2") +
    facet_wrap(~drugclass, nrow = 1, labeller = as_labeller(c("GLP1" = "GLP-1RA", "DPP4" = "DPP4i", "SGLT2" = "SGLT2i", "TZD" = "TZD", "SU" = "SU"))) +
    theme_bw() +
    guides(fill = "none") +
    theme(
      plot.title.position = "plot",
      legend.position = "bottom",
      axis.text.x = element_blank(),
      axis.title.x = element_blank(),
      strip.text.x = element_blank()
    ),
  
  # HbA1c
  discontinuation_per_drug_not_combined %>%
    mutate(
      perc_12m = perc_12m-perc_6m,
      perc_6m = perc_6m-perc_3m
    ) %>%
    gather("key", "value", -drugclass, -analysis_var, -group, -variable) %>%
    mutate(
      key = factor(key, levels = c("perc_12m", "perc_6m", "perc_3m"))
    ) %>%
    filter(variable == "hba1c") %>%
    ggplot(aes(y = group, x = value, fill = drugclass, alpha = key)) +
    geom_bar(stat = "identity", colour = "black") +
    scale_fill_manual(values = c("other" = "black", "SGLT2" = "#E69F00", "GLP1" = "#56B4E9", "SU" = "#CC79A7", "DPP4" = "#0072B2", "TZD" = "#D55E00", "MFN" = "grey"), breaks = c("other"), labels = c("Other glucose-lowering therapies"), name = " ", guide = guide_legend(reverse = FALSE)) +
    scale_y_discrete(limits = rev, labels = c("86+", "75-86", "64-75", "53-64")) +
    scale_alpha_discrete(limits = rev, range = c(1, 0.4), breaks = c("perc_12m", "perc_6m", "perc_3m"), labels = c("12-month", "6-month", "3-month"), name = "Discontinuation", guide = guide_legend(reverse = TRUE, row = 1)) +
    scale_x_continuous("Proportion of discontinuation", labels = scales::percent, limits = c(0, 0.6)) +
    ggtitle("HbA1c, mmol/mol") +
    facet_wrap(~drugclass, nrow = 1, labeller = as_labeller(c("GLP1" = "GLP-1RA", "DPP4" = "DPP4i", "SGLT2" = "SGLT2i", "TZD" = "TZD", "SU" = "SU"))) +
    theme_bw() +
    guides(fill = "none") +
    theme(
      plot.title.position = "plot",
      legend.position = "bottom",
      axis.text.x = element_blank(),
      axis.title.x = element_blank(),
      strip.text.x = element_blank()
    ),
  
  # eGFR
  discontinuation_per_drug_not_combined %>%
    mutate(
      perc_12m = perc_12m-perc_6m,
      perc_6m = perc_6m-perc_3m
    ) %>%
    gather("key", "value", -drugclass, -analysis_var, -group, -variable) %>%
    mutate(
      key = factor(key, levels = c("perc_12m", "perc_6m", "perc_3m"))
    ) %>%
    filter(variable == "preegfr") %>%
    ggplot(aes(y = group, x = value, fill = drugclass, alpha = key)) +
    geom_bar(stat = "identity", colour = "black") +
    scale_fill_manual(values = c("other" = "black", "SGLT2" = "#E69F00", "GLP1" = "#56B4E9", "SU" = "#CC79A7", "DPP4" = "#0072B2", "TZD" = "#D55E00", "MFN" = "grey"), breaks = c("other"), labels = c("Other glucose-lowering therapies"), name = " ", guide = guide_legend(reverse = FALSE)) +
    scale_y_discrete(limits = rev, labels = c("90+", "60-90", "<60")) +
    scale_alpha_discrete(limits = rev, range = c(1, 0.4), breaks = c("perc_12m", "perc_6m", "perc_3m"), labels = c("12-month", "6-month", "3-month"), name = "Discontinuation", guide = guide_legend(reverse = TRUE, row = 1)) +
    scale_x_continuous("Proportion of discontinuation", labels = scales::percent, limits = c(0, 0.6)) +
    ggtitle("eGFR, ml/min per 1.73 m2") +
    facet_wrap(~drugclass, nrow = 1, labeller = as_labeller(c("GLP1" = "GLP-1RA", "DPP4" = "DPP4i", "SGLT2" = "SGLT2i", "TZD" = "TZD", "SU" = "SU"))) +
    theme_bw() +
    guides(fill = "none") +
    theme(
      plot.title.position = "plot",
      legend.position = "bottom",
      axis.text.x = element_blank(),
      axis.title.x = element_blank(),
      strip.text.x = element_blank()
    ),
  
  # Statins
  discontinuation_per_drug_not_combined %>%
    mutate(
      perc_12m = perc_12m-perc_6m,
      perc_6m = perc_6m-perc_3m
    ) %>%
    gather("key", "value", -drugclass, -analysis_var, -group, -variable) %>%
    mutate(
      key = factor(key, levels = c("perc_12m", "perc_6m", "perc_3m"))
    ) %>%
    filter(variable == "statins") %>%
    ggplot(aes(y = group, x = value, fill = drugclass, alpha = key)) +
    geom_bar(stat = "identity", colour = "black") +
    scale_fill_manual(values = c("other" = "black", "SGLT2" = "#E69F00", "GLP1" = "#56B4E9", "SU" = "#CC79A7", "DPP4" = "#0072B2", "TZD" = "#D55E00", "MFN" = "grey"), breaks = c("other"), labels = c("Other glucose-lowering therapies"), name = " ", guide = guide_legend(reverse = FALSE)) +
    scale_y_discrete(limits = rev, breaks = c(0, 1), labels = c("No", "Yes")) +
    scale_alpha_discrete(limits = rev, range = c(1, 0.4), breaks = c("perc_12m", "perc_6m", "perc_3m"), labels = c("12-month", "6-month", "3-month"), name = "Discontinuation", guide = guide_legend(reverse = TRUE, row = 1)) +
    scale_x_continuous("Proportion of discontinuation", labels = scales::percent, limits = c(0, 0.6)) +
    ggtitle("Statins") +
    facet_wrap(~drugclass, nrow = 1, labeller = as_labeller(c("GLP1" = "GLP-1RA", "DPP4" = "DPP4i", "SGLT2" = "SGLT2i", "TZD" = "TZD", "SU" = "SU"))) +
    theme_bw() +
    guides(fill = "none") +
    theme(
      plot.title.position = "plot",
      legend.position = "bottom",
      axis.text.x = element_blank(),
      axis.title.x = element_blank(),
      strip.text.x = element_blank()
    ),
  
  # Blood pressure medication
  discontinuation_per_drug_not_combined %>%
    mutate(
      perc_12m = perc_12m-perc_6m,
      perc_6m = perc_6m-perc_3m
    ) %>%
    gather("key", "value", -drugclass, -analysis_var, -group, -variable) %>%
    mutate(
      key = factor(key, levels = c("perc_12m", "perc_6m", "perc_3m"))
    ) %>%
    filter(variable == "bloodmed") %>%
    ggplot(aes(y = group, x = value, fill = drugclass, alpha = key)) +
    geom_bar(stat = "identity", colour = "black") +
    scale_fill_manual(values = c("other" = "black", "SGLT2" = "#E69F00", "GLP1" = "#56B4E9", "SU" = "#CC79A7", "DPP4" = "#0072B2", "TZD" = "#D55E00", "MFN" = "grey"), breaks = c("other"), labels = c("Other glucose-lowering therapies"), name = " ", guide = guide_legend(reverse = FALSE)) +
    scale_y_discrete(limits = rev, breaks = c(0, 1), labels = c("No", "Yes")) +
    scale_alpha_discrete(limits = rev, range = c(1, 0.4), breaks = c("perc_12m", "perc_6m", "perc_3m"), labels = c("12-month", "6-month", "3-month"), name = "Discontinuation", guide = guide_legend(reverse = TRUE, row = 1)) +
    scale_x_continuous("Proportion of discontinuation", labels = scales::percent, limits = c(0, 0.6)) +
    ggtitle("Blood pressure medication") +
    facet_wrap(~drugclass, nrow = 1, labeller = as_labeller(c("GLP1" = "GLP-1RA", "DPP4" = "DPP4i", "SGLT2" = "SGLT2i", "TZD" = "TZD", "SU" = "SU"))) +
    theme_bw() +
    guides(fill = "none") +
    theme(
      plot.title.position = "plot",
      legend.position = "bottom",
      axis.text.x = element_blank(),
      axis.title.x = element_blank(),
      strip.text.x = element_blank()
    ),
  
  # Frailty
  discontinuation_per_drug_not_combined %>%
    mutate(
      perc_12m = perc_12m-perc_6m,
      perc_6m = perc_6m-perc_3m
    ) %>%
    gather("key", "value", -drugclass, -analysis_var, -group, -variable) %>%
    mutate(
      key = factor(key, levels = c("perc_12m", "perc_6m", "perc_3m"))
    ) %>%
    filter(variable == "frailty") %>%
    ggplot(aes(y = group, x = value, fill = drugclass, alpha = key)) +
    geom_bar(stat = "identity", colour = "black") +
    scale_fill_manual(values = c("other" = "black", "SGLT2" = "#E69F00", "GLP1" = "#56B4E9", "SU" = "#CC79A7", "DPP4" = "#0072B2", "TZD" = "#D55E00", "MFN" = "grey"), breaks = c("other"), labels = c("Other glucose-lowering therapies"), name = " ", guide = guide_legend(reverse = FALSE)) +
    scale_y_discrete(limits = rev, breaks = c(0, 1), labels = c("No", "Yes")) +
    scale_alpha_discrete(limits = rev, range = c(1, 0.4), breaks = c("perc_12m", "perc_6m", "perc_3m"), labels = c("12-month", "6-month", "3-month"), name = "Discontinuation", guide = guide_legend(reverse = TRUE, row = 1)) +
    scale_x_continuous("Proportion of discontinuation", labels = scales::percent, limits = c(0, 0.6)) +
    ggtitle("Frailty") +
    facet_wrap(~drugclass, nrow = 1, labeller = as_labeller(c("GLP1" = "GLP-1RA", "DPP4" = "DPP4i", "SGLT2" = "SGLT2i", "TZD" = "TZD", "SU" = "SU"))) +
    theme_bw() +
    guides(fill = "none") +
    theme(
      plot.title.position = "plot",
      legend.position = "bottom",
      axis.text.x = element_blank(),
      axis.title.x = element_blank(),
      strip.text.x = element_blank()
    ),
  
  # stopping metformin
  discontinuation_per_drug_not_combined %>%
    mutate(
      perc_12m = perc_12m-perc_6m,
      perc_6m = perc_6m-perc_3m
    ) %>%
    gather("key", "value", -drugclass, -analysis_var, -group, -variable) %>%
    mutate(
      key = factor(key, levels = c("perc_12m", "perc_6m", "perc_3m"))
    ) %>%
    filter(variable == "stopdrug_MFN") %>%
    ggplot(aes(y = group, x = value, fill = drugclass, alpha = key)) +
    geom_bar(stat = "identity", colour = "black") +
    scale_fill_manual(values = c("other" = "black", "SGLT2" = "#E69F00", "GLP1" = "#56B4E9", "SU" = "#CC79A7", "DPP4" = "#0072B2", "TZD" = "#D55E00", "MFN" = "grey"), breaks = c("other"), labels = c("Other glucose-lowering therapies"), name = " ", guide = guide_legend(reverse = FALSE)) +
    scale_y_discrete(limits = rev, breaks = c(0, 1), labels = c("No", "Yes")) +
    scale_alpha_discrete(limits = rev, range = c(1, 0.4), breaks = c("perc_12m", "perc_6m", "perc_3m"), labels = c("12-month", "6-month", "3-month"), name = "Discontinuation", guide = guide_legend(reverse = TRUE, row = 1)) +
    scale_x_continuous("Proportion of discontinuation", labels = scales::percent, limits = c(0, 0.6)) +
    ggtitle("History of Metformin discontinuation within 3-months") +
    facet_wrap(~drugclass, nrow = 1, labeller = as_labeller(c("GLP1" = "GLP-1RA", "DPP4" = "DPP4i", "SGLT2" = "SGLT2i", "TZD" = "TZD", "SU" = "SU"))) +
    theme_bw() +
    guides(fill = "none") +
    theme(
      plot.title.position = "plot",
      legend.position = "bottom",
      axis.text.x = element_blank(),
      axis.title.x = element_blank(),
      strip.text.x = element_blank()
    ),
  
  # Numdrugs
  discontinuation_per_drug_not_combined %>%
    mutate(
      perc_12m = perc_12m-perc_6m,
      perc_6m = perc_6m-perc_3m
    ) %>%
    gather("key", "value", -drugclass, -analysis_var, -group, -variable) %>%
    mutate(
      key = factor(key, levels = c("perc_12m", "perc_6m", "perc_3m"))
    ) %>%
    filter(variable == "numdrugs") %>%
    ggplot(aes(y = group, x = value, fill = drugclass, alpha = key)) +
    geom_bar(stat = "identity", colour = "black") +
    scale_fill_manual(values = c("other" = "black", "SGLT2" = "#E69F00", "GLP1" = "#56B4E9", "SU" = "#CC79A7", "DPP4" = "#0072B2", "TZD" = "#D55E00", "MFN" = "grey"), breaks = c("other"), labels = c("Other glucose-lowering therapies"), name = " ", guide = guide_legend(reverse = FALSE)) +
    scale_y_discrete(limits = rev, labels = c("2+", "1", "0")) +
    scale_alpha_discrete(limits = rev, range = c(1, 0.4), breaks = c("perc_12m", "perc_6m", "perc_3m"), labels = c("12-month", "6-month", "3-month"), name = "Discontinuation", guide = guide_legend(reverse = TRUE, row = 1)) +
    scale_x_continuous("Proportion of discontinuation", labels = scales::percent, limits = c(0, 0.6)) +
    ggtitle("Number of other current glucose-lowering therapies") +
    facet_wrap(~drugclass, nrow = 1, labeller = as_labeller(c("GLP1" = "GLP-1RA", "DPP4" = "DPP4i", "SGLT2" = "SGLT2i", "TZD" = "TZD", "SU" = "SU"))) +
    theme_bw() +
    guides(fill = "none") +
    theme(
      plot.title.position = "plot",
      legend.position = "bottom",
      axis.text.x = element_blank(),
      axis.title.x = element_blank(),
      strip.text.x = element_blank()
    ),
  
  # Drugline
  discontinuation_per_drug_not_combined %>%
    mutate(
      perc_12m = perc_12m-perc_6m,
      perc_6m = perc_6m-perc_3m
    ) %>%
    gather("key", "value", -drugclass, -analysis_var, -group, -variable) %>%
    mutate(
      key = factor(key, levels = c("perc_12m", "perc_6m", "perc_3m"))
    ) %>%
    filter(variable == "drugline") %>%
    ggplot(aes(y = group, x = value, fill = drugclass, alpha = key)) +
    geom_bar(stat = "identity", colour = "black") +
    scale_fill_manual(values = c("other" = "black", "SGLT2" = "#E69F00", "GLP1" = "#56B4E9", "SU" = "#CC79A7", "DPP4" = "#0072B2", "TZD" = "#D55E00", "MFN" = "grey"), breaks = c("other"), labels = c("Other glucose-lowering therapies"), name = " ", guide = guide_legend(reverse = FALSE)) +
    scale_y_discrete(limits = rev) +
    scale_alpha_discrete(limits = rev, range = c(1, 0.4), breaks = c("perc_12m", "perc_6m", "perc_3m"), labels = c("12-month", "6-month", "3-month"), name = "Discontinuation", guide = guide_legend(reverse = TRUE, row = 1)) +
    scale_x_continuous("Proportion of discontinuation", labels = scales::percent, limits = c(0, 0.6)) +
    ggtitle("Number of therapy classes ever prescribed") +
    facet_wrap(~drugclass, nrow = 1, labeller = as_labeller(c("GLP1" = "GLP-1RA", "DPP4" = "DPP4i", "SGLT2" = "SGLT2i", "TZD" = "TZD", "SU" = "SU"))) +
    theme_bw() +
    guides(fill = "none") +
    theme(
      plot.title.position = "plot",
      legend.position = "bottom",
      axis.text.x = element_text(size = 10),
      axis.title.x = element_text(size = 11),
      strip.text.x = element_blank()
    )
  
  , ncol = 1
  
) + 
  plot_layout(
    guides = "collect",
    height = c(4, 5, 2, 3, 4, 5, 4, 4, 3, 2, 2, 2, 2, 3, 3),
    # axis_titles = "collect" # not working, manually done in the code
  ) +
  plot_layout(
    axis_titles = "collect"
  ) &
  theme(
    legend.position = "bottom",
    plot.title.position = "plot",
    axis.title.y = element_blank(),
    axis.text.y = element_text(size = 10),
    panel.spacing = unit(0.3, "cm", data = NULL)
  )


dev.off()









###############################################################################
###############################################################################


# Start plot by variables


pdf("results/figures/03.3m_discontinuation_by_variable.pdf", width = 7, height = 15)

patchwork::wrap_plots(
  
  # Overall
  discontinuation_per_drug_combined %>%
    gather("key", "value", -drugclass, -analysis_var, -group, -variable) %>%
    mutate(
      key = factor(key, levels = c("perc_12m", "perc_6m", "perc_3m"))
    ) %>%
    filter(variable == "overall" & key == "perc_3m") %>%
    ggplot(aes(y = analysis_var, x = value, fill = drugclass)) +
    geom_col() +
    scale_fill_manual(values = c("other" = "black", "SGLT2" = "#E69F00", "GLP1" = "#56B4E9", "SU" = "#CC79A7", "DPP4" = "#0072B2", "TZD" = "#D55E00", "MFN" = "grey"), breaks = c("other"), labels = c("Other glucose-lowering therapies"), name = " ", guide = guide_legend(reverse = FALSE)) +
    scale_y_discrete(breaks = c(1), labels = c("Overall")) +
    scale_x_continuous("Proportion of discontinuation", labels = scales::percent, limits = c(0, 0.3)) +
    facet_wrap(~drugclass, nrow = 1, labeller = as_labeller(c("other" = "Other glucose-lowering therapies"))) +
    theme_bw() +
    guides(fill = guide_legend(nrow = 1, byrow = TRUE)) +
    theme(
      legend.position = "bottom",
      axis.text.x = element_blank(),
      axis.title.x = element_blank()
    ),
  
  # Age
  discontinuation_per_drug_combined %>%
    gather("key", "value", -drugclass, -analysis_var, -group, -variable) %>%
    mutate(
      key = factor(key, levels = c("perc_12m", "perc_6m", "perc_3m"))
    ) %>%
    filter(variable == "dstartdate_age" & key == "perc_3m") %>%
    ggplot(aes(y = group, x = value, fill = drugclass)) +
    geom_col() +
    scale_fill_manual(values = c("other" = "black", "SGLT2" = "#E69F00", "GLP1" = "#56B4E9", "SU" = "#CC79A7", "DPP4" = "#0072B2", "TZD" = "#D55E00", "MFN" = "grey"), breaks = c("other"), labels = c("Other glucose-lowering therapies"), name = " ", guide = guide_legend(reverse = FALSE)) +
    scale_y_discrete(limits = rev, labels = c("70+", "60-69", "50-59", "<50")) +
    scale_x_continuous("Proportion of discontinuation", labels = scales::percent, limits = c(0, 0.3)) +
    ggtitle("Age group, years") +
    facet_wrap(~drugclass, nrow = 1, labeller = as_labeller(c("other" = "Other glucose-lowering therapies"))) +
    theme_bw() +
    guides(fill = guide_legend(nrow = 1, byrow = TRUE)) +
    theme(
      plot.title.position = "plot",
      legend.position = "bottom",
      axis.text.x = element_blank(),
      axis.title.x = element_blank(),
      strip.text.x = element_blank()
    ),
  
  # Diabetes duration
  discontinuation_per_drug_combined %>%
    gather("key", "value", -drugclass, -analysis_var, -group, -variable) %>%
    mutate(
      key = factor(key, levels = c("perc_12m", "perc_6m", "perc_3m"))
    ) %>%
    filter(variable == "dstartdate_dm_dur" & key == "perc_3m") %>%
    ggplot(aes(y = group, x = value, fill = drugclass)) +
    geom_col() +
    scale_fill_manual(values = c("other" = "black", "SGLT2" = "#E69F00", "GLP1" = "#56B4E9", "SU" = "#CC79A7", "DPP4" = "#0072B2", "TZD" = "#D55E00", "MFN" = "grey"), breaks = c("other"), labels = c("Other glucose-lowering therapies"), name = " ", guide = guide_legend(reverse = FALSE)) +
    scale_y_discrete(limits = rev, labels = c("10+", "6-9", "3-5", "1-2", "<1")) +
    ggtitle("Duration of diabetes, years") +
    scale_x_continuous("Proportion of discontinuation", labels = scales::percent, limits = c(0, 0.3)) +
    facet_wrap(~drugclass, nrow = 1, labeller = as_labeller(c("other" = "Other glucose-lowering therapies"))) +
    theme_bw() +
    guides(fill = guide_legend(nrow = 1, byrow = TRUE)) +
    theme(
      plot.title.position = "plot",
      legend.position = "bottom",
      axis.text.x = element_blank(),
      axis.title.x = element_blank(),
      strip.text.x = element_blank()
    ),
  
  # Sex
  discontinuation_per_drug_combined %>%
    gather("key", "value", -drugclass, -analysis_var, -group, -variable) %>%
    mutate(
      key = factor(key, levels = c("perc_12m", "perc_6m", "perc_3m"))
    ) %>%
    filter(variable == "gender" & key == "perc_3m") %>%
    ggplot(aes(y = group, x = value, fill = drugclass)) +
    geom_col() +
    scale_fill_manual(values = c("other" = "black", "SGLT2" = "#E69F00", "GLP1" = "#56B4E9", "SU" = "#CC79A7", "DPP4" = "#0072B2", "TZD" = "#D55E00", "MFN" = "grey"), breaks = c("other"), labels = c("Other glucose-lowering therapies"), name = " ", guide = guide_legend(reverse = FALSE)) +
    scale_y_discrete(limits = rev, labels = c("Female", "Male")) +
    scale_x_continuous("Proportion of discontinuation", labels = scales::percent, limits = c(0, 0.3)) +
    ggtitle("Sex") +
    facet_wrap(~drugclass, nrow = 1, labeller = as_labeller(c("other" = "Other glucose-lowering therapies"))) +
    theme_bw() +
    guides(fill = guide_legend(nrow = 1, byrow = TRUE)) +
    theme(
      plot.title.position = "plot",
      legend.position = "bottom",
      axis.text.x = element_blank(),
      axis.title.x = element_blank(),
      strip.text.x = element_blank()
    ),
  
  # Smoking
  discontinuation_per_drug_combined %>%
    gather("key", "value", -drugclass, -analysis_var, -group, -variable) %>%
    mutate(
      key = factor(key, levels = c("perc_12m", "perc_6m", "perc_3m"))
    ) %>%
    filter(variable == "smoking" & key == "perc_3m") %>%
    ggplot(aes(y = group, x = value, fill = drugclass)) +
    geom_col() +
    scale_fill_manual(values = c("other" = "black", "SGLT2" = "#E69F00", "GLP1" = "#56B4E9", "SU" = "#CC79A7", "DPP4" = "#0072B2", "TZD" = "#D55E00", "MFN" = "grey"), breaks = c("other"), labels = c("Other glucose-lowering therapies"), name = " ", guide = guide_legend(reverse = FALSE)) +
    ggtitle("Smoking") +
    scale_x_continuous("Proportion of discontinuation", labels = scales::percent, limits = c(0, 0.3)) +
    facet_wrap(~drugclass, nrow = 1, labeller = as_labeller(c("other" = "Other glucose-lowering therapies"))) +
    theme_bw() +
    guides(fill = guide_legend(nrow = 1, byrow = TRUE)) +
    theme(
      plot.title.position = "plot",
      legend.position = "bottom",
      axis.text.x = element_blank(),
      axis.title.x = element_blank(),
      strip.text.x = element_blank()
    ),
  
  # Ethnicity
  discontinuation_per_drug_combined %>%
    gather("key", "value", -drugclass, -analysis_var, -group, -variable) %>%
    mutate(
      key = factor(key, levels = c("perc_12m", "perc_6m", "perc_3m"))
    ) %>%
    filter(variable == "ethnicity" & key == "perc_3m") %>%
    ggplot(aes(y = group, x = value, fill = drugclass)) +
    geom_col() +
    scale_fill_manual(values = c("other" = "black", "SGLT2" = "#E69F00", "GLP1" = "#56B4E9", "SU" = "#CC79A7", "DPP4" = "#0072B2", "TZD" = "#D55E00", "MFN" = "grey"), breaks = c("other"), labels = c("Other glucose-lowering therapies"), name = " ", guide = guide_legend(reverse = FALSE)) +
    ggtitle("Ethnicity") +
    scale_x_continuous("Proportion of discontinuation", labels = scales::percent, limits = c(0, 0.3)) +
    facet_wrap(~drugclass, nrow = 1, labeller = as_labeller(c("other" = "Other glucose-lowering therapies"))) +
    theme_bw() +
    guides(fill = guide_legend(nrow = 1, byrow = TRUE)) +
    theme(
      plot.title.position = "plot",
      legend.position = "bottom",
      axis.text.x = element_blank(),
      axis.title.x = element_blank(),
      strip.text.x = element_blank()
    ),
  
  # Quintile
  discontinuation_per_drug_combined %>%
    gather("key", "value", -drugclass, -analysis_var, -group, -variable) %>%
    mutate(
      key = factor(key, levels = c("perc_12m", "perc_6m", "perc_3m"))
    ) %>%
    filter(variable == "imd2015" & key == "perc_3m") %>%
    ggplot(aes(y = group, x = value, fill = drugclass)) +
    geom_col() +
    scale_fill_manual(values = c("other" = "black", "SGLT2" = "#E69F00", "GLP1" = "#56B4E9", "SU" = "#CC79A7", "DPP4" = "#0072B2", "TZD" = "#D55E00", "MFN" = "grey"), breaks = c("other"), labels = c("Other glucose-lowering therapies"), name = " ", guide = guide_legend(reverse = FALSE)) +
    scale_y_discrete(limits = rev, labels = c("5 (most deprived)", "4", "3", "2", "1 (least deprived)")) +
    scale_x_continuous("Proportion of discontinuation", labels = scales::percent, limits = c(0, 0.3)) +
    ggtitle("Index of multiple deprivation") +
    facet_wrap(~drugclass, nrow = 1, labeller = as_labeller(c("other" = "Other glucose-lowering therapies"))) +
    theme_bw() +
    guides(fill = guide_legend(nrow = 1, byrow = TRUE)) +
    theme(
      plot.title.position = "plot",
      legend.position = "bottom",
      axis.text.x = element_blank(),
      axis.title.x = element_blank(),
      strip.text.x = element_blank()
    ),
  
  # BMI
  discontinuation_per_drug_combined %>%
    gather("key", "value", -drugclass, -analysis_var, -group, -variable) %>%
    mutate(
      key = factor(key, levels = c("perc_12m", "perc_6m", "perc_3m"))
    ) %>%
    filter(variable == "bmi" & key == "perc_3m") %>%
    ggplot(aes(y = group, x = value, fill = drugclass)) +
    geom_col() +
    scale_fill_manual(values = c("other" = "black", "SGLT2" = "#E69F00", "GLP1" = "#56B4E9", "SU" = "#CC79A7", "DPP4" = "#0072B2", "TZD" = "#D55E00", "MFN" = "grey"), breaks = c("other"), labels = c("Other glucose-lowering therapies"), name = " ", guide = guide_legend(reverse = FALSE)) +
    scale_y_discrete(limits = rev, labels = c("35+", "30-34.9", "25-29.9", "<25")) +
    ggtitle("BMI, kg/m2") +
    scale_x_continuous("Proportion of discontinuation", labels = scales::percent, limits = c(0, 0.3)) +
    facet_wrap(~drugclass, nrow = 1, labeller = as_labeller(c("other" = "Other glucose-lowering therapies"))) +
    theme_bw() +
    guides(fill = guide_legend(nrow = 1, byrow = TRUE)) +
    theme(
      plot.title.position = "plot",
      legend.position = "bottom",
      axis.text.x = element_blank(),
      axis.title.x = element_blank(),
      strip.text.x = element_blank()
    ),
  
  # HbA1c
  discontinuation_per_drug_combined %>%
    gather("key", "value", -drugclass, -analysis_var, -group, -variable) %>%
    mutate(
      key = factor(key, levels = c("perc_12m", "perc_6m", "perc_3m"))
    ) %>%
    filter(variable == "hba1c" & key == "perc_3m") %>%
    ggplot(aes(y = group, x = value, fill = drugclass)) +
    geom_col() +
    scale_fill_manual(values = c("other" = "black", "SGLT2" = "#E69F00", "GLP1" = "#56B4E9", "SU" = "#CC79A7", "DPP4" = "#0072B2", "TZD" = "#D55E00", "MFN" = "grey"), breaks = c("other"), labels = c("Other glucose-lowering therapies"), name = " ", guide = guide_legend(reverse = FALSE)) +
    scale_y_discrete(limits = rev, labels = c("86+", "75-86", "64-75", "53-64")) +
    ggtitle("HbA1c, mmol/mol") +
    scale_x_continuous("Proportion of discontinuation", labels = scales::percent, limits = c(0, 0.3)) +
    facet_wrap(~drugclass, nrow = 1, labeller = as_labeller(c("other" = "Other glucose-lowering therapies"))) +
    theme_bw() +
    guides(fill = guide_legend(nrow = 1, byrow = TRUE)) +
    theme(
      plot.title.position = "plot",
      legend.position = "bottom",
      axis.text.x = element_blank(),
      axis.title.x = element_blank(),
      strip.text.x = element_blank()
    ),
  
  # eGFR
  discontinuation_per_drug_combined %>%
    gather("key", "value", -drugclass, -analysis_var, -group, -variable) %>%
    mutate(
      key = factor(key, levels = c("perc_12m", "perc_6m", "perc_3m"))
    ) %>%
    filter(variable == "preegfr" & key == "perc_3m") %>%
    ggplot(aes(y = group, x = value, fill = drugclass)) +
    geom_col() +
    scale_fill_manual(values = c("other" = "black", "SGLT2" = "#E69F00", "GLP1" = "#56B4E9", "SU" = "#CC79A7", "DPP4" = "#0072B2", "TZD" = "#D55E00", "MFN" = "grey"), breaks = c("other"), labels = c("Other glucose-lowering therapies"), name = " ", guide = guide_legend(reverse = FALSE)) +
    scale_y_discrete(limits = rev, labels = c("90+", "60-90", "<60")) +
    ggtitle("eGFR, ml/min per 1.73 m2") +
    scale_x_continuous("Proportion of discontinuation", labels = scales::percent, limits = c(0, 0.3)) +
    facet_wrap(~drugclass, nrow = 1, labeller = as_labeller(c("other" = "Other glucose-lowering therapies"))) +
    theme_bw() +
    guides(fill = guide_legend(nrow = 1, byrow = TRUE)) +
    theme(
      plot.title.position = "plot",
      legend.position = "bottom",
      axis.text.x = element_blank(),
      axis.title.x = element_blank(),
      strip.text.x = element_blank()
    ),
  
  # Statins
  discontinuation_per_drug_combined %>%
    gather("key", "value", -drugclass, -analysis_var, -group, -variable) %>%
    mutate(
      key = factor(key, levels = c("perc_12m", "perc_6m", "perc_3m"))
    ) %>%
    filter(variable == "statins" & key == "perc_3m") %>%
    ggplot(aes(y = group, x = value, fill = drugclass)) +
    geom_col() +
    scale_fill_manual(values = c("other" = "black", "SGLT2" = "#E69F00", "GLP1" = "#56B4E9", "SU" = "#CC79A7", "DPP4" = "#0072B2", "TZD" = "#D55E00", "MFN" = "grey"), breaks = c("other"), labels = c("Other glucose-lowering therapies"), name = " ", guide = guide_legend(reverse = FALSE)) +
    scale_y_discrete(limits = rev, breaks = c(0, 1), labels = c("No", "Yes")) +
    ggtitle("Statins") +
    scale_x_continuous("Proportion of discontinuation", labels = scales::percent, limits = c(0, 0.3)) +
    facet_wrap(~drugclass, nrow = 1, labeller = as_labeller(c("other" = "Other glucose-lowering therapies"))) +
    theme_bw() +
    guides(fill = guide_legend(nrow = 1, byrow = TRUE)) +
    theme(
      plot.title.position = "plot",
      legend.position = "bottom",
      axis.text.x = element_blank(),
      axis.title.x = element_blank(),
      strip.text.x = element_blank()
    ),
  
  # Bloodmed
  discontinuation_per_drug_combined %>%
    gather("key", "value", -drugclass, -analysis_var, -group, -variable) %>%
    mutate(
      key = factor(key, levels = c("perc_12m", "perc_6m", "perc_3m"))
    ) %>%
    filter(variable == "bloodmed" & key == "perc_3m") %>%
    ggplot(aes(y = group, x = value, fill = drugclass)) +
    geom_col() +
    scale_fill_manual(values = c("other" = "black", "SGLT2" = "#E69F00", "GLP1" = "#56B4E9", "SU" = "#CC79A7", "DPP4" = "#0072B2", "TZD" = "#D55E00", "MFN" = "grey"), breaks = c("other"), labels = c("Other glucose-lowering therapies"), name = " ", guide = guide_legend(reverse = FALSE)) +
    scale_y_discrete(limits = rev, breaks = c(0, 1), labels = c("No", "Yes")) +
    ggtitle("Blood pressure medication") +
    scale_x_continuous("Proportion of discontinuation", labels = scales::percent, limits = c(0, 0.3)) +
    facet_wrap(~drugclass, nrow = 1, labeller = as_labeller(c("other" = "Other glucose-lowering therapies"))) +
    theme_bw() +
    guides(fill = guide_legend(nrow = 1, byrow = TRUE)) +
    theme(
      plot.title.position = "plot",
      legend.position = "bottom",
      axis.text.x = element_blank(),
      axis.title.x = element_blank(),
      strip.text.x = element_blank()
    ),
  
  # Frailty
  discontinuation_per_drug_combined %>%
    gather("key", "value", -drugclass, -analysis_var, -group, -variable) %>%
    mutate(
      key = factor(key, levels = c("perc_12m", "perc_6m", "perc_3m"))
    ) %>%
    filter(variable == "frailty" & key == "perc_3m") %>%
    ggplot(aes(y = group, x = value, fill = drugclass)) +
    geom_col() +
    scale_fill_manual(values = c("other" = "black", "SGLT2" = "#E69F00", "GLP1" = "#56B4E9", "SU" = "#CC79A7", "DPP4" = "#0072B2", "TZD" = "#D55E00", "MFN" = "grey"), breaks = c("other"), labels = c("Other glucose-lowering therapies"), name = " ", guide = guide_legend(reverse = FALSE)) +
    scale_y_discrete(limits = rev, breaks = c(0, 1), labels = c("No", "Yes")) +
    ggtitle("Frailty") +
    scale_x_continuous("Proportion of discontinuation", labels = scales::percent, limits = c(0, 0.3)) +
    facet_wrap(~drugclass, nrow = 1, labeller = as_labeller(c("other" = "Other glucose-lowering therapies"))) +
    theme_bw() +
    guides(fill = guide_legend(nrow = 1, byrow = TRUE)) +
    theme(
      plot.title.position = "plot",
      legend.position = "bottom",
      axis.text.x = element_blank(),
      axis.title.x = element_blank(),
      strip.text.x = element_blank()
    ),
  
  # stopping metformin
  discontinuation_per_drug_combined %>%
    gather("key", "value", -drugclass, -analysis_var, -group, -variable) %>%
    mutate(
      key = factor(key, levels = c("perc_12m", "perc_6m", "perc_3m"))
    ) %>%
    filter(variable == "stopdrug_MFN" & key == "perc_3m") %>%
    ggplot(aes(y = group, x = value, fill = drugclass)) +
    geom_col() +
    scale_fill_manual(values = c("other" = "black", "SGLT2" = "#E69F00", "GLP1" = "#56B4E9", "SU" = "#CC79A7", "DPP4" = "#0072B2", "TZD" = "#D55E00", "MFN" = "grey"), breaks = c("other"), labels = c("Other glucose-lowering therapies"), name = " ", guide = guide_legend(reverse = FALSE)) +
    scale_y_discrete(limits = rev, breaks = c(0, 1), labels = c("No", "Yes")) +
    ggtitle("History of Metformin discontinuation within 3-months") +
    scale_x_continuous("Proportion of discontinuation", labels = scales::percent, limits = c(0, 0.3)) +
    facet_wrap(~drugclass, nrow = 1, labeller = as_labeller(c("other" = "Other glucose-lowering therapies"))) +
    theme_bw() +
    guides(fill = guide_legend(nrow = 1, byrow = TRUE)) +
    theme(
      plot.title.position = "plot",
      legend.position = "bottom",
      axis.text.x = element_blank(),
      axis.title.x = element_blank(),
      strip.text.x = element_blank()
    ),
  
  # Numdrugs
  discontinuation_per_drug_combined %>%
    gather("key", "value", -drugclass, -analysis_var, -group, -variable) %>%
    mutate(
      key = factor(key, levels = c("perc_12m", "perc_6m", "perc_3m"))
    ) %>%
    filter(variable == "numdrugs" & key == "perc_3m") %>%
    ggplot(aes(y = group, x = value, fill = drugclass)) +
    geom_col() +
    scale_fill_manual(values = c("other" = "black", "SGLT2" = "#E69F00", "GLP1" = "#56B4E9", "SU" = "#CC79A7", "DPP4" = "#0072B2", "TZD" = "#D55E00", "MFN" = "grey"), breaks = c("other"), labels = c("Other glucose-lowering therapies"), name = " ", guide = guide_legend(reverse = FALSE)) +
    scale_y_discrete(limits = rev, labels = c("2+", "1", "0")) +
    ggtitle("Number of other current glucose-lowering therapies") +
    scale_x_continuous("Proportion of discontinuation", labels = scales::percent, limits = c(0, 0.3)) +
    facet_wrap(~drugclass, nrow = 1, labeller = as_labeller(c("other" = "Other glucose-lowering therapies"))) +
    theme_bw() +
    guides(fill = guide_legend(nrow = 1, byrow = TRUE)) +
    theme(
      plot.title.position = "plot",
      legend.position = "bottom",
      axis.text.x = element_blank(),
      axis.title.x = element_blank(),
      strip.text.x = element_blank()
    ),
  
  # drugline
  discontinuation_per_drug_combined %>%
    gather("key", "value", -drugclass, -analysis_var, -group, -variable) %>%
    mutate(
      key = factor(key, levels = c("perc_12m", "perc_6m", "perc_3m"))
    ) %>%
    filter(variable == "drugline" & key == "perc_3m") %>%
    ggplot(aes(y = group, x = value, fill = drugclass)) +
    geom_col() +
    scale_fill_manual(values = c("other" = "black", "SGLT2" = "#E69F00", "GLP1" = "#56B4E9", "SU" = "#CC79A7", "DPP4" = "#0072B2", "TZD" = "#D55E00", "MFN" = "grey"), breaks = c("other"), labels = c("Other glucose-lowering therapies"), name = " ", guide = guide_legend(reverse = FALSE)) +
    scale_y_discrete(limits = rev) +
    ggtitle("Number of therapy classes ever prescribed") +
    scale_x_continuous("Proportion of discontinuation", labels = scales::percent, limits = c(0, 0.3)) +
    facet_wrap(~drugclass, nrow = 1, labeller = as_labeller(c("other" = "Other glucose-lowering therapies"))) +
    theme_bw() +
    guides(fill = guide_legend(nrow = 1, byrow = TRUE)) +
    theme(
      plot.title.position = "plot",
      legend.position = "bottom",
      axis.text.x = element_text(size = 10),
      axis.title.x = element_text(size = 11),
      strip.text.x = element_blank()
    )
  
  , ncol = 1
  
) + 
  plot_layout(
    guides = "collect",
    height = c(1, 4, 5, 2, 3, 4, 5, 4, 4, 3, 2, 2, 2, 2, 3, 3),
    # axis_titles = "collect" # not working, manually done in the code
  ) +
  plot_layout(
    axis_titles = "collect"
  ) &
  theme(
    legend.position = "none",
    plot.title.position = "plot",
    axis.title.y = element_blank(),
    axis.text.y = element_text(size = 10),
    panel.spacing = unit(0.3, "cm", data = NULL)
  )


dev.off()


pdf("results/figures/03.3m_discontinuation_by_variable_by_drug.pdf", width = 9, height = 15)

patchwork::wrap_plots(
  
  # Age - Metformin
  discontinuation_per_drug_not_combined %>%
    gather("key", "value", -drugclass, -analysis_var, -group, -variable) %>%
    mutate(
      key = factor(key, levels = c("perc_12m", "perc_6m", "perc_3m"))
    ) %>%
    filter(variable == "dstartdate_age" & key == "perc_3m") %>%
    ggplot(aes(y = group, x = value, fill = drugclass)) +
    geom_col(position = "dodge") +
    scale_fill_manual(values = c("Pooled" = "black", "SGLT2" = "#E69F00", "GLP1" = "#56B4E9", "SU" = "#CC79A7", "DPP4" = "#0072B2", "TZD" = "#D55E00", "MFN" = "grey"), breaks = c("GLP1", "DPP4", "SGLT2", "TZD", "SU"), labels = c("GLP-1RA", "DPP4i", "SGLT2i", "TZD", "SU"), name = " ", guide = guide_legend(reverse = TRUE)) +
    scale_y_discrete(limits = rev, labels = c("70+", "60-69", "50-59", "<50")) +
    ggtitle("Age group, years") +
    scale_x_continuous("Proportion of discontinuation", labels = scales::percent, limits = c(0, 0.3), breaks = seq(0, 0.3, 0.1)) +
    facet_wrap(~drugclass, nrow = 1, labeller = as_labeller(c("GLP1" = "GLP-1RA", "DPP4" = "DPP4i", "SGLT2" = "SGLT2i", "TZD" = "TZD", "SU" = "SU"))) +
    theme_bw() +
    guides(fill = guide_legend(nrow = 1, byrow = TRUE)) +
    theme(
      legend.position = "bottom",
      axis.text.x = element_blank(),
      axis.title.x = element_blank()
    ),
  
  # Diabetes duration
  discontinuation_per_drug_not_combined %>%
    gather("key", "value", -drugclass, -analysis_var, -group, -variable) %>%
    mutate(
      key = factor(key, levels = c("perc_12m", "perc_6m", "perc_3m"))
    ) %>%
    filter(variable == "dstartdate_dm_dur" & key == "perc_3m") %>%
    ggplot(aes(y = group, x = value, fill = drugclass)) +
    geom_col(position = "dodge") +
    scale_fill_manual(values = c("Pooled" = "black", "SGLT2" = "#E69F00", "GLP1" = "#56B4E9", "SU" = "#CC79A7", "DPP4" = "#0072B2", "TZD" = "#D55E00", "MFN" = "grey"), breaks = c("GLP1", "DPP4", "SGLT2", "TZD", "SU"), labels = c("GLP-1RA", "DPP4i", "SGLT2i", "TZD", "SU"), name = " ", guide = guide_legend(reverse = TRUE)) +
    scale_y_discrete(limits = rev, labels = c("10+", "6-9", "3-5", "1-2", "<1")) +
    ggtitle("Duration of diabetes, years") +
    scale_x_continuous("Proportion of discontinuation", labels = scales::percent, limits = c(0, 0.3), breaks = seq(0, 0.3, 0.1)) +
    facet_wrap(~drugclass, nrow = 1, labeller = as_labeller(c("GLP1" = "GLP-1RA", "DPP4" = "DPP4i", "SGLT2" = "SGLT2i", "TZD" = "TZD", "SU" = "SU"))) +
    theme_bw() +
    guides(fill = guide_legend(nrow = 1, byrow = TRUE)) +
    theme(
      legend.position = "bottom",
      axis.text.x = element_blank(),
      axis.title.x = element_blank(),
      strip.text.x = element_blank()
    ),
  
  # Sex 
  discontinuation_per_drug_not_combined %>%
    gather("key", "value", -drugclass, -analysis_var, -group, -variable) %>%
    mutate(
      key = factor(key, levels = c("perc_12m", "perc_6m", "perc_3m"))
    ) %>%
    filter(variable == "gender" & key == "perc_3m") %>%
    ggplot(aes(y = group, x = value, fill = drugclass)) +
    geom_col(position = "dodge") +
    scale_fill_manual(values = c("Pooled" = "black", "SGLT2" = "#E69F00", "GLP1" = "#56B4E9", "SU" = "#CC79A7", "DPP4" = "#0072B2", "TZD" = "#D55E00", "MFN" = "grey"), breaks = c("GLP1", "DPP4", "SGLT2", "TZD", "SU"), labels = c("GLP-1RA", "DPP4i", "SGLT2i", "TZD", "SU"), name = " ", guide = guide_legend(reverse = TRUE)) +
    scale_y_discrete(limits = rev, labels = c("Female", "Male")) +
    ggtitle("Sex") +
    scale_x_continuous("Proportion of discontinuation", labels = scales::percent, limits = c(0, 0.3), breaks = seq(0, 0.3, 0.1)) +
    facet_wrap(~drugclass, nrow = 1, labeller = as_labeller(c("GLP1" = "GLP-1RA", "DPP4" = "DPP4i", "SGLT2" = "SGLT2i", "TZD" = "TZD", "SU" = "SU"))) +
    theme_bw() +
    guides(fill = guide_legend(nrow = 1, byrow = TRUE)) +
    theme(
      legend.position = "bottom",
      axis.text.x = element_blank(),
      axis.title.x = element_blank(),
      strip.text.x = element_blank()
    ),
  
  
  # Smoking 
  discontinuation_per_drug_not_combined %>%
    gather("key", "value", -drugclass, -analysis_var, -group, -variable) %>%
    mutate(
      key = factor(key, levels = c("perc_12m", "perc_6m", "perc_3m"))
    ) %>%
    filter(variable == "smoking" & key == "perc_3m") %>%
    ggplot(aes(y = group, x = value, fill = drugclass)) +
    geom_col(position = "dodge") +
    scale_fill_manual(values = c("Pooled" = "black", "SGLT2" = "#E69F00", "GLP1" = "#56B4E9", "SU" = "#CC79A7", "DPP4" = "#0072B2", "TZD" = "#D55E00", "MFN" = "grey"), breaks = c("GLP1", "DPP4", "SGLT2", "TZD", "SU"), labels = c("GLP-1RA", "DPP4i", "SGLT2i", "TZD", "SU"), name = " ", guide = guide_legend(reverse = TRUE)) +
    ggtitle("Smoking") +
    scale_x_continuous("Proportion of discontinuation", labels = scales::percent, limits = c(0, 0.3), breaks = seq(0, 0.3, 0.1)) +
    facet_wrap(~drugclass, nrow = 1, labeller = as_labeller(c("GLP1" = "GLP-1RA", "DPP4" = "DPP4i", "SGLT2" = "SGLT2i", "TZD" = "TZD", "SU" = "SU"))) +
    theme_bw() +
    guides(fill = guide_legend(nrow = 1, byrow = TRUE)) +
    theme(
      legend.position = "bottom",
      axis.text.x = element_blank(),
      axis.title.x = element_blank(),
      strip.text.x = element_blank()
    ),
  
  # Ethnicity
  discontinuation_per_drug_not_combined %>%
    gather("key", "value", -drugclass, -analysis_var, -group, -variable) %>%
    mutate(
      key = factor(key, levels = c("perc_12m", "perc_6m", "perc_3m"))
    ) %>%
    filter(variable == "ethnicity" & key == "perc_3m") %>%
    ggplot(aes(y = group, x = value, fill = drugclass)) +
    geom_col(position = "dodge") +
    scale_fill_manual(values = c("Pooled" = "black", "SGLT2" = "#E69F00", "GLP1" = "#56B4E9", "SU" = "#CC79A7", "DPP4" = "#0072B2", "TZD" = "#D55E00", "MFN" = "grey"), breaks = c("GLP1", "DPP4", "SGLT2", "TZD", "SU"), labels = c("GLP-1RA", "DPP4i", "SGLT2i", "TZD", "SU"), name = " ", guide = guide_legend(reverse = TRUE)) +
    ggtitle("Ethnicity") +
    scale_x_continuous("Proportion of discontinuation", labels = scales::percent, limits = c(0, 0.3), breaks = seq(0, 0.3, 0.1)) +
    facet_wrap(~drugclass, nrow = 1, labeller = as_labeller(c("GLP1" = "GLP-1RA", "DPP4" = "DPP4i", "SGLT2" = "SGLT2i", "TZD" = "TZD", "SU" = "SU"))) +
    theme_bw() +
    guides(fill = guide_legend(nrow = 1, byrow = TRUE)) +
    theme(
      legend.position = "bottom",
      axis.text.x = element_blank(),
      axis.title.x = element_blank(),
      strip.text.x = element_blank()
    ),
  
  # Quintile
  discontinuation_per_drug_not_combined %>%
    gather("key", "value", -drugclass, -analysis_var, -group, -variable) %>%
    mutate(
      key = factor(key, levels = c("perc_12m", "perc_6m", "perc_3m"))
    ) %>%
    filter(variable == "imd2015" & key == "perc_3m") %>%
    ggplot(aes(y = group, x = value, fill = drugclass)) +
    geom_col(position = "dodge") +
    scale_fill_manual(values = c("Pooled" = "black", "SGLT2" = "#E69F00", "GLP1" = "#56B4E9", "SU" = "#CC79A7", "DPP4" = "#0072B2", "TZD" = "#D55E00", "MFN" = "grey"), breaks = c("GLP1", "DPP4", "SGLT2", "TZD", "SU"), labels = c("GLP-1RA", "DPP4i", "SGLT2i", "TZD", "SU"), name = " ", guide = guide_legend(reverse = TRUE)) +
    scale_y_discrete(limits = rev, labels = c("5 (most deprived)", "4", "3", "2", "1 (least deprived)")) +
    ggtitle("Index of multiple deprivation") +
    scale_x_continuous("Proportion of discontinuation", labels = scales::percent, limits = c(0, 0.3), breaks = seq(0, 0.3, 0.1)) +
    facet_wrap(~drugclass, nrow = 1, labeller = as_labeller(c("GLP1" = "GLP-1RA", "DPP4" = "DPP4i", "SGLT2" = "SGLT2i", "TZD" = "TZD", "SU" = "SU"))) +
    theme_bw() +
    guides(fill = guide_legend(nrow = 1, byrow = TRUE)) +
    theme(
      legend.position = "bottom",
      axis.text.x = element_blank(),
      axis.title.x = element_blank(),
      strip.text.x = element_blank()
    ),
  
  # BMI
  discontinuation_per_drug_not_combined %>%
    gather("key", "value", -drugclass, -analysis_var, -group, -variable) %>%
    mutate(
      key = factor(key, levels = c("perc_12m", "perc_6m", "perc_3m"))
    ) %>%
    filter(variable == "bmi" & key == "perc_3m") %>%
    ggplot(aes(y = group, x = value, fill = drugclass)) +
    geom_col(position = "dodge") +
    scale_fill_manual(values = c("Pooled" = "black", "SGLT2" = "#E69F00", "GLP1" = "#56B4E9", "SU" = "#CC79A7", "DPP4" = "#0072B2", "TZD" = "#D55E00", "MFN" = "grey"), breaks = c("GLP1", "DPP4", "SGLT2", "TZD", "SU"), labels = c("GLP-1RA", "DPP4i", "SGLT2i", "TZD", "SU"), name = " ", guide = guide_legend(reverse = TRUE)) +
    scale_y_discrete(limits = rev, labels = c("35+", "30-34.9", "25-29.9", "<25")) +
    ggtitle("BMI, kg/m2") +
    scale_x_continuous("Proportion of discontinuation", labels = scales::percent, limits = c(0, 0.3), breaks = seq(0, 0.3, 0.1)) +
    facet_wrap(~drugclass, nrow = 1, labeller = as_labeller(c("GLP1" = "GLP-1RA", "DPP4" = "DPP4i", "SGLT2" = "SGLT2i", "TZD" = "TZD", "SU" = "SU"))) +
    theme_bw() +
    guides(fill = guide_legend(nrow = 1, byrow = TRUE)) +
    theme(
      legend.position = "bottom",
      axis.text.x = element_blank(),
      axis.title.x = element_blank(),
      strip.text.x = element_blank()
    ),
  
  # HbA1c
  discontinuation_per_drug_not_combined %>%
    gather("key", "value", -drugclass, -analysis_var, -group, -variable) %>%
    mutate(
      key = factor(key, levels = c("perc_12m", "perc_6m", "perc_3m"))
    ) %>%
    filter(variable == "hba1c" & key == "perc_3m") %>%
    ggplot(aes(y = group, x = value, fill = drugclass)) +
    geom_col(position = "dodge") +
    scale_fill_manual(values = c("Pooled" = "black", "SGLT2" = "#E69F00", "GLP1" = "#56B4E9", "SU" = "#CC79A7", "DPP4" = "#0072B2", "TZD" = "#D55E00", "MFN" = "grey"), breaks = c("GLP1", "DPP4", "SGLT2", "TZD", "SU"), labels = c("GLP-1RA", "DPP4i", "SGLT2i", "TZD", "SU"), name = " ", guide = guide_legend(reverse = TRUE)) +
    scale_y_discrete(limits = rev, labels = c("86+", "75-86", "64-75", "53-64")) +
    ggtitle("HbA1c, mmol/mol") +
    scale_x_continuous("Proportion of discontinuation", labels = scales::percent, limits = c(0, 0.3), breaks = seq(0, 0.3, 0.1)) +
    facet_wrap(~drugclass, nrow = 1, labeller = as_labeller(c("GLP1" = "GLP-1RA", "DPP4" = "DPP4i", "SGLT2" = "SGLT2i", "TZD" = "TZD", "SU" = "SU"))) +
    theme_bw() +
    guides(fill = guide_legend(nrow = 1, byrow = TRUE)) +
    theme(
      legend.position = "bottom",
      axis.text.x = element_blank(),
      axis.title.x = element_blank(),
      strip.text.x = element_blank()
    ),
  
  # eGFR
  discontinuation_per_drug_not_combined %>%
    gather("key", "value", -drugclass, -analysis_var, -group, -variable) %>%
    mutate(
      key = factor(key, levels = c("perc_12m", "perc_6m", "perc_3m"))
    ) %>%
    filter(variable == "preegfr" & key == "perc_3m") %>%
    ggplot(aes(y = group, x = value, fill = drugclass)) +
    geom_col(position = "dodge") +
    scale_fill_manual(values = c("Pooled" = "black", "SGLT2" = "#E69F00", "GLP1" = "#56B4E9", "SU" = "#CC79A7", "DPP4" = "#0072B2", "TZD" = "#D55E00", "MFN" = "grey"), breaks = c("GLP1", "DPP4", "SGLT2", "TZD", "SU"), labels = c("GLP-1RA", "DPP4i", "SGLT2i", "TZD", "SU"), name = " ", guide = guide_legend(reverse = TRUE)) +
    scale_y_discrete(limits = rev, labels = c("90+", "60-90", "<60")) +
    ggtitle("eGFR, ml/min per 1.73 m2") +
    scale_x_continuous("Proportion of discontinuation", labels = scales::percent, limits = c(0, 0.3), breaks = seq(0, 0.3, 0.1)) +
    facet_wrap(~drugclass, nrow = 1, labeller = as_labeller(c("GLP1" = "GLP-1RA", "DPP4" = "DPP4i", "SGLT2" = "SGLT2i", "TZD" = "TZD", "SU" = "SU"))) +
    theme_bw() +
    guides(fill = guide_legend(nrow = 1, byrow = TRUE)) +
    theme(
      legend.position = "bottom",
      axis.text.x = element_blank(),
      axis.title.x = element_blank(),
      strip.text.x = element_blank()
    ),
  
  # Statins
  discontinuation_per_drug_not_combined %>%
    gather("key", "value", -drugclass, -analysis_var, -group, -variable) %>%
    mutate(
      key = factor(key, levels = c("perc_12m", "perc_6m", "perc_3m"))
    ) %>%
    filter(variable == "statins" & key == "perc_3m") %>%
    ggplot(aes(y = group, x = value, fill = drugclass)) +
    geom_col(position = "dodge") +
    scale_fill_manual(values = c("Pooled" = "black", "SGLT2" = "#E69F00", "GLP1" = "#56B4E9", "SU" = "#CC79A7", "DPP4" = "#0072B2", "TZD" = "#D55E00", "MFN" = "grey"), breaks = c("GLP1", "DPP4", "SGLT2", "TZD", "SU"), labels = c("GLP-1RA", "DPP4i", "SGLT2i", "TZD", "SU"), name = " ", guide = guide_legend(reverse = TRUE)) +
    scale_y_discrete(limits = rev, breaks = c(0, 1), labels = c("No", "Yes")) +
    ggtitle("Statins") +
    scale_x_continuous("Proportion of discontinuation", labels = scales::percent, limits = c(0, 0.3), breaks = seq(0, 0.3, 0.1)) +
    facet_wrap(~drugclass, nrow = 1, labeller = as_labeller(c("GLP1" = "GLP-1RA", "DPP4" = "DPP4i", "SGLT2" = "SGLT2i", "TZD" = "TZD", "SU" = "SU"))) +
    theme_bw() +
    guides(fill = guide_legend(nrow = 1, byrow = TRUE)) +
    theme(
      legend.position = "bottom",
      axis.text.x = element_blank(),
      axis.title.x = element_blank(),
      strip.text.x = element_blank()
    ),
  
  # Blood pressure medication
  discontinuation_per_drug_not_combined %>%
    gather("key", "value", -drugclass, -analysis_var, -group, -variable) %>%
    mutate(
      key = factor(key, levels = c("perc_12m", "perc_6m", "perc_3m"))
    ) %>%
    filter(variable == "statins" & key == "perc_3m") %>%
    ggplot(aes(y = group, x = value, fill = drugclass)) +
    geom_col(position = "dodge") +
    scale_fill_manual(values = c("Pooled" = "black", "SGLT2" = "#E69F00", "GLP1" = "#56B4E9", "SU" = "#CC79A7", "DPP4" = "#0072B2", "TZD" = "#D55E00", "MFN" = "grey"), breaks = c("GLP1", "DPP4", "SGLT2", "TZD", "SU"), labels = c("GLP-1RA", "DPP4i", "SGLT2i", "TZD", "SU"), name = " ", guide = guide_legend(reverse = TRUE)) +
    scale_y_discrete(limits = rev, breaks = c(0, 1), labels = c("No", "Yes")) +
    ggtitle("Blood pressure medication") +
    scale_x_continuous("Proportion of discontinuation", labels = scales::percent, limits = c(0, 0.3), breaks = seq(0, 0.3, 0.1)) +
    facet_wrap(~drugclass, nrow = 1, labeller = as_labeller(c("GLP1" = "GLP-1RA", "DPP4" = "DPP4i", "SGLT2" = "SGLT2i", "TZD" = "TZD", "SU" = "SU"))) +
    theme_bw() +
    guides(fill = guide_legend(nrow = 1, byrow = TRUE)) +
    theme(
      legend.position = "bottom",
      axis.text.x = element_blank(),
      axis.title.x = element_blank(),
      strip.text.x = element_blank()
    ),
  
  # Frailty
  discontinuation_per_drug_not_combined %>%
    gather("key", "value", -drugclass, -analysis_var, -group, -variable) %>%
    mutate(
      key = factor(key, levels = c("perc_12m", "perc_6m", "perc_3m"))
    ) %>%
    filter(variable == "frailty" & key == "perc_3m") %>%
    ggplot(aes(y = group, x = value, fill = drugclass)) +
    geom_col(position = "dodge") +
    scale_fill_manual(values = c("Pooled" = "black", "SGLT2" = "#E69F00", "GLP1" = "#56B4E9", "SU" = "#CC79A7", "DPP4" = "#0072B2", "TZD" = "#D55E00", "MFN" = "grey"), breaks = c("GLP1", "DPP4", "SGLT2", "TZD", "SU"), labels = c("GLP-1RA", "DPP4i", "SGLT2i", "TZD", "SU"), name = " ", guide = guide_legend(reverse = TRUE)) +
    scale_y_discrete(limits = rev, breaks = c(0, 1), labels = c("No", "Yes")) +
    ggtitle("Frailty") +
    scale_x_continuous("Proportion of discontinuation", labels = scales::percent, limits = c(0, 0.3), breaks = seq(0, 0.3, 0.1)) +
    facet_wrap(~drugclass, nrow = 1, labeller = as_labeller(c("GLP1" = "GLP-1RA", "DPP4" = "DPP4i", "SGLT2" = "SGLT2i", "TZD" = "TZD", "SU" = "SU"))) +
    theme_bw() +
    guides(fill = guide_legend(nrow = 1, byrow = TRUE)) +
    theme(
      legend.position = "bottom",
      axis.text.x = element_blank(),
      axis.title.x = element_blank(),
      strip.text.x = element_blank()
    ),
  
  # stopping metformin
  discontinuation_per_drug_not_combined %>%
    gather("key", "value", -drugclass, -analysis_var, -group, -variable) %>%
    mutate(
      key = factor(key, levels = c("perc_12m", "perc_6m", "perc_3m"))
    ) %>%
    filter(variable == "stopdrug_MFN" & key == "perc_3m") %>%
    ggplot(aes(y = group, x = value, fill = drugclass)) +
    geom_col() +
    scale_fill_manual(values = c("Pooled" = "black", "SGLT2" = "#E69F00", "GLP1" = "#56B4E9", "SU" = "#CC79A7", "DPP4" = "#0072B2", "TZD" = "#D55E00", "MFN" = "grey"), breaks = c("GLP1", "DPP4", "SGLT2", "TZD", "SU"), labels = c("GLP-1RA", "DPP4i", "SGLT2i", "TZD", "SU"), name = " ", guide = guide_legend(reverse = TRUE)) +
    scale_y_discrete(limits = rev, breaks = c(0, 1), labels = c("No", "Yes")) +
    ggtitle("History of Metformin discontinuation within 3-months") +
    scale_x_continuous("Proportion of discontinuation", labels = scales::percent, limits = c(0, 0.3), breaks = seq(0, 0.3, 0.1)) +
    facet_wrap(~drugclass, nrow = 1, labeller = as_labeller(c("GLP1" = "GLP-1RA", "DPP4" = "DPP4i", "SGLT2" = "SGLT2i", "TZD" = "TZD", "SU" = "SU"))) +
    theme_bw() +
    guides(fill = guide_legend(nrow = 1, byrow = TRUE)) +
    theme(
      legend.position = "bottom",
      axis.text.x = element_blank(),
      axis.title.x = element_blank(),
      strip.text.x = element_blank()
    ),
  
  # numdrugs
  discontinuation_per_drug_not_combined %>%
    gather("key", "value", -drugclass, -analysis_var, -group, -variable) %>%
    mutate(
      key = factor(key, levels = c("perc_12m", "perc_6m", "perc_3m"))
    ) %>%
    filter(variable == "numdrugs" & key == "perc_3m") %>%
    ggplot(aes(y = group, x = value, fill = drugclass)) +
    geom_col() +
    scale_fill_manual(values = c("Pooled" = "black", "SGLT2" = "#E69F00", "GLP1" = "#56B4E9", "SU" = "#CC79A7", "DPP4" = "#0072B2", "TZD" = "#D55E00", "MFN" = "grey"), breaks = c("GLP1", "DPP4", "SGLT2", "TZD", "SU"), labels = c("GLP-1RA", "DPP4i", "SGLT2i", "TZD", "SU"), name = " ", guide = guide_legend(reverse = TRUE)) +
    scale_y_discrete(limits = rev, labels = c("2+", "1", "0")) +
    ggtitle("Number of other current glucose-lowering therapies") +
    scale_x_continuous("Proportion of discontinuation", labels = scales::percent, limits = c(0, 0.3), breaks = seq(0, 0.3, 0.1)) +
    facet_wrap(~drugclass, nrow = 1, labeller = as_labeller(c("GLP1" = "GLP-1RA", "DPP4" = "DPP4i", "SGLT2" = "SGLT2i", "TZD" = "TZD", "SU" = "SU"))) +
    theme_bw() +
    guides(fill = guide_legend(nrow = 1, byrow = TRUE)) +
    theme(
      legend.position = "bottom",
      axis.text.x = element_blank(),
      axis.title.x = element_blank(),
      strip.text.x = element_blank()
    ),
  
  # drugline
  discontinuation_per_drug_not_combined %>%
    gather("key", "value", -drugclass, -analysis_var, -group, -variable) %>%
    mutate(
      key = factor(key, levels = c("perc_12m", "perc_6m", "perc_3m"))
    ) %>%
    filter(variable == "drugline" & key == "perc_3m") %>%
    ggplot(aes(y = group, x = value, fill = drugclass)) +
    geom_col(position = "dodge") +
    scale_fill_manual(values = c("Pooled" = "black", "SGLT2" = "#E69F00", "GLP1" = "#56B4E9", "SU" = "#CC79A7", "DPP4" = "#0072B2", "TZD" = "#D55E00", "MFN" = "grey"), breaks = c("GLP1", "DPP4", "SGLT2", "TZD", "SU"), labels = c("GLP-1RA", "DPP4i", "SGLT2i", "TZD", "SU"), name = " ", guide = guide_legend(reverse = TRUE)) +
    scale_y_discrete(limits = rev) +
    ggtitle("Number of therapy classes ever prescribed") +
    scale_x_continuous("Proportion of discontinuation", labels = scales::percent, limits = c(0, 0.3), breaks = seq(0, 0.3, 0.1)) +
    facet_wrap(~drugclass, nrow = 1, labeller = as_labeller(c("GLP1" = "GLP-1RA", "DPP4" = "DPP4i", "SGLT2" = "SGLT2i", "TZD" = "TZD", "SU" = "SU"))) +
    theme_bw() +
    guides(fill = guide_legend(nrow = 1, byrow = TRUE)) +
    theme(
      legend.position = "bottom",
      axis.text.x = element_text(size = 10),
      axis.title.x = element_text(size = 11),
      strip.text.x = element_blank()
    )
  
  , ncol = 1
  
) + 
  plot_layout(
    guides = "collect",
    height = c(4, 5, 2, 3, 4, 5, 4, 4, 3, 2, 2, 2, 2, 3, 3),
    # axis_titles = "collect" # not working, manually done in the code
  ) +
  plot_layout(
    axis_titles = "collect"
  ) &
  theme(
    legend.position = "none",
    plot.title.position = "plot",
    axis.title.y = element_blank(),
    axis.text.y = element_text(size = 10),
    panel.spacing = unit(0.4, "cm", data = NULL)
  )


dev.off()





# load dataset
cprd_dataset <- set_up_data(
  raw_data = "20240814_t2d_1stinstance",
  diagnosis = FALSE,
  therapies = c("DPP4", "GLP1", "SGLT2", "SU", "TZD"),
  dataset = "6m.disc.dataset",
  follow_up = "6-months",
  full_prescribing_history = TRUE
) %>%
  drop_na(-prehdl, -stopdrug_12m_6mFU)



# Start plot by variables


pdf("results/figures/03.6m_discontinuation_by_variable.pdf", width = 7, height = 15)

patchwork::wrap_plots(
  
  # Overall
  discontinuation_per_drug_combined %>%
    gather("key", "value", -drugclass, -analysis_var, -group, -variable) %>%
    mutate(
      key = factor(key, levels = c("perc_12m", "perc_6m", "perc_3m"))
    ) %>%
    filter(variable == "overall" & key == "perc_6m") %>%
    ggplot(aes(y = analysis_var, x = value, fill = drugclass)) +
    geom_col() +
    scale_fill_manual(values = c("other" = "black", "SGLT2" = "#E69F00", "GLP1" = "#56B4E9", "SU" = "#CC79A7", "DPP4" = "#0072B2", "TZD" = "#D55E00", "MFN" = "grey"), breaks = c("other"), labels = c("Other glucose-lowering therapies"), name = " ", guide = guide_legend(reverse = FALSE)) +
    scale_y_discrete(breaks = c(1), labels = c("Overall")) +
    scale_x_continuous("Proportion of discontinuation", labels = scales::percent, limits = c(0, 0.3)) +
    facet_wrap(~drugclass, nrow = 1, labeller = as_labeller(c("other" = "Other glucose-lowering therapies"))) +
    theme_bw() +
    guides(fill = guide_legend(nrow = 1, byrow = TRUE)) +
    theme(
      legend.position = "bottom",
      axis.text.x = element_blank(),
      axis.title.x = element_blank()
    ),
  
  # Age
  discontinuation_per_drug_combined %>%
    gather("key", "value", -drugclass, -analysis_var, -group, -variable) %>%
    mutate(
      key = factor(key, levels = c("perc_12m", "perc_6m", "perc_3m"))
    ) %>%
    filter(variable == "dstartdate_age" & key == "perc_6m") %>%
    ggplot(aes(y = group, x = value, fill = drugclass)) +
    geom_col() +
    scale_fill_manual(values = c("other" = "black", "SGLT2" = "#E69F00", "GLP1" = "#56B4E9", "SU" = "#CC79A7", "DPP4" = "#0072B2", "TZD" = "#D55E00", "MFN" = "grey"), breaks = c("other"), labels = c("Other glucose-lowering therapies"), name = " ", guide = guide_legend(reverse = FALSE)) +
    scale_y_discrete(limits = rev, labels = c("70+", "60-69", "50-59", "<50")) +
    scale_x_continuous("Proportion of discontinuation", labels = scales::percent, limits = c(0, 0.3)) +
    ggtitle("Age group, years") +
    facet_wrap(~drugclass, nrow = 1, labeller = as_labeller(c("other" = "Other glucose-lowering therapies"))) +
    theme_bw() +
    guides(fill = guide_legend(nrow = 1, byrow = TRUE)) +
    theme(
      plot.title.position = "plot",
      legend.position = "bottom",
      axis.text.x = element_blank(),
      axis.title.x = element_blank(),
      strip.text.x = element_blank()
    ),
  
  # Diabetes duration
  discontinuation_per_drug_combined %>%
    gather("key", "value", -drugclass, -analysis_var, -group, -variable) %>%
    mutate(
      key = factor(key, levels = c("perc_12m", "perc_6m", "perc_3m"))
    ) %>%
    filter(variable == "dstartdate_dm_dur" & key == "perc_6m") %>%
    ggplot(aes(y = group, x = value, fill = drugclass)) +
    geom_col() +
    scale_fill_manual(values = c("other" = "black", "SGLT2" = "#E69F00", "GLP1" = "#56B4E9", "SU" = "#CC79A7", "DPP4" = "#0072B2", "TZD" = "#D55E00", "MFN" = "grey"), breaks = c("other"), labels = c("Other glucose-lowering therapies"), name = " ", guide = guide_legend(reverse = FALSE)) +
    scale_y_discrete(limits = rev, labels = c("10+", "6-9", "3-5", "1-2", "<1")) +
    ggtitle("Duration of diabetes, years") +
    scale_x_continuous("Proportion of discontinuation", labels = scales::percent, limits = c(0, 0.3)) +
    facet_wrap(~drugclass, nrow = 1, labeller = as_labeller(c("other" = "Other glucose-lowering therapies"))) +
    theme_bw() +
    guides(fill = guide_legend(nrow = 1, byrow = TRUE)) +
    theme(
      plot.title.position = "plot",
      legend.position = "bottom",
      axis.text.x = element_blank(),
      axis.title.x = element_blank(),
      strip.text.x = element_blank()
    ),
  
  # Sex
  discontinuation_per_drug_combined %>%
    gather("key", "value", -drugclass, -analysis_var, -group, -variable) %>%
    mutate(
      key = factor(key, levels = c("perc_12m", "perc_6m", "perc_3m"))
    ) %>%
    filter(variable == "gender" & key == "perc_6m") %>%
    ggplot(aes(y = group, x = value, fill = drugclass)) +
    geom_col() +
    scale_fill_manual(values = c("other" = "black", "SGLT2" = "#E69F00", "GLP1" = "#56B4E9", "SU" = "#CC79A7", "DPP4" = "#0072B2", "TZD" = "#D55E00", "MFN" = "grey"), breaks = c("other"), labels = c("Other glucose-lowering therapies"), name = " ", guide = guide_legend(reverse = FALSE)) +
    scale_y_discrete(limits = rev, labels = c("Female", "Male")) +
    scale_x_continuous("Proportion of discontinuation", labels = scales::percent, limits = c(0, 0.3)) +
    ggtitle("Sex") +
    facet_wrap(~drugclass, nrow = 1, labeller = as_labeller(c("other" = "Other glucose-lowering therapies"))) +
    theme_bw() +
    guides(fill = guide_legend(nrow = 1, byrow = TRUE)) +
    theme(
      plot.title.position = "plot",
      legend.position = "bottom",
      axis.text.x = element_blank(),
      axis.title.x = element_blank(),
      strip.text.x = element_blank()
    ),
  
  # Smoking
  discontinuation_per_drug_combined %>%
    gather("key", "value", -drugclass, -analysis_var, -group, -variable) %>%
    mutate(
      key = factor(key, levels = c("perc_12m", "perc_6m", "perc_3m"))
    ) %>%
    filter(variable == "smoking" & key == "perc_6m") %>%
    ggplot(aes(y = group, x = value, fill = drugclass)) +
    geom_col() +
    scale_fill_manual(values = c("other" = "black", "SGLT2" = "#E69F00", "GLP1" = "#56B4E9", "SU" = "#CC79A7", "DPP4" = "#0072B2", "TZD" = "#D55E00", "MFN" = "grey"), breaks = c("other"), labels = c("Other glucose-lowering therapies"), name = " ", guide = guide_legend(reverse = FALSE)) +
    ggtitle("Smoking") +
    scale_x_continuous("Proportion of discontinuation", labels = scales::percent, limits = c(0, 0.3)) +
    facet_wrap(~drugclass, nrow = 1, labeller = as_labeller(c("other" = "Other glucose-lowering therapies"))) +
    theme_bw() +
    guides(fill = guide_legend(nrow = 1, byrow = TRUE)) +
    theme(
      plot.title.position = "plot",
      legend.position = "bottom",
      axis.text.x = element_blank(),
      axis.title.x = element_blank(),
      strip.text.x = element_blank()
    ),
  
  # Ethnicity
  discontinuation_per_drug_combined %>%
    gather("key", "value", -drugclass, -analysis_var, -group, -variable) %>%
    mutate(
      key = factor(key, levels = c("perc_12m", "perc_6m", "perc_3m"))
    ) %>%
    filter(variable == "ethnicity" & key == "perc_6m") %>%
    ggplot(aes(y = group, x = value, fill = drugclass)) +
    geom_col() +
    scale_fill_manual(values = c("other" = "black", "SGLT2" = "#E69F00", "GLP1" = "#56B4E9", "SU" = "#CC79A7", "DPP4" = "#0072B2", "TZD" = "#D55E00", "MFN" = "grey"), breaks = c("other"), labels = c("Other glucose-lowering therapies"), name = " ", guide = guide_legend(reverse = FALSE)) +
    ggtitle("Ethnicity") +
    scale_x_continuous("Proportion of discontinuation", labels = scales::percent, limits = c(0, 0.3)) +
    facet_wrap(~drugclass, nrow = 1, labeller = as_labeller(c("other" = "Other glucose-lowering therapies"))) +
    theme_bw() +
    guides(fill = guide_legend(nrow = 1, byrow = TRUE)) +
    theme(
      plot.title.position = "plot",
      legend.position = "bottom",
      axis.text.x = element_blank(),
      axis.title.x = element_blank(),
      strip.text.x = element_blank()
    ),
  
  # Quintile
  discontinuation_per_drug_combined %>%
    gather("key", "value", -drugclass, -analysis_var, -group, -variable) %>%
    mutate(
      key = factor(key, levels = c("perc_12m", "perc_6m", "perc_3m"))
    ) %>%
    filter(variable == "imd2015" & key == "perc_6m") %>%
    ggplot(aes(y = group, x = value, fill = drugclass)) +
    geom_col() +
    scale_fill_manual(values = c("other" = "black", "SGLT2" = "#E69F00", "GLP1" = "#56B4E9", "SU" = "#CC79A7", "DPP4" = "#0072B2", "TZD" = "#D55E00", "MFN" = "grey"), breaks = c("other"), labels = c("Other glucose-lowering therapies"), name = " ", guide = guide_legend(reverse = FALSE)) +
    scale_y_discrete(limits = rev, labels = c("5 (most deprived)", "4", "3", "2", "1 (least deprived)")) +
    scale_x_continuous("Proportion of discontinuation", labels = scales::percent, limits = c(0, 0.3)) +
    ggtitle("Index of multiple deprivation") +
    facet_wrap(~drugclass, nrow = 1, labeller = as_labeller(c("other" = "Other glucose-lowering therapies"))) +
    theme_bw() +
    guides(fill = guide_legend(nrow = 1, byrow = TRUE)) +
    theme(
      plot.title.position = "plot",
      legend.position = "bottom",
      axis.text.x = element_blank(),
      axis.title.x = element_blank(),
      strip.text.x = element_blank()
    ),
  
  # BMI
  discontinuation_per_drug_combined %>%
    gather("key", "value", -drugclass, -analysis_var, -group, -variable) %>%
    mutate(
      key = factor(key, levels = c("perc_12m", "perc_6m", "perc_3m"))
    ) %>%
    filter(variable == "bmi" & key == "perc_6m") %>%
    ggplot(aes(y = group, x = value, fill = drugclass)) +
    geom_col() +
    scale_fill_manual(values = c("other" = "black", "SGLT2" = "#E69F00", "GLP1" = "#56B4E9", "SU" = "#CC79A7", "DPP4" = "#0072B2", "TZD" = "#D55E00", "MFN" = "grey"), breaks = c("other"), labels = c("Other glucose-lowering therapies"), name = " ", guide = guide_legend(reverse = FALSE)) +
    scale_y_discrete(limits = rev, labels = c("35+", "30-34.9", "25-29.9", "<25")) +
    ggtitle("BMI, kg/m2") +
    scale_x_continuous("Proportion of discontinuation", labels = scales::percent, limits = c(0, 0.3)) +
    facet_wrap(~drugclass, nrow = 1, labeller = as_labeller(c("other" = "Other glucose-lowering therapies"))) +
    theme_bw() +
    guides(fill = guide_legend(nrow = 1, byrow = TRUE)) +
    theme(
      plot.title.position = "plot",
      legend.position = "bottom",
      axis.text.x = element_blank(),
      axis.title.x = element_blank(),
      strip.text.x = element_blank()
    ),
  
  # HbA1c
  discontinuation_per_drug_combined %>%
    gather("key", "value", -drugclass, -analysis_var, -group, -variable) %>%
    mutate(
      key = factor(key, levels = c("perc_12m", "perc_6m", "perc_3m"))
    ) %>%
    filter(variable == "hba1c" & key == "perc_6m") %>%
    ggplot(aes(y = group, x = value, fill = drugclass)) +
    geom_col() +
    scale_fill_manual(values = c("other" = "black", "SGLT2" = "#E69F00", "GLP1" = "#56B4E9", "SU" = "#CC79A7", "DPP4" = "#0072B2", "TZD" = "#D55E00", "MFN" = "grey"), breaks = c("other"), labels = c("Other glucose-lowering therapies"), name = " ", guide = guide_legend(reverse = FALSE)) +
    scale_y_discrete(limits = rev, labels = c("86+", "75-86", "64-75", "53-64")) +
    ggtitle("HbA1c, mmol/mol") +
    scale_x_continuous("Proportion of discontinuation", labels = scales::percent, limits = c(0, 0.3)) +
    facet_wrap(~drugclass, nrow = 1, labeller = as_labeller(c("other" = "Other glucose-lowering therapies"))) +
    theme_bw() +
    guides(fill = guide_legend(nrow = 1, byrow = TRUE)) +
    theme(
      plot.title.position = "plot",
      legend.position = "bottom",
      axis.text.x = element_blank(),
      axis.title.x = element_blank(),
      strip.text.x = element_blank()
    ),
  
  # eGFR
  discontinuation_per_drug_combined %>%
    gather("key", "value", -drugclass, -analysis_var, -group, -variable) %>%
    mutate(
      key = factor(key, levels = c("perc_12m", "perc_6m", "perc_3m"))
    ) %>%
    filter(variable == "preegfr" & key == "perc_6m") %>%
    ggplot(aes(y = group, x = value, fill = drugclass)) +
    geom_col() +
    scale_fill_manual(values = c("other" = "black", "SGLT2" = "#E69F00", "GLP1" = "#56B4E9", "SU" = "#CC79A7", "DPP4" = "#0072B2", "TZD" = "#D55E00", "MFN" = "grey"), breaks = c("other"), labels = c("Other glucose-lowering therapies"), name = " ", guide = guide_legend(reverse = FALSE)) +
    scale_y_discrete(limits = rev, labels = c("90+", "60-90", "<60")) +
    ggtitle("eGFR, ml/min per 1.73 m2") +
    scale_x_continuous("Proportion of discontinuation", labels = scales::percent, limits = c(0, 0.3)) +
    facet_wrap(~drugclass, nrow = 1, labeller = as_labeller(c("other" = "Other glucose-lowering therapies"))) +
    theme_bw() +
    guides(fill = guide_legend(nrow = 1, byrow = TRUE)) +
    theme(
      plot.title.position = "plot",
      legend.position = "bottom",
      axis.text.x = element_blank(),
      axis.title.x = element_blank(),
      strip.text.x = element_blank()
    ),
  
  # Statins
  discontinuation_per_drug_combined %>%
    gather("key", "value", -drugclass, -analysis_var, -group, -variable) %>%
    mutate(
      key = factor(key, levels = c("perc_12m", "perc_6m", "perc_3m"))
    ) %>%
    filter(variable == "statins" & key == "perc_6m") %>%
    ggplot(aes(y = group, x = value, fill = drugclass)) +
    geom_col() +
    scale_fill_manual(values = c("other" = "black", "SGLT2" = "#E69F00", "GLP1" = "#56B4E9", "SU" = "#CC79A7", "DPP4" = "#0072B2", "TZD" = "#D55E00", "MFN" = "grey"), breaks = c("other"), labels = c("Other glucose-lowering therapies"), name = " ", guide = guide_legend(reverse = FALSE)) +
    scale_y_discrete(limits = rev, breaks = c(0, 1), labels = c("No", "Yes")) +
    ggtitle("Statins") +
    scale_x_continuous("Proportion of discontinuation", labels = scales::percent, limits = c(0, 0.3)) +
    facet_wrap(~drugclass, nrow = 1, labeller = as_labeller(c("other" = "Other glucose-lowering therapies"))) +
    theme_bw() +
    guides(fill = guide_legend(nrow = 1, byrow = TRUE)) +
    theme(
      plot.title.position = "plot",
      legend.position = "bottom",
      axis.text.x = element_blank(),
      axis.title.x = element_blank(),
      strip.text.x = element_blank()
    ),
  
  # Bloodmed
  discontinuation_per_drug_combined %>%
    gather("key", "value", -drugclass, -analysis_var, -group, -variable) %>%
    mutate(
      key = factor(key, levels = c("perc_12m", "perc_6m", "perc_3m"))
    ) %>%
    filter(variable == "bloodmed" & key == "perc_6m") %>%
    ggplot(aes(y = group, x = value, fill = drugclass)) +
    geom_col() +
    scale_fill_manual(values = c("other" = "black", "SGLT2" = "#E69F00", "GLP1" = "#56B4E9", "SU" = "#CC79A7", "DPP4" = "#0072B2", "TZD" = "#D55E00", "MFN" = "grey"), breaks = c("other"), labels = c("Other glucose-lowering therapies"), name = " ", guide = guide_legend(reverse = FALSE)) +
    scale_y_discrete(limits = rev, breaks = c(0, 1), labels = c("No", "Yes")) +
    ggtitle("Blood pressure medication") +
    scale_x_continuous("Proportion of discontinuation", labels = scales::percent, limits = c(0, 0.3)) +
    facet_wrap(~drugclass, nrow = 1, labeller = as_labeller(c("other" = "Other glucose-lowering therapies"))) +
    theme_bw() +
    guides(fill = guide_legend(nrow = 1, byrow = TRUE)) +
    theme(
      plot.title.position = "plot",
      legend.position = "bottom",
      axis.text.x = element_blank(),
      axis.title.x = element_blank(),
      strip.text.x = element_blank()
    ),
  
  # Frailty
  discontinuation_per_drug_combined %>%
    gather("key", "value", -drugclass, -analysis_var, -group, -variable) %>%
    mutate(
      key = factor(key, levels = c("perc_12m", "perc_6m", "perc_3m"))
    ) %>%
    filter(variable == "frailty" & key == "perc_6m") %>%
    ggplot(aes(y = group, x = value, fill = drugclass)) +
    geom_col() +
    scale_fill_manual(values = c("other" = "black", "SGLT2" = "#E69F00", "GLP1" = "#56B4E9", "SU" = "#CC79A7", "DPP4" = "#0072B2", "TZD" = "#D55E00", "MFN" = "grey"), breaks = c("other"), labels = c("Other glucose-lowering therapies"), name = " ", guide = guide_legend(reverse = FALSE)) +
    scale_y_discrete(limits = rev, breaks = c(0, 1), labels = c("No", "Yes")) +
    ggtitle("Frailty") +
    scale_x_continuous("Proportion of discontinuation", labels = scales::percent, limits = c(0, 0.3)) +
    facet_wrap(~drugclass, nrow = 1, labeller = as_labeller(c("other" = "Other glucose-lowering therapies"))) +
    theme_bw() +
    guides(fill = guide_legend(nrow = 1, byrow = TRUE)) +
    theme(
      plot.title.position = "plot",
      legend.position = "bottom",
      axis.text.x = element_blank(),
      axis.title.x = element_blank(),
      strip.text.x = element_blank()
    ),
  
  # stopping metformin
  discontinuation_per_drug_combined %>%
    gather("key", "value", -drugclass, -analysis_var, -group, -variable) %>%
    mutate(
      key = factor(key, levels = c("perc_12m", "perc_6m", "perc_3m"))
    ) %>%
    filter(variable == "stopdrug_MFN" & key == "perc_6m") %>%
    ggplot(aes(y = group, x = value, fill = drugclass)) +
    geom_col() +
    scale_fill_manual(values = c("other" = "black", "SGLT2" = "#E69F00", "GLP1" = "#56B4E9", "SU" = "#CC79A7", "DPP4" = "#0072B2", "TZD" = "#D55E00", "MFN" = "grey"), breaks = c("other"), labels = c("Other glucose-lowering therapies"), name = " ", guide = guide_legend(reverse = FALSE)) +
    scale_y_discrete(limits = rev, breaks = c(0, 1), labels = c("No", "Yes")) +
    ggtitle("History of Metformin discontinuation within 3-months") +
    scale_x_continuous("Proportion of discontinuation", labels = scales::percent, limits = c(0, 0.3)) +
    facet_wrap(~drugclass, nrow = 1, labeller = as_labeller(c("other" = "Other glucose-lowering therapies"))) +
    theme_bw() +
    guides(fill = guide_legend(nrow = 1, byrow = TRUE)) +
    theme(
      plot.title.position = "plot",
      legend.position = "bottom",
      axis.text.x = element_blank(),
      axis.title.x = element_blank(),
      strip.text.x = element_blank()
    ),
  
  # Numdrugs
  discontinuation_per_drug_combined %>%
    gather("key", "value", -drugclass, -analysis_var, -group, -variable) %>%
    mutate(
      key = factor(key, levels = c("perc_12m", "perc_6m", "perc_3m"))
    ) %>%
    filter(variable == "numdrugs" & key == "perc_6m") %>%
    ggplot(aes(y = group, x = value, fill = drugclass)) +
    geom_col() +
    scale_fill_manual(values = c("other" = "black", "SGLT2" = "#E69F00", "GLP1" = "#56B4E9", "SU" = "#CC79A7", "DPP4" = "#0072B2", "TZD" = "#D55E00", "MFN" = "grey"), breaks = c("other"), labels = c("Other glucose-lowering therapies"), name = " ", guide = guide_legend(reverse = FALSE)) +
    scale_y_discrete(limits = rev, labels = c("2+", "1", "0")) +
    ggtitle("Number of other current glucose-lowering therapies") +
    scale_x_continuous("Proportion of discontinuation", labels = scales::percent, limits = c(0, 0.3)) +
    facet_wrap(~drugclass, nrow = 1, labeller = as_labeller(c("other" = "Other glucose-lowering therapies"))) +
    theme_bw() +
    guides(fill = guide_legend(nrow = 1, byrow = TRUE)) +
    theme(
      plot.title.position = "plot",
      legend.position = "bottom",
      axis.text.x = element_blank(),
      axis.title.x = element_blank(),
      strip.text.x = element_blank()
    ),
  
  # drugline
  discontinuation_per_drug_combined %>%
    gather("key", "value", -drugclass, -analysis_var, -group, -variable) %>%
    mutate(
      key = factor(key, levels = c("perc_12m", "perc_6m", "perc_3m"))
    ) %>%
    filter(variable == "drugline" & key == "perc_6m") %>%
    ggplot(aes(y = group, x = value, fill = drugclass)) +
    geom_col() +
    scale_fill_manual(values = c("other" = "black", "SGLT2" = "#E69F00", "GLP1" = "#56B4E9", "SU" = "#CC79A7", "DPP4" = "#0072B2", "TZD" = "#D55E00", "MFN" = "grey"), breaks = c("other"), labels = c("Other glucose-lowering therapies"), name = " ", guide = guide_legend(reverse = FALSE)) +
    scale_y_discrete(limits = rev) +
    ggtitle("Number of therapy classes ever prescribed") +
    scale_x_continuous("Proportion of discontinuation", labels = scales::percent, limits = c(0, 0.3)) +
    facet_wrap(~drugclass, nrow = 1, labeller = as_labeller(c("other" = "Other glucose-lowering therapies"))) +
    theme_bw() +
    guides(fill = guide_legend(nrow = 1, byrow = TRUE)) +
    theme(
      plot.title.position = "plot",
      legend.position = "bottom",
      axis.text.x = element_text(size = 10),
      axis.title.x = element_text(size = 11),
      strip.text.x = element_blank()
    )
  
  , ncol = 1
  
) + 
  plot_layout(
    guides = "collect",
    height = c(1, 4, 5, 2, 3, 4, 5, 4, 4, 3, 2, 2, 2, 2, 3, 3),
    # axis_titles = "collect" # not working, manually done in the code
  ) +
  plot_layout(
    axis_titles = "collect"
  ) &
  theme(
    legend.position = "none",
    plot.title.position = "plot",
    axis.title.y = element_blank(),
    axis.text.y = element_text(size = 10),
    panel.spacing = unit(0.3, "cm", data = NULL)
  )


dev.off()



pdf("results/figures/03.6m_discontinuation_by_variable_by_drug.pdf", width = 14, height = 15)

patchwork::wrap_plots(
  
  # Age - Metformin
  discontinuation_per_drug_not_combined %>%
    gather("key", "value", -drugclass, -analysis_var, -group, -variable) %>%
    mutate(
      key = factor(key, levels = c("perc_12m", "perc_6m", "perc_3m"))
    ) %>%
    filter(variable == "dstartdate_age" & key == "perc_6m") %>%
    ggplot(aes(y = group, x = value, fill = drugclass)) +
    geom_col(position = "dodge") +
    scale_fill_manual(values = c("Pooled" = "black", "SGLT2" = "#E69F00", "GLP1" = "#56B4E9", "SU" = "#CC79A7", "DPP4" = "#0072B2", "TZD" = "#D55E00", "MFN" = "grey"), breaks = c("GLP1", "DPP4", "SGLT2", "TZD", "SU"), labels = c("GLP-1RA", "DPP4i", "SGLT2i", "TZD", "SU"), name = " ", guide = guide_legend(reverse = TRUE)) +
    scale_y_discrete(limits = rev, labels = c("70+", "60-69", "50-59", "<50")) +
    ggtitle("Age group, years") +
    scale_x_continuous("Proportion of discontinuation", labels = scales::percent, limits = c(0, 0.5), breaks = seq(0, 0.5, 0.1)) +
    facet_wrap(~drugclass, nrow = 1, labeller = as_labeller(c("GLP1" = "GLP-1RA", "DPP4" = "DPP4i", "SGLT2" = "SGLT2i", "TZD" = "TZD", "SU" = "SU"))) +
    theme_bw() +
    guides(fill = guide_legend(nrow = 1, byrow = TRUE)) +
    theme(
      legend.position = "bottom",
      axis.text.x = element_blank(),
      axis.title.x = element_blank()
    ),
  
  # Diabetes duration
  discontinuation_per_drug_not_combined %>%
    gather("key", "value", -drugclass, -analysis_var, -group, -variable) %>%
    mutate(
      key = factor(key, levels = c("perc_12m", "perc_6m", "perc_3m"))
    ) %>%
    filter(variable == "dstartdate_dm_dur" & key == "perc_6m") %>%
    ggplot(aes(y = group, x = value, fill = drugclass)) +
    geom_col(position = "dodge") +
    scale_fill_manual(values = c("Pooled" = "black", "SGLT2" = "#E69F00", "GLP1" = "#56B4E9", "SU" = "#CC79A7", "DPP4" = "#0072B2", "TZD" = "#D55E00", "MFN" = "grey"), breaks = c("GLP1", "DPP4", "SGLT2", "TZD", "SU"), labels = c("GLP-1RA", "DPP4i", "SGLT2i", "TZD", "SU"), name = " ", guide = guide_legend(reverse = TRUE)) +
    scale_y_discrete(limits = rev, labels = c("10+", "6-9", "3-5", "1-2", "<1")) +
    ggtitle("Duration of diabetes, years") +
    scale_x_continuous("Proportion of discontinuation", labels = scales::percent, limits = c(0, 0.5), breaks = seq(0, 0.5, 0.1)) +
    facet_wrap(~drugclass, nrow = 1, labeller = as_labeller(c("GLP1" = "GLP-1RA", "DPP4" = "DPP4i", "SGLT2" = "SGLT2i", "TZD" = "TZD", "SU" = "SU"))) +
    theme_bw() +
    guides(fill = guide_legend(nrow = 1, byrow = TRUE)) +
    theme(
      legend.position = "bottom",
      axis.text.x = element_blank(),
      axis.title.x = element_blank(),
      strip.text.x = element_blank()
    ),
  
  # Sex 
  discontinuation_per_drug_not_combined %>%
    gather("key", "value", -drugclass, -analysis_var, -group, -variable) %>%
    mutate(
      key = factor(key, levels = c("perc_12m", "perc_6m", "perc_3m"))
    ) %>%
    filter(variable == "gender" & key == "perc_6m") %>%
    ggplot(aes(y = group, x = value, fill = drugclass)) +
    geom_col(position = "dodge") +
    scale_fill_manual(values = c("Pooled" = "black", "SGLT2" = "#E69F00", "GLP1" = "#56B4E9", "SU" = "#CC79A7", "DPP4" = "#0072B2", "TZD" = "#D55E00", "MFN" = "grey"), breaks = c("GLP1", "DPP4", "SGLT2", "TZD", "SU"), labels = c("GLP-1RA", "DPP4i", "SGLT2i", "TZD", "SU"), name = " ", guide = guide_legend(reverse = TRUE)) +
    scale_y_discrete(limits = rev, labels = c("Female", "Male")) +
    ggtitle("Sex") +
    scale_x_continuous("Proportion of discontinuation", labels = scales::percent, limits = c(0, 0.5), breaks = seq(0, 0.5, 0.1)) +
    facet_wrap(~drugclass, nrow = 1, labeller = as_labeller(c("GLP1" = "GLP-1RA", "DPP4" = "DPP4i", "SGLT2" = "SGLT2i", "TZD" = "TZD", "SU" = "SU"))) +
    theme_bw() +
    guides(fill = guide_legend(nrow = 1, byrow = TRUE)) +
    theme(
      legend.position = "bottom",
      axis.text.x = element_blank(),
      axis.title.x = element_blank(),
      strip.text.x = element_blank()
    ),
  
  
  # Smoking 
  discontinuation_per_drug_not_combined %>%
    gather("key", "value", -drugclass, -analysis_var, -group, -variable) %>%
    mutate(
      key = factor(key, levels = c("perc_12m", "perc_6m", "perc_3m"))
    ) %>%
    filter(variable == "smoking" & key == "perc_6m") %>%
    ggplot(aes(y = group, x = value, fill = drugclass)) +
    geom_col(position = "dodge") +
    scale_fill_manual(values = c("Pooled" = "black", "SGLT2" = "#E69F00", "GLP1" = "#56B4E9", "SU" = "#CC79A7", "DPP4" = "#0072B2", "TZD" = "#D55E00", "MFN" = "grey"), breaks = c("GLP1", "DPP4", "SGLT2", "TZD", "SU"), labels = c("GLP-1RA", "DPP4i", "SGLT2i", "TZD", "SU"), name = " ", guide = guide_legend(reverse = TRUE)) +
    ggtitle("Smoking") +
    scale_x_continuous("Proportion of discontinuation", labels = scales::percent, limits = c(0, 0.5), breaks = seq(0, 0.5, 0.1)) +
    facet_wrap(~drugclass, nrow = 1, labeller = as_labeller(c("GLP1" = "GLP-1RA", "DPP4" = "DPP4i", "SGLT2" = "SGLT2i", "TZD" = "TZD", "SU" = "SU"))) +
    theme_bw() +
    guides(fill = guide_legend(nrow = 1, byrow = TRUE)) +
    theme(
      legend.position = "bottom",
      axis.text.x = element_blank(),
      axis.title.x = element_blank(),
      strip.text.x = element_blank()
    ),
  
  # Ethnicity
  discontinuation_per_drug_not_combined %>%
    gather("key", "value", -drugclass, -analysis_var, -group, -variable) %>%
    mutate(
      key = factor(key, levels = c("perc_12m", "perc_6m", "perc_3m"))
    ) %>%
    filter(variable == "ethnicity" & key == "perc_6m") %>%
    ggplot(aes(y = group, x = value, fill = drugclass)) +
    geom_col(position = "dodge") +
    scale_fill_manual(values = c("Pooled" = "black", "SGLT2" = "#E69F00", "GLP1" = "#56B4E9", "SU" = "#CC79A7", "DPP4" = "#0072B2", "TZD" = "#D55E00", "MFN" = "grey"), breaks = c("GLP1", "DPP4", "SGLT2", "TZD", "SU"), labels = c("GLP-1RA", "DPP4i", "SGLT2i", "TZD", "SU"), name = " ", guide = guide_legend(reverse = TRUE)) +
    ggtitle("Ethnicity") +
    scale_x_continuous("Proportion of discontinuation", labels = scales::percent, limits = c(0, 0.5), breaks = seq(0, 0.5, 0.1)) +
    facet_wrap(~drugclass, nrow = 1, labeller = as_labeller(c("GLP1" = "GLP-1RA", "DPP4" = "DPP4i", "SGLT2" = "SGLT2i", "TZD" = "TZD", "SU" = "SU"))) +
    theme_bw() +
    guides(fill = guide_legend(nrow = 1, byrow = TRUE)) +
    theme(
      legend.position = "bottom",
      axis.text.x = element_blank(),
      axis.title.x = element_blank(),
      strip.text.x = element_blank()
    ),
  
  # Quintile
  discontinuation_per_drug_not_combined %>%
    gather("key", "value", -drugclass, -analysis_var, -group, -variable) %>%
    mutate(
      key = factor(key, levels = c("perc_12m", "perc_6m", "perc_3m"))
    ) %>%
    filter(variable == "imd2015" & key == "perc_6m") %>%
    ggplot(aes(y = group, x = value, fill = drugclass)) +
    geom_col(position = "dodge") +
    scale_fill_manual(values = c("Pooled" = "black", "SGLT2" = "#E69F00", "GLP1" = "#56B4E9", "SU" = "#CC79A7", "DPP4" = "#0072B2", "TZD" = "#D55E00", "MFN" = "grey"), breaks = c("GLP1", "DPP4", "SGLT2", "TZD", "SU"), labels = c("GLP-1RA", "DPP4i", "SGLT2i", "TZD", "SU"), name = " ", guide = guide_legend(reverse = TRUE)) +
    scale_y_discrete(limits = rev, labels = c("5 (most deprived)", "4", "3", "2", "1 (least deprived)")) +
    ggtitle("Index of multiple deprivation") +
    scale_x_continuous("Proportion of discontinuation", labels = scales::percent, limits = c(0, 0.5), breaks = seq(0, 0.5, 0.1)) +
    facet_wrap(~drugclass, nrow = 1, labeller = as_labeller(c("GLP1" = "GLP-1RA", "DPP4" = "DPP4i", "SGLT2" = "SGLT2i", "TZD" = "TZD", "SU" = "SU"))) +
    theme_bw() +
    guides(fill = guide_legend(nrow = 1, byrow = TRUE)) +
    theme(
      legend.position = "bottom",
      axis.text.x = element_blank(),
      axis.title.x = element_blank(),
      strip.text.x = element_blank()
    ),
  
  # BMI
  discontinuation_per_drug_not_combined %>%
    gather("key", "value", -drugclass, -analysis_var, -group, -variable) %>%
    mutate(
      key = factor(key, levels = c("perc_12m", "perc_6m", "perc_3m"))
    ) %>%
    filter(variable == "bmi" & key == "perc_6m") %>%
    ggplot(aes(y = group, x = value, fill = drugclass)) +
    geom_col(position = "dodge") +
    scale_fill_manual(values = c("Pooled" = "black", "SGLT2" = "#E69F00", "GLP1" = "#56B4E9", "SU" = "#CC79A7", "DPP4" = "#0072B2", "TZD" = "#D55E00", "MFN" = "grey"), breaks = c("GLP1", "DPP4", "SGLT2", "TZD", "SU"), labels = c("GLP-1RA", "DPP4i", "SGLT2i", "TZD", "SU"), name = " ", guide = guide_legend(reverse = TRUE)) +
    scale_y_discrete(limits = rev, labels = c("35+", "30-34.9", "25-29.9", "<25")) +
    ggtitle("BMI, kg/m2") +
    scale_x_continuous("Proportion of discontinuation", labels = scales::percent, limits = c(0, 0.5), breaks = seq(0, 0.5, 0.1)) +
    facet_wrap(~drugclass, nrow = 1, labeller = as_labeller(c("GLP1" = "GLP-1RA", "DPP4" = "DPP4i", "SGLT2" = "SGLT2i", "TZD" = "TZD", "SU" = "SU"))) +
    theme_bw() +
    guides(fill = guide_legend(nrow = 1, byrow = TRUE)) +
    theme(
      legend.position = "bottom",
      axis.text.x = element_blank(),
      axis.title.x = element_blank(),
      strip.text.x = element_blank()
    ),
  
  # HbA1c
  discontinuation_per_drug_not_combined %>%
    gather("key", "value", -drugclass, -analysis_var, -group, -variable) %>%
    mutate(
      key = factor(key, levels = c("perc_12m", "perc_6m", "perc_3m"))
    ) %>%
    filter(variable == "hba1c" & key == "perc_6m") %>%
    ggplot(aes(y = group, x = value, fill = drugclass)) +
    geom_col(position = "dodge") +
    scale_fill_manual(values = c("Pooled" = "black", "SGLT2" = "#E69F00", "GLP1" = "#56B4E9", "SU" = "#CC79A7", "DPP4" = "#0072B2", "TZD" = "#D55E00", "MFN" = "grey"), breaks = c("GLP1", "DPP4", "SGLT2", "TZD", "SU"), labels = c("GLP-1RA", "DPP4i", "SGLT2i", "TZD", "SU"), name = " ", guide = guide_legend(reverse = TRUE)) +
    scale_y_discrete(limits = rev, labels = c("86+", "75-86", "64-75", "53-64")) +
    ggtitle("HbA1c, mmol/mol") +
    scale_x_continuous("Proportion of discontinuation", labels = scales::percent, limits = c(0, 0.5), breaks = seq(0, 0.5, 0.1)) +
    facet_wrap(~drugclass, nrow = 1, labeller = as_labeller(c("GLP1" = "GLP-1RA", "DPP4" = "DPP4i", "SGLT2" = "SGLT2i", "TZD" = "TZD", "SU" = "SU"))) +
    theme_bw() +
    guides(fill = guide_legend(nrow = 1, byrow = TRUE)) +
    theme(
      legend.position = "bottom",
      axis.text.x = element_blank(),
      axis.title.x = element_blank(),
      strip.text.x = element_blank()
    ),
  
  # eGFR
  discontinuation_per_drug_not_combined %>%
    gather("key", "value", -drugclass, -analysis_var, -group, -variable) %>%
    mutate(
      key = factor(key, levels = c("perc_12m", "perc_6m", "perc_3m"))
    ) %>%
    filter(variable == "preegfr" & key == "perc_6m") %>%
    ggplot(aes(y = group, x = value, fill = drugclass)) +
    geom_col(position = "dodge") +
    scale_fill_manual(values = c("Pooled" = "black", "SGLT2" = "#E69F00", "GLP1" = "#56B4E9", "SU" = "#CC79A7", "DPP4" = "#0072B2", "TZD" = "#D55E00", "MFN" = "grey"), breaks = c("GLP1", "DPP4", "SGLT2", "TZD", "SU"), labels = c("GLP-1RA", "DPP4i", "SGLT2i", "TZD", "SU"), name = " ", guide = guide_legend(reverse = TRUE)) +
    scale_y_discrete(limits = rev, labels = c("90+", "60-90", "<60")) +
    ggtitle("eGFR, ml/min per 1.73 m2") +
    scale_x_continuous("Proportion of discontinuation", labels = scales::percent, limits = c(0, 0.5), breaks = seq(0, 0.5, 0.1)) +
    facet_wrap(~drugclass, nrow = 1, labeller = as_labeller(c("GLP1" = "GLP-1RA", "DPP4" = "DPP4i", "SGLT2" = "SGLT2i", "TZD" = "TZD", "SU" = "SU"))) +
    theme_bw() +
    guides(fill = guide_legend(nrow = 1, byrow = TRUE)) +
    theme(
      legend.position = "bottom",
      axis.text.x = element_blank(),
      axis.title.x = element_blank(),
      strip.text.x = element_blank()
    ),
  
  # Statins
  discontinuation_per_drug_not_combined %>%
    gather("key", "value", -drugclass, -analysis_var, -group, -variable) %>%
    mutate(
      key = factor(key, levels = c("perc_12m", "perc_6m", "perc_3m"))
    ) %>%
    filter(variable == "statins" & key == "perc_6m") %>%
    ggplot(aes(y = group, x = value, fill = drugclass)) +
    geom_col(position = "dodge") +
    scale_fill_manual(values = c("Pooled" = "black", "SGLT2" = "#E69F00", "GLP1" = "#56B4E9", "SU" = "#CC79A7", "DPP4" = "#0072B2", "TZD" = "#D55E00", "MFN" = "grey"), breaks = c("GLP1", "DPP4", "SGLT2", "TZD", "SU"), labels = c("GLP-1RA", "DPP4i", "SGLT2i", "TZD", "SU"), name = " ", guide = guide_legend(reverse = TRUE)) +
    scale_y_discrete(limits = rev, breaks = c(0, 1), labels = c("No", "Yes")) +
    ggtitle("Statins") +
    scale_x_continuous("Proportion of discontinuation", labels = scales::percent, limits = c(0, 0.5), breaks = seq(0, 0.5, 0.1)) +
    facet_wrap(~drugclass, nrow = 1, labeller = as_labeller(c("GLP1" = "GLP-1RA", "DPP4" = "DPP4i", "SGLT2" = "SGLT2i", "TZD" = "TZD", "SU" = "SU"))) +
    theme_bw() +
    guides(fill = guide_legend(nrow = 1, byrow = TRUE)) +
    theme(
      legend.position = "bottom",
      axis.text.x = element_blank(),
      axis.title.x = element_blank(),
      strip.text.x = element_blank()
    ),
  
  # Blood pressure medication
  discontinuation_per_drug_not_combined %>%
    gather("key", "value", -drugclass, -analysis_var, -group, -variable) %>%
    mutate(
      key = factor(key, levels = c("perc_12m", "perc_6m", "perc_3m"))
    ) %>%
    filter(variable == "statins" & key == "perc_6m") %>%
    ggplot(aes(y = group, x = value, fill = drugclass)) +
    geom_col(position = "dodge") +
    scale_fill_manual(values = c("Pooled" = "black", "SGLT2" = "#E69F00", "GLP1" = "#56B4E9", "SU" = "#CC79A7", "DPP4" = "#0072B2", "TZD" = "#D55E00", "MFN" = "grey"), breaks = c("GLP1", "DPP4", "SGLT2", "TZD", "SU"), labels = c("GLP-1RA", "DPP4i", "SGLT2i", "TZD", "SU"), name = " ", guide = guide_legend(reverse = TRUE)) +
    scale_y_discrete(limits = rev, breaks = c(0, 1), labels = c("No", "Yes")) +
    ggtitle("Blood pressure medication") +
    scale_x_continuous("Proportion of discontinuation", labels = scales::percent, limits = c(0, 0.5), breaks = seq(0, 0.5, 0.1)) +
    facet_wrap(~drugclass, nrow = 1, labeller = as_labeller(c("GLP1" = "GLP-1RA", "DPP4" = "DPP4i", "SGLT2" = "SGLT2i", "TZD" = "TZD", "SU" = "SU"))) +
    theme_bw() +
    guides(fill = guide_legend(nrow = 1, byrow = TRUE)) +
    theme(
      legend.position = "bottom",
      axis.text.x = element_blank(),
      axis.title.x = element_blank(),
      strip.text.x = element_blank()
    ),
  
  # Frailty
  discontinuation_per_drug_not_combined %>%
    gather("key", "value", -drugclass, -analysis_var, -group, -variable) %>%
    mutate(
      key = factor(key, levels = c("perc_12m", "perc_6m", "perc_3m"))
    ) %>%
    filter(variable == "frailty" & key == "perc_6m") %>%
    ggplot(aes(y = group, x = value, fill = drugclass)) +
    geom_col(position = "dodge") +
    scale_fill_manual(values = c("Pooled" = "black", "SGLT2" = "#E69F00", "GLP1" = "#56B4E9", "SU" = "#CC79A7", "DPP4" = "#0072B2", "TZD" = "#D55E00", "MFN" = "grey"), breaks = c("GLP1", "DPP4", "SGLT2", "TZD", "SU"), labels = c("GLP-1RA", "DPP4i", "SGLT2i", "TZD", "SU"), name = " ", guide = guide_legend(reverse = TRUE)) +
    scale_y_discrete(limits = rev, breaks = c(0, 1), labels = c("No", "Yes")) +
    ggtitle("Frailty") +
    scale_x_continuous("Proportion of discontinuation", labels = scales::percent, limits = c(0, 0.5), breaks = seq(0, 0.5, 0.1)) +
    facet_wrap(~drugclass, nrow = 1, labeller = as_labeller(c("GLP1" = "GLP-1RA", "DPP4" = "DPP4i", "SGLT2" = "SGLT2i", "TZD" = "TZD", "SU" = "SU"))) +
    theme_bw() +
    guides(fill = guide_legend(nrow = 1, byrow = TRUE)) +
    theme(
      legend.position = "bottom",
      axis.text.x = element_blank(),
      axis.title.x = element_blank(),
      strip.text.x = element_blank()
    ),
  
  # stopping metformin
  discontinuation_per_drug_not_combined %>%
    gather("key", "value", -drugclass, -analysis_var, -group, -variable) %>%
    mutate(
      key = factor(key, levels = c("perc_12m", "perc_6m", "perc_3m"))
    ) %>%
    filter(variable == "stopdrug_MFN" & key == "perc_6m") %>%
    ggplot(aes(y = group, x = value, fill = drugclass)) +
    geom_col() +
    scale_fill_manual(values = c("Pooled" = "black", "SGLT2" = "#E69F00", "GLP1" = "#56B4E9", "SU" = "#CC79A7", "DPP4" = "#0072B2", "TZD" = "#D55E00", "MFN" = "grey"), breaks = c("GLP1", "DPP4", "SGLT2", "TZD", "SU"), labels = c("GLP-1RA", "DPP4i", "SGLT2i", "TZD", "SU"), name = " ", guide = guide_legend(reverse = TRUE)) +
    scale_y_discrete(limits = rev, breaks = c(0, 1), labels = c("No", "Yes")) +
    ggtitle("History of Metformin discontinuation within 3-months") +
    scale_x_continuous("Proportion of discontinuation", labels = scales::percent, limits = c(0, 0.5), breaks = seq(0, 0.5, 0.1)) +
    facet_wrap(~drugclass, nrow = 1, labeller = as_labeller(c("GLP1" = "GLP-1RA", "DPP4" = "DPP4i", "SGLT2" = "SGLT2i", "TZD" = "TZD", "SU" = "SU"))) +
    theme_bw() +
    guides(fill = guide_legend(nrow = 1, byrow = TRUE)) +
    theme(
      legend.position = "bottom",
      axis.text.x = element_blank(),
      axis.title.x = element_blank(),
      strip.text.x = element_blank()
    ),
  
  # numdrugs
  discontinuation_per_drug_not_combined %>%
    gather("key", "value", -drugclass, -analysis_var, -group, -variable) %>%
    mutate(
      key = factor(key, levels = c("perc_12m", "perc_6m", "perc_3m"))
    ) %>%
    filter(variable == "numdrugs" & key == "perc_6m") %>%
    ggplot(aes(y = group, x = value, fill = drugclass)) +
    geom_col() +
    scale_fill_manual(values = c("Pooled" = "black", "SGLT2" = "#E69F00", "GLP1" = "#56B4E9", "SU" = "#CC79A7", "DPP4" = "#0072B2", "TZD" = "#D55E00", "MFN" = "grey"), breaks = c("GLP1", "DPP4", "SGLT2", "TZD", "SU"), labels = c("GLP-1RA", "DPP4i", "SGLT2i", "TZD", "SU"), name = " ", guide = guide_legend(reverse = TRUE)) +
    scale_y_discrete(limits = rev, labels = c("2+", "1", "0")) +
    ggtitle("Number of other current glucose-lowering therapies") +
    scale_x_continuous("Proportion of discontinuation", labels = scales::percent, limits = c(0, 0.5), breaks = seq(0, 0.5, 0.1)) +
    facet_wrap(~drugclass, nrow = 1, labeller = as_labeller(c("GLP1" = "GLP-1RA", "DPP4" = "DPP4i", "SGLT2" = "SGLT2i", "TZD" = "TZD", "SU" = "SU"))) +
    theme_bw() +
    guides(fill = guide_legend(nrow = 1, byrow = TRUE)) +
    theme(
      legend.position = "bottom",
      axis.text.x = element_blank(),
      axis.title.x = element_blank(),
      strip.text.x = element_blank()
    ),
  
  # drugline
  discontinuation_per_drug_not_combined %>%
    gather("key", "value", -drugclass, -analysis_var, -group, -variable) %>%
    mutate(
      key = factor(key, levels = c("perc_12m", "perc_6m", "perc_3m"))
    ) %>%
    filter(variable == "drugline" & key == "perc_6m") %>%
    ggplot(aes(y = group, x = value, fill = drugclass)) +
    geom_col(position = "dodge") +
    scale_fill_manual(values = c("Pooled" = "black", "SGLT2" = "#E69F00", "GLP1" = "#56B4E9", "SU" = "#CC79A7", "DPP4" = "#0072B2", "TZD" = "#D55E00", "MFN" = "grey"), breaks = c("GLP1", "DPP4", "SGLT2", "TZD", "SU"), labels = c("GLP-1RA", "DPP4i", "SGLT2i", "TZD", "SU"), name = " ", guide = guide_legend(reverse = TRUE)) +
    scale_y_discrete(limits = rev) +
    ggtitle("Number of therapy classes ever prescribed") +
    scale_x_continuous("Proportion of discontinuation", labels = scales::percent, limits = c(0, 0.5), breaks = seq(0, 0.5, 0.1)) +
    facet_wrap(~drugclass, nrow = 1, labeller = as_labeller(c("GLP1" = "GLP-1RA", "DPP4" = "DPP4i", "SGLT2" = "SGLT2i", "TZD" = "TZD", "SU" = "SU"))) +
    theme_bw() +
    guides(fill = guide_legend(nrow = 1, byrow = TRUE)) +
    theme(
      legend.position = "bottom",
      axis.text.x = element_text(size = 10),
      axis.title.x = element_text(size = 11),
      strip.text.x = element_blank()
    )
  
  , ncol = 1
  
) + 
  plot_layout(
    guides = "collect",
    height = c(4, 5, 2, 3, 4, 5, 4, 4, 3, 2, 2, 2, 2, 3, 3),
    # axis_titles = "collect" # not working, manually done in the code
  ) +
  plot_layout(
    axis_titles = "collect"
  ) &
  theme(
    legend.position = "none",
    plot.title.position = "plot",
    axis.title.y = element_blank(),
    axis.text.y = element_text(size = 10),
    panel.spacing = unit(0.3, "cm", data = NULL)
  )

dev.off()



# load dataset
cprd_dataset <- set_up_data(
  raw_data = "20240814_t2d_1stinstance",
  diagnosis = FALSE,
  therapies = c("DPP4", "GLP1", "SGLT2", "SU", "TZD"),
  dataset = "6m.disc.dataset",
  follow_up = "6-months",
  full_prescribing_history = TRUE
) %>%
  drop_na(-prehdl)




# Start plot by variables


pdf("results/figures/03.12m_discontinuation_by_variable.pdf", width = 7, height = 15)

patchwork::wrap_plots(
  
  # Overall
  discontinuation_per_drug_combined %>%
    gather("key", "value", -drugclass, -analysis_var, -group, -variable) %>%
    mutate(
      key = factor(key, levels = c("perc_12m", "perc_6m", "perc_3m"))
    ) %>%
    filter(variable == "overall" & key == "perc_12m") %>%
    ggplot(aes(y = analysis_var, x = value, fill = drugclass)) +
    geom_col() +
    scale_fill_manual(values = c("other" = "black", "SGLT2" = "#E69F00", "GLP1" = "#56B4E9", "SU" = "#CC79A7", "DPP4" = "#0072B2", "TZD" = "#D55E00", "MFN" = "grey"), breaks = c("other"), labels = c("Other glucose-lowering therapies"), name = " ", guide = guide_legend(reverse = FALSE)) +
    scale_y_discrete(breaks = c(1), labels = c("Overall")) +
    scale_x_continuous("Proportion of discontinuation", labels = scales::percent, limits = c(0, 0.3)) +
    facet_wrap(~drugclass, nrow = 1, labeller = as_labeller(c("other" = "Other glucose-lowering therapies"))) +
    theme_bw() +
    guides(fill = guide_legend(nrow = 1, byrow = TRUE)) +
    theme(
      legend.position = "bottom",
      axis.text.x = element_blank(),
      axis.title.x = element_blank()
    ),
  
  # Age
  discontinuation_per_drug_combined %>%
    gather("key", "value", -drugclass, -analysis_var, -group, -variable) %>%
    mutate(
      key = factor(key, levels = c("perc_12m", "perc_6m", "perc_3m"))
    ) %>%
    filter(variable == "dstartdate_age" & key == "perc_12m") %>%
    ggplot(aes(y = group, x = value, fill = drugclass)) +
    geom_col() +
    scale_fill_manual(values = c("other" = "black", "SGLT2" = "#E69F00", "GLP1" = "#56B4E9", "SU" = "#CC79A7", "DPP4" = "#0072B2", "TZD" = "#D55E00", "MFN" = "grey"), breaks = c("other"), labels = c("Other glucose-lowering therapies"), name = " ", guide = guide_legend(reverse = FALSE)) +
    scale_y_discrete(limits = rev, labels = c("70+", "60-69", "50-59", "<50")) +
    scale_x_continuous("Proportion of discontinuation", labels = scales::percent, limits = c(0, 0.3)) +
    ggtitle("Age group, years") +
    facet_wrap(~drugclass, nrow = 1, labeller = as_labeller(c("other" = "Other glucose-lowering therapies"))) +
    theme_bw() +
    guides(fill = guide_legend(nrow = 1, byrow = TRUE)) +
    theme(
      plot.title.position = "plot",
      legend.position = "bottom",
      axis.text.x = element_blank(),
      axis.title.x = element_blank(),
      strip.text.x = element_blank()
    ),
  
  # Diabetes duration
  discontinuation_per_drug_combined %>%
    gather("key", "value", -drugclass, -analysis_var, -group, -variable) %>%
    mutate(
      key = factor(key, levels = c("perc_12m", "perc_6m", "perc_3m"))
    ) %>%
    filter(variable == "dstartdate_dm_dur" & key == "perc_12m") %>%
    ggplot(aes(y = group, x = value, fill = drugclass)) +
    geom_col() +
    scale_fill_manual(values = c("other" = "black", "SGLT2" = "#E69F00", "GLP1" = "#56B4E9", "SU" = "#CC79A7", "DPP4" = "#0072B2", "TZD" = "#D55E00", "MFN" = "grey"), breaks = c("other"), labels = c("Other glucose-lowering therapies"), name = " ", guide = guide_legend(reverse = FALSE)) +
    scale_y_discrete(limits = rev, labels = c("10+", "6-9", "3-5", "1-2", "<1")) +
    ggtitle("Duration of diabetes, years") +
    scale_x_continuous("Proportion of discontinuation", labels = scales::percent, limits = c(0, 0.3)) +
    facet_wrap(~drugclass, nrow = 1, labeller = as_labeller(c("other" = "Other glucose-lowering therapies"))) +
    theme_bw() +
    guides(fill = guide_legend(nrow = 1, byrow = TRUE)) +
    theme(
      plot.title.position = "plot",
      legend.position = "bottom",
      axis.text.x = element_blank(),
      axis.title.x = element_blank(),
      strip.text.x = element_blank()
    ),
  
  # Sex
  discontinuation_per_drug_combined %>%
    gather("key", "value", -drugclass, -analysis_var, -group, -variable) %>%
    mutate(
      key = factor(key, levels = c("perc_12m", "perc_6m", "perc_3m"))
    ) %>%
    filter(variable == "gender" & key == "perc_12m") %>%
    ggplot(aes(y = group, x = value, fill = drugclass)) +
    geom_col() +
    scale_fill_manual(values = c("other" = "black", "SGLT2" = "#E69F00", "GLP1" = "#56B4E9", "SU" = "#CC79A7", "DPP4" = "#0072B2", "TZD" = "#D55E00", "MFN" = "grey"), breaks = c("other"), labels = c("Other glucose-lowering therapies"), name = " ", guide = guide_legend(reverse = FALSE)) +
    scale_y_discrete(limits = rev, labels = c("Female", "Male")) +
    scale_x_continuous("Proportion of discontinuation", labels = scales::percent, limits = c(0, 0.3)) +
    ggtitle("Sex") +
    facet_wrap(~drugclass, nrow = 1, labeller = as_labeller(c("other" = "Other glucose-lowering therapies"))) +
    theme_bw() +
    guides(fill = guide_legend(nrow = 1, byrow = TRUE)) +
    theme(
      plot.title.position = "plot",
      legend.position = "bottom",
      axis.text.x = element_blank(),
      axis.title.x = element_blank(),
      strip.text.x = element_blank()
    ),
  
  # Smoking
  discontinuation_per_drug_combined %>%
    gather("key", "value", -drugclass, -analysis_var, -group, -variable) %>%
    mutate(
      key = factor(key, levels = c("perc_12m", "perc_6m", "perc_3m"))
    ) %>%
    filter(variable == "smoking" & key == "perc_12m") %>%
    ggplot(aes(y = group, x = value, fill = drugclass)) +
    geom_col() +
    scale_fill_manual(values = c("other" = "black", "SGLT2" = "#E69F00", "GLP1" = "#56B4E9", "SU" = "#CC79A7", "DPP4" = "#0072B2", "TZD" = "#D55E00", "MFN" = "grey"), breaks = c("other"), labels = c("Other glucose-lowering therapies"), name = " ", guide = guide_legend(reverse = FALSE)) +
    ggtitle("Smoking") +
    scale_x_continuous("Proportion of discontinuation", labels = scales::percent, limits = c(0, 0.3)) +
    facet_wrap(~drugclass, nrow = 1, labeller = as_labeller(c("other" = "Other glucose-lowering therapies"))) +
    theme_bw() +
    guides(fill = guide_legend(nrow = 1, byrow = TRUE)) +
    theme(
      plot.title.position = "plot",
      legend.position = "bottom",
      axis.text.x = element_blank(),
      axis.title.x = element_blank(),
      strip.text.x = element_blank()
    ),
  
  # Ethnicity
  discontinuation_per_drug_combined %>%
    gather("key", "value", -drugclass, -analysis_var, -group, -variable) %>%
    mutate(
      key = factor(key, levels = c("perc_12m", "perc_6m", "perc_3m"))
    ) %>%
    filter(variable == "ethnicity" & key == "perc_12m") %>%
    ggplot(aes(y = group, x = value, fill = drugclass)) +
    geom_col() +
    scale_fill_manual(values = c("other" = "black", "SGLT2" = "#E69F00", "GLP1" = "#56B4E9", "SU" = "#CC79A7", "DPP4" = "#0072B2", "TZD" = "#D55E00", "MFN" = "grey"), breaks = c("other"), labels = c("Other glucose-lowering therapies"), name = " ", guide = guide_legend(reverse = FALSE)) +
    ggtitle("Ethnicity") +
    scale_x_continuous("Proportion of discontinuation", labels = scales::percent, limits = c(0, 0.3)) +
    facet_wrap(~drugclass, nrow = 1, labeller = as_labeller(c("other" = "Other glucose-lowering therapies"))) +
    theme_bw() +
    guides(fill = guide_legend(nrow = 1, byrow = TRUE)) +
    theme(
      plot.title.position = "plot",
      legend.position = "bottom",
      axis.text.x = element_blank(),
      axis.title.x = element_blank(),
      strip.text.x = element_blank()
    ),
  
  # Quintile
  discontinuation_per_drug_combined %>%
    gather("key", "value", -drugclass, -analysis_var, -group, -variable) %>%
    mutate(
      key = factor(key, levels = c("perc_12m", "perc_6m", "perc_3m"))
    ) %>%
    filter(variable == "imd2015" & key == "perc_12m") %>%
    ggplot(aes(y = group, x = value, fill = drugclass)) +
    geom_col() +
    scale_fill_manual(values = c("other" = "black", "SGLT2" = "#E69F00", "GLP1" = "#56B4E9", "SU" = "#CC79A7", "DPP4" = "#0072B2", "TZD" = "#D55E00", "MFN" = "grey"), breaks = c("other"), labels = c("Other glucose-lowering therapies"), name = " ", guide = guide_legend(reverse = FALSE)) +
    scale_y_discrete(limits = rev, labels = c("5 (most deprived)", "4", "3", "2", "1 (least deprived)")) +
    scale_x_continuous("Proportion of discontinuation", labels = scales::percent, limits = c(0, 0.3)) +
    ggtitle("Index of multiple deprivation") +
    facet_wrap(~drugclass, nrow = 1, labeller = as_labeller(c("other" = "Other glucose-lowering therapies"))) +
    theme_bw() +
    guides(fill = guide_legend(nrow = 1, byrow = TRUE)) +
    theme(
      plot.title.position = "plot",
      legend.position = "bottom",
      axis.text.x = element_blank(),
      axis.title.x = element_blank(),
      strip.text.x = element_blank()
    ),
  
  # BMI
  discontinuation_per_drug_combined %>%
    gather("key", "value", -drugclass, -analysis_var, -group, -variable) %>%
    mutate(
      key = factor(key, levels = c("perc_12m", "perc_6m", "perc_3m"))
    ) %>%
    filter(variable == "bmi" & key == "perc_12m") %>%
    ggplot(aes(y = group, x = value, fill = drugclass)) +
    geom_col() +
    scale_fill_manual(values = c("other" = "black", "SGLT2" = "#E69F00", "GLP1" = "#56B4E9", "SU" = "#CC79A7", "DPP4" = "#0072B2", "TZD" = "#D55E00", "MFN" = "grey"), breaks = c("other"), labels = c("Other glucose-lowering therapies"), name = " ", guide = guide_legend(reverse = FALSE)) +
    scale_y_discrete(limits = rev, labels = c("35+", "30-34.9", "25-29.9", "<25")) +
    ggtitle("BMI, kg/m2") +
    scale_x_continuous("Proportion of discontinuation", labels = scales::percent, limits = c(0, 0.3)) +
    facet_wrap(~drugclass, nrow = 1, labeller = as_labeller(c("other" = "Other glucose-lowering therapies"))) +
    theme_bw() +
    guides(fill = guide_legend(nrow = 1, byrow = TRUE)) +
    theme(
      plot.title.position = "plot",
      legend.position = "bottom",
      axis.text.x = element_blank(),
      axis.title.x = element_blank(),
      strip.text.x = element_blank()
    ),
  
  # HbA1c
  discontinuation_per_drug_combined %>%
    gather("key", "value", -drugclass, -analysis_var, -group, -variable) %>%
    mutate(
      key = factor(key, levels = c("perc_12m", "perc_6m", "perc_3m"))
    ) %>%
    filter(variable == "hba1c" & key == "perc_12m") %>%
    ggplot(aes(y = group, x = value, fill = drugclass)) +
    geom_col() +
    scale_fill_manual(values = c("other" = "black", "SGLT2" = "#E69F00", "GLP1" = "#56B4E9", "SU" = "#CC79A7", "DPP4" = "#0072B2", "TZD" = "#D55E00", "MFN" = "grey"), breaks = c("other"), labels = c("Other glucose-lowering therapies"), name = " ", guide = guide_legend(reverse = FALSE)) +
    scale_y_discrete(limits = rev, labels = c("86+", "75-86", "64-75", "53-64")) +
    ggtitle("HbA1c, mmol/mol") +
    scale_x_continuous("Proportion of discontinuation", labels = scales::percent, limits = c(0, 0.3)) +
    facet_wrap(~drugclass, nrow = 1, labeller = as_labeller(c("other" = "Other glucose-lowering therapies"))) +
    theme_bw() +
    guides(fill = guide_legend(nrow = 1, byrow = TRUE)) +
    theme(
      plot.title.position = "plot",
      legend.position = "bottom",
      axis.text.x = element_blank(),
      axis.title.x = element_blank(),
      strip.text.x = element_blank()
    ),
  
  # eGFR
  discontinuation_per_drug_combined %>%
    gather("key", "value", -drugclass, -analysis_var, -group, -variable) %>%
    mutate(
      key = factor(key, levels = c("perc_12m", "perc_6m", "perc_3m"))
    ) %>%
    filter(variable == "preegfr" & key == "perc_12m") %>%
    ggplot(aes(y = group, x = value, fill = drugclass)) +
    geom_col() +
    scale_fill_manual(values = c("other" = "black", "SGLT2" = "#E69F00", "GLP1" = "#56B4E9", "SU" = "#CC79A7", "DPP4" = "#0072B2", "TZD" = "#D55E00", "MFN" = "grey"), breaks = c("other"), labels = c("Other glucose-lowering therapies"), name = " ", guide = guide_legend(reverse = FALSE)) +
    scale_y_discrete(limits = rev, labels = c("90+", "60-90", "<60")) +
    ggtitle("eGFR, ml/min per 1.73 m2") +
    scale_x_continuous("Proportion of discontinuation", labels = scales::percent, limits = c(0, 0.3)) +
    facet_wrap(~drugclass, nrow = 1, labeller = as_labeller(c("other" = "Other glucose-lowering therapies"))) +
    theme_bw() +
    guides(fill = guide_legend(nrow = 1, byrow = TRUE)) +
    theme(
      plot.title.position = "plot",
      legend.position = "bottom",
      axis.text.x = element_blank(),
      axis.title.x = element_blank(),
      strip.text.x = element_blank()
    ),
  
  # Statins
  discontinuation_per_drug_combined %>%
    gather("key", "value", -drugclass, -analysis_var, -group, -variable) %>%
    mutate(
      key = factor(key, levels = c("perc_12m", "perc_6m", "perc_3m"))
    ) %>%
    filter(variable == "statins" & key == "perc_12m") %>%
    ggplot(aes(y = group, x = value, fill = drugclass)) +
    geom_col() +
    scale_fill_manual(values = c("other" = "black", "SGLT2" = "#E69F00", "GLP1" = "#56B4E9", "SU" = "#CC79A7", "DPP4" = "#0072B2", "TZD" = "#D55E00", "MFN" = "grey"), breaks = c("other"), labels = c("Other glucose-lowering therapies"), name = " ", guide = guide_legend(reverse = FALSE)) +
    scale_y_discrete(limits = rev, breaks = c(0, 1), labels = c("No", "Yes")) +
    ggtitle("Statins") +
    scale_x_continuous("Proportion of discontinuation", labels = scales::percent, limits = c(0, 0.3)) +
    facet_wrap(~drugclass, nrow = 1, labeller = as_labeller(c("other" = "Other glucose-lowering therapies"))) +
    theme_bw() +
    guides(fill = guide_legend(nrow = 1, byrow = TRUE)) +
    theme(
      plot.title.position = "plot",
      legend.position = "bottom",
      axis.text.x = element_blank(),
      axis.title.x = element_blank(),
      strip.text.x = element_blank()
    ),
  
  # Bloodmed
  discontinuation_per_drug_combined %>%
    gather("key", "value", -drugclass, -analysis_var, -group, -variable) %>%
    mutate(
      key = factor(key, levels = c("perc_12m", "perc_6m", "perc_3m"))
    ) %>%
    filter(variable == "bloodmed" & key == "perc_12m") %>%
    ggplot(aes(y = group, x = value, fill = drugclass)) +
    geom_col() +
    scale_fill_manual(values = c("other" = "black", "SGLT2" = "#E69F00", "GLP1" = "#56B4E9", "SU" = "#CC79A7", "DPP4" = "#0072B2", "TZD" = "#D55E00", "MFN" = "grey"), breaks = c("other"), labels = c("Other glucose-lowering therapies"), name = " ", guide = guide_legend(reverse = FALSE)) +
    scale_y_discrete(limits = rev, breaks = c(0, 1), labels = c("No", "Yes")) +
    ggtitle("Blood pressure medication") +
    scale_x_continuous("Proportion of discontinuation", labels = scales::percent, limits = c(0, 0.3)) +
    facet_wrap(~drugclass, nrow = 1, labeller = as_labeller(c("other" = "Other glucose-lowering therapies"))) +
    theme_bw() +
    guides(fill = guide_legend(nrow = 1, byrow = TRUE)) +
    theme(
      plot.title.position = "plot",
      legend.position = "bottom",
      axis.text.x = element_blank(),
      axis.title.x = element_blank(),
      strip.text.x = element_blank()
    ),
  
  # Frailty
  discontinuation_per_drug_combined %>%
    gather("key", "value", -drugclass, -analysis_var, -group, -variable) %>%
    mutate(
      key = factor(key, levels = c("perc_12m", "perc_6m", "perc_3m"))
    ) %>%
    filter(variable == "frailty" & key == "perc_12m") %>%
    ggplot(aes(y = group, x = value, fill = drugclass)) +
    geom_col() +
    scale_fill_manual(values = c("other" = "black", "SGLT2" = "#E69F00", "GLP1" = "#56B4E9", "SU" = "#CC79A7", "DPP4" = "#0072B2", "TZD" = "#D55E00", "MFN" = "grey"), breaks = c("other"), labels = c("Other glucose-lowering therapies"), name = " ", guide = guide_legend(reverse = FALSE)) +
    scale_y_discrete(limits = rev, breaks = c(0, 1), labels = c("No", "Yes")) +
    ggtitle("Frailty") +
    scale_x_continuous("Proportion of discontinuation", labels = scales::percent, limits = c(0, 0.3)) +
    facet_wrap(~drugclass, nrow = 1, labeller = as_labeller(c("other" = "Other glucose-lowering therapies"))) +
    theme_bw() +
    guides(fill = guide_legend(nrow = 1, byrow = TRUE)) +
    theme(
      plot.title.position = "plot",
      legend.position = "bottom",
      axis.text.x = element_blank(),
      axis.title.x = element_blank(),
      strip.text.x = element_blank()
    ),
  
  # stopping metformin
  discontinuation_per_drug_combined %>%
    gather("key", "value", -drugclass, -analysis_var, -group, -variable) %>%
    mutate(
      key = factor(key, levels = c("perc_12m", "perc_6m", "perc_3m"))
    ) %>%
    filter(variable == "stopdrug_MFN" & key == "perc_12m") %>%
    ggplot(aes(y = group, x = value, fill = drugclass)) +
    geom_col() +
    scale_fill_manual(values = c("other" = "black", "SGLT2" = "#E69F00", "GLP1" = "#56B4E9", "SU" = "#CC79A7", "DPP4" = "#0072B2", "TZD" = "#D55E00", "MFN" = "grey"), breaks = c("other"), labels = c("Other glucose-lowering therapies"), name = " ", guide = guide_legend(reverse = FALSE)) +
    scale_y_discrete(limits = rev, breaks = c(0, 1), labels = c("No", "Yes")) +
    ggtitle("History of Metformin discontinuation within 3-months") +
    scale_x_continuous("Proportion of discontinuation", labels = scales::percent, limits = c(0, 0.3)) +
    facet_wrap(~drugclass, nrow = 1, labeller = as_labeller(c("other" = "Other glucose-lowering therapies"))) +
    theme_bw() +
    guides(fill = guide_legend(nrow = 1, byrow = TRUE)) +
    theme(
      plot.title.position = "plot",
      legend.position = "bottom",
      axis.text.x = element_blank(),
      axis.title.x = element_blank(),
      strip.text.x = element_blank()
    ),
  
  # Numdrugs
  discontinuation_per_drug_combined %>%
    gather("key", "value", -drugclass, -analysis_var, -group, -variable) %>%
    mutate(
      key = factor(key, levels = c("perc_12m", "perc_6m", "perc_3m"))
    ) %>%
    filter(variable == "numdrugs" & key == "perc_12m") %>%
    ggplot(aes(y = group, x = value, fill = drugclass)) +
    geom_col() +
    scale_fill_manual(values = c("other" = "black", "SGLT2" = "#E69F00", "GLP1" = "#56B4E9", "SU" = "#CC79A7", "DPP4" = "#0072B2", "TZD" = "#D55E00", "MFN" = "grey"), breaks = c("other"), labels = c("Other glucose-lowering therapies"), name = " ", guide = guide_legend(reverse = FALSE)) +
    scale_y_discrete(limits = rev, labels = c("2+", "1", "0")) +
    ggtitle("Number of other current glucose-lowering therapies") +
    scale_x_continuous("Proportion of discontinuation", labels = scales::percent, limits = c(0, 0.3)) +
    facet_wrap(~drugclass, nrow = 1, labeller = as_labeller(c("other" = "Other glucose-lowering therapies"))) +
    theme_bw() +
    guides(fill = guide_legend(nrow = 1, byrow = TRUE)) +
    theme(
      plot.title.position = "plot",
      legend.position = "bottom",
      axis.text.x = element_blank(),
      axis.title.x = element_blank(),
      strip.text.x = element_blank()
    ),
  
  # drugline
  discontinuation_per_drug_combined %>%
    gather("key", "value", -drugclass, -analysis_var, -group, -variable) %>%
    mutate(
      key = factor(key, levels = c("perc_12m", "perc_6m", "perc_3m"))
    ) %>%
    filter(variable == "drugline" & key == "perc_12m") %>%
    ggplot(aes(y = group, x = value, fill = drugclass)) +
    geom_col() +
    scale_fill_manual(values = c("other" = "black", "SGLT2" = "#E69F00", "GLP1" = "#56B4E9", "SU" = "#CC79A7", "DPP4" = "#0072B2", "TZD" = "#D55E00", "MFN" = "grey"), breaks = c("other"), labels = c("Other glucose-lowering therapies"), name = " ", guide = guide_legend(reverse = FALSE)) +
    scale_y_discrete(limits = rev) +
    ggtitle("Number of therapy classes ever prescribed") +
    scale_x_continuous("Proportion of discontinuation", labels = scales::percent, limits = c(0, 0.3)) +
    facet_wrap(~drugclass, nrow = 1, labeller = as_labeller(c("other" = "Other glucose-lowering therapies"))) +
    theme_bw() +
    guides(fill = guide_legend(nrow = 1, byrow = TRUE)) +
    theme(
      plot.title.position = "plot",
      legend.position = "bottom",
      axis.text.x = element_text(size = 10),
      axis.title.x = element_text(size = 11),
      strip.text.x = element_blank()
    )
  
  , ncol = 1
  
) + 
  plot_layout(
    guides = "collect",
    height = c(1, 4, 5, 2, 3, 4, 5, 4, 4, 3, 2, 2, 2, 2, 3, 3),
    # axis_titles = "collect" # not working, manually done in the code
  ) +
  plot_layout(
    axis_titles = "collect"
  ) &
  theme(
    legend.position = "none",
    plot.title.position = "plot",
    axis.title.y = element_blank(),
    axis.text.y = element_text(size = 10),
    panel.spacing = unit(0.3, "cm", data = NULL)
  )


dev.off()


pdf("results/figures/03.12m_discontinuation_by_variable_by_drug.pdf", width = 14, height = 15)

patchwork::wrap_plots(
  
  # Age - Metformin
  discontinuation_per_drug_not_combined %>%
    gather("key", "value", -drugclass, -analysis_var, -group, -variable) %>%
    mutate(
      key = factor(key, levels = c("perc_12m", "perc_6m", "perc_3m"))
    ) %>%
    filter(variable == "dstartdate_age" & key == "perc_12m") %>%
    ggplot(aes(y = group, x = value, fill = drugclass)) +
    geom_col(position = "dodge") +
    scale_fill_manual(values = c("Pooled" = "black", "SGLT2" = "#E69F00", "GLP1" = "#56B4E9", "SU" = "#CC79A7", "DPP4" = "#0072B2", "TZD" = "#D55E00", "MFN" = "grey"), breaks = c("GLP1", "DPP4", "SGLT2", "TZD", "SU"), labels = c("GLP-1RA", "DPP4i", "SGLT2i", "TZD", "SU"), name = " ", guide = guide_legend(reverse = TRUE)) +
    scale_y_discrete(limits = rev, labels = c("70+", "60-69", "50-59", "<50")) +
    ggtitle("Age group, years") +
    scale_x_continuous("Proportion of discontinuation", labels = scales::percent, limits = c(0, 0.5), breaks = seq(0, 0.5, 0.1)) +
    facet_wrap(~drugclass, nrow = 1, labeller = as_labeller(c("GLP1" = "GLP-1RA", "DPP4" = "DPP4i", "SGLT2" = "SGLT2i", "TZD" = "TZD", "SU" = "SU"))) +
    theme_bw() +
    guides(fill = guide_legend(nrow = 1, byrow = TRUE)) +
    theme(
      legend.position = "bottom",
      axis.text.x = element_blank(),
      axis.title.x = element_blank()
    ),
  
  # Diabetes duration
  discontinuation_per_drug_not_combined %>%
    gather("key", "value", -drugclass, -analysis_var, -group, -variable) %>%
    mutate(
      key = factor(key, levels = c("perc_12m", "perc_6m", "perc_3m"))
    ) %>%
    filter(variable == "dstartdate_dm_dur" & key == "perc_12m") %>%
    ggplot(aes(y = group, x = value, fill = drugclass)) +
    geom_col(position = "dodge") +
    scale_fill_manual(values = c("Pooled" = "black", "SGLT2" = "#E69F00", "GLP1" = "#56B4E9", "SU" = "#CC79A7", "DPP4" = "#0072B2", "TZD" = "#D55E00", "MFN" = "grey"), breaks = c("GLP1", "DPP4", "SGLT2", "TZD", "SU"), labels = c("GLP-1RA", "DPP4i", "SGLT2i", "TZD", "SU"), name = " ", guide = guide_legend(reverse = TRUE)) +
    scale_y_discrete(limits = rev, labels = c("10+", "6-9", "3-5", "1-2", "<1")) +
    ggtitle("Duration of diabetes, years") +
    scale_x_continuous("Proportion of discontinuation", labels = scales::percent, limits = c(0, 0.5), breaks = seq(0, 0.5, 0.1)) +
    facet_wrap(~drugclass, nrow = 1, labeller = as_labeller(c("GLP1" = "GLP-1RA", "DPP4" = "DPP4i", "SGLT2" = "SGLT2i", "TZD" = "TZD", "SU" = "SU"))) +
    theme_bw() +
    guides(fill = guide_legend(nrow = 1, byrow = TRUE)) +
    theme(
      legend.position = "bottom",
      axis.text.x = element_blank(),
      axis.title.x = element_blank(),
      strip.text.x = element_blank()
    ),
  
  # Sex 
  discontinuation_per_drug_not_combined %>%
    gather("key", "value", -drugclass, -analysis_var, -group, -variable) %>%
    mutate(
      key = factor(key, levels = c("perc_12m", "perc_6m", "perc_3m"))
    ) %>%
    filter(variable == "gender" & key == "perc_12m") %>%
    ggplot(aes(y = group, x = value, fill = drugclass)) +
    geom_col(position = "dodge") +
    scale_fill_manual(values = c("Pooled" = "black", "SGLT2" = "#E69F00", "GLP1" = "#56B4E9", "SU" = "#CC79A7", "DPP4" = "#0072B2", "TZD" = "#D55E00", "MFN" = "grey"), breaks = c("GLP1", "DPP4", "SGLT2", "TZD", "SU"), labels = c("GLP-1RA", "DPP4i", "SGLT2i", "TZD", "SU"), name = " ", guide = guide_legend(reverse = TRUE)) +
    scale_y_discrete(limits = rev, labels = c("Female", "Male")) +
    ggtitle("Sex") +
    scale_x_continuous("Proportion of discontinuation", labels = scales::percent, limits = c(0, 0.5), breaks = seq(0, 0.5, 0.1)) +
    facet_wrap(~drugclass, nrow = 1, labeller = as_labeller(c("GLP1" = "GLP-1RA", "DPP4" = "DPP4i", "SGLT2" = "SGLT2i", "TZD" = "TZD", "SU" = "SU"))) +
    theme_bw() +
    guides(fill = guide_legend(nrow = 1, byrow = TRUE)) +
    theme(
      legend.position = "bottom",
      axis.text.x = element_blank(),
      axis.title.x = element_blank(),
      strip.text.x = element_blank()
    ),
  
  
  # Smoking 
  discontinuation_per_drug_not_combined %>%
    gather("key", "value", -drugclass, -analysis_var, -group, -variable) %>%
    mutate(
      key = factor(key, levels = c("perc_12m", "perc_6m", "perc_3m"))
    ) %>%
    filter(variable == "smoking" & key == "perc_12m") %>%
    ggplot(aes(y = group, x = value, fill = drugclass)) +
    geom_col(position = "dodge") +
    scale_fill_manual(values = c("Pooled" = "black", "SGLT2" = "#E69F00", "GLP1" = "#56B4E9", "SU" = "#CC79A7", "DPP4" = "#0072B2", "TZD" = "#D55E00", "MFN" = "grey"), breaks = c("GLP1", "DPP4", "SGLT2", "TZD", "SU"), labels = c("GLP-1RA", "DPP4i", "SGLT2i", "TZD", "SU"), name = " ", guide = guide_legend(reverse = TRUE)) +
    ggtitle("Smoking") +
    scale_x_continuous("Proportion of discontinuation", labels = scales::percent, limits = c(0, 0.5), breaks = seq(0, 0.5, 0.1)) +
    facet_wrap(~drugclass, nrow = 1, labeller = as_labeller(c("GLP1" = "GLP-1RA", "DPP4" = "DPP4i", "SGLT2" = "SGLT2i", "TZD" = "TZD", "SU" = "SU"))) +
    theme_bw() +
    guides(fill = guide_legend(nrow = 1, byrow = TRUE)) +
    theme(
      legend.position = "bottom",
      axis.text.x = element_blank(),
      axis.title.x = element_blank(),
      strip.text.x = element_blank()
    ),
  
  # Ethnicity
  discontinuation_per_drug_not_combined %>%
    gather("key", "value", -drugclass, -analysis_var, -group, -variable) %>%
    mutate(
      key = factor(key, levels = c("perc_12m", "perc_6m", "perc_3m"))
    ) %>%
    filter(variable == "ethnicity" & key == "perc_12m") %>%
    ggplot(aes(y = group, x = value, fill = drugclass)) +
    geom_col(position = "dodge") +
    scale_fill_manual(values = c("Pooled" = "black", "SGLT2" = "#E69F00", "GLP1" = "#56B4E9", "SU" = "#CC79A7", "DPP4" = "#0072B2", "TZD" = "#D55E00", "MFN" = "grey"), breaks = c("GLP1", "DPP4", "SGLT2", "TZD", "SU"), labels = c("GLP-1RA", "DPP4i", "SGLT2i", "TZD", "SU"), name = " ", guide = guide_legend(reverse = TRUE)) +
    ggtitle("Ethnicity") +
    scale_x_continuous("Proportion of discontinuation", labels = scales::percent, limits = c(0, 0.5), breaks = seq(0, 0.5, 0.1)) +
    facet_wrap(~drugclass, nrow = 1, labeller = as_labeller(c("GLP1" = "GLP-1RA", "DPP4" = "DPP4i", "SGLT2" = "SGLT2i", "TZD" = "TZD", "SU" = "SU"))) +
    theme_bw() +
    guides(fill = guide_legend(nrow = 1, byrow = TRUE)) +
    theme(
      legend.position = "bottom",
      axis.text.x = element_blank(),
      axis.title.x = element_blank(),
      strip.text.x = element_blank()
    ),
  
  # Quintile
  discontinuation_per_drug_not_combined %>%
    gather("key", "value", -drugclass, -analysis_var, -group, -variable) %>%
    mutate(
      key = factor(key, levels = c("perc_12m", "perc_6m", "perc_3m"))
    ) %>%
    filter(variable == "imd2015" & key == "perc_12m") %>%
    ggplot(aes(y = group, x = value, fill = drugclass)) +
    geom_col(position = "dodge") +
    scale_fill_manual(values = c("Pooled" = "black", "SGLT2" = "#E69F00", "GLP1" = "#56B4E9", "SU" = "#CC79A7", "DPP4" = "#0072B2", "TZD" = "#D55E00", "MFN" = "grey"), breaks = c("GLP1", "DPP4", "SGLT2", "TZD", "SU"), labels = c("GLP-1RA", "DPP4i", "SGLT2i", "TZD", "SU"), name = " ", guide = guide_legend(reverse = TRUE)) +
    scale_y_discrete(limits = rev, labels = c("5 (most deprived)", "4", "3", "2", "1 (least deprived)")) +
    ggtitle("Index of multiple deprivation") +
    scale_x_continuous("Proportion of discontinuation", labels = scales::percent, limits = c(0, 0.5), breaks = seq(0, 0.5, 0.1)) +
    facet_wrap(~drugclass, nrow = 1, labeller = as_labeller(c("GLP1" = "GLP-1RA", "DPP4" = "DPP4i", "SGLT2" = "SGLT2i", "TZD" = "TZD", "SU" = "SU"))) +
    theme_bw() +
    guides(fill = guide_legend(nrow = 1, byrow = TRUE)) +
    theme(
      legend.position = "bottom",
      axis.text.x = element_blank(),
      axis.title.x = element_blank(),
      strip.text.x = element_blank()
    ),
  
  # BMI
  discontinuation_per_drug_not_combined %>%
    gather("key", "value", -drugclass, -analysis_var, -group, -variable) %>%
    mutate(
      key = factor(key, levels = c("perc_12m", "perc_6m", "perc_3m"))
    ) %>%
    filter(variable == "bmi" & key == "perc_12m") %>%
    ggplot(aes(y = group, x = value, fill = drugclass)) +
    geom_col(position = "dodge") +
    scale_fill_manual(values = c("Pooled" = "black", "SGLT2" = "#E69F00", "GLP1" = "#56B4E9", "SU" = "#CC79A7", "DPP4" = "#0072B2", "TZD" = "#D55E00", "MFN" = "grey"), breaks = c("GLP1", "DPP4", "SGLT2", "TZD", "SU"), labels = c("GLP-1RA", "DPP4i", "SGLT2i", "TZD", "SU"), name = " ", guide = guide_legend(reverse = TRUE)) +
    scale_y_discrete(limits = rev, labels = c("35+", "30-34.9", "25-29.9", "<25")) +
    ggtitle("BMI, kg/m2") +
    scale_x_continuous("Proportion of discontinuation", labels = scales::percent, limits = c(0, 0.5), breaks = seq(0, 0.5, 0.1)) +
    facet_wrap(~drugclass, nrow = 1, labeller = as_labeller(c("GLP1" = "GLP-1RA", "DPP4" = "DPP4i", "SGLT2" = "SGLT2i", "TZD" = "TZD", "SU" = "SU"))) +
    theme_bw() +
    guides(fill = guide_legend(nrow = 1, byrow = TRUE)) +
    theme(
      legend.position = "bottom",
      axis.text.x = element_blank(),
      axis.title.x = element_blank(),
      strip.text.x = element_blank()
    ),
  
  # HbA1c
  discontinuation_per_drug_not_combined %>%
    gather("key", "value", -drugclass, -analysis_var, -group, -variable) %>%
    mutate(
      key = factor(key, levels = c("perc_12m", "perc_6m", "perc_3m"))
    ) %>%
    filter(variable == "hba1c" & key == "perc_12m") %>%
    ggplot(aes(y = group, x = value, fill = drugclass)) +
    geom_col(position = "dodge") +
    scale_fill_manual(values = c("Pooled" = "black", "SGLT2" = "#E69F00", "GLP1" = "#56B4E9", "SU" = "#CC79A7", "DPP4" = "#0072B2", "TZD" = "#D55E00", "MFN" = "grey"), breaks = c("GLP1", "DPP4", "SGLT2", "TZD", "SU"), labels = c("GLP-1RA", "DPP4i", "SGLT2i", "TZD", "SU"), name = " ", guide = guide_legend(reverse = TRUE)) +
    scale_y_discrete(limits = rev, labels = c("86+", "75-86", "64-75", "53-64")) +
    ggtitle("HbA1c, mmol/mol") +
    scale_x_continuous("Proportion of discontinuation", labels = scales::percent, limits = c(0, 0.5), breaks = seq(0, 0.5, 0.1)) +
    facet_wrap(~drugclass, nrow = 1, labeller = as_labeller(c("GLP1" = "GLP-1RA", "DPP4" = "DPP4i", "SGLT2" = "SGLT2i", "TZD" = "TZD", "SU" = "SU"))) +
    theme_bw() +
    guides(fill = guide_legend(nrow = 1, byrow = TRUE)) +
    theme(
      legend.position = "bottom",
      axis.text.x = element_blank(),
      axis.title.x = element_blank(),
      strip.text.x = element_blank()
    ),
  
  # eGFR
  discontinuation_per_drug_not_combined %>%
    gather("key", "value", -drugclass, -analysis_var, -group, -variable) %>%
    mutate(
      key = factor(key, levels = c("perc_12m", "perc_6m", "perc_3m"))
    ) %>%
    filter(variable == "preegfr" & key == "perc_12m") %>%
    ggplot(aes(y = group, x = value, fill = drugclass)) +
    geom_col(position = "dodge") +
    scale_fill_manual(values = c("Pooled" = "black", "SGLT2" = "#E69F00", "GLP1" = "#56B4E9", "SU" = "#CC79A7", "DPP4" = "#0072B2", "TZD" = "#D55E00", "MFN" = "grey"), breaks = c("GLP1", "DPP4", "SGLT2", "TZD", "SU"), labels = c("GLP-1RA", "DPP4i", "SGLT2i", "TZD", "SU"), name = " ", guide = guide_legend(reverse = TRUE)) +
    scale_y_discrete(limits = rev, labels = c("90+", "60-90", "<60")) +
    ggtitle("eGFR, ml/min per 1.73 m2") +
    scale_x_continuous("Proportion of discontinuation", labels = scales::percent, limits = c(0, 0.5), breaks = seq(0, 0.5, 0.1)) +
    facet_wrap(~drugclass, nrow = 1, labeller = as_labeller(c("GLP1" = "GLP-1RA", "DPP4" = "DPP4i", "SGLT2" = "SGLT2i", "TZD" = "TZD", "SU" = "SU"))) +
    theme_bw() +
    guides(fill = guide_legend(nrow = 1, byrow = TRUE)) +
    theme(
      legend.position = "bottom",
      axis.text.x = element_blank(),
      axis.title.x = element_blank(),
      strip.text.x = element_blank()
    ),
  
  # Statins
  discontinuation_per_drug_not_combined %>%
    gather("key", "value", -drugclass, -analysis_var, -group, -variable) %>%
    mutate(
      key = factor(key, levels = c("perc_12m", "perc_6m", "perc_3m"))
    ) %>%
    filter(variable == "statins" & key == "perc_12m") %>%
    ggplot(aes(y = group, x = value, fill = drugclass)) +
    geom_col(position = "dodge") +
    scale_fill_manual(values = c("Pooled" = "black", "SGLT2" = "#E69F00", "GLP1" = "#56B4E9", "SU" = "#CC79A7", "DPP4" = "#0072B2", "TZD" = "#D55E00", "MFN" = "grey"), breaks = c("GLP1", "DPP4", "SGLT2", "TZD", "SU"), labels = c("GLP-1RA", "DPP4i", "SGLT2i", "TZD", "SU"), name = " ", guide = guide_legend(reverse = TRUE)) +
    scale_y_discrete(limits = rev, breaks = c(0, 1), labels = c("No", "Yes")) +
    ggtitle("Statins") +
    scale_x_continuous("Proportion of discontinuation", labels = scales::percent, limits = c(0, 0.5), breaks = seq(0, 0.5, 0.1)) +
    facet_wrap(~drugclass, nrow = 1, labeller = as_labeller(c("GLP1" = "GLP-1RA", "DPP4" = "DPP4i", "SGLT2" = "SGLT2i", "TZD" = "TZD", "SU" = "SU"))) +
    theme_bw() +
    guides(fill = guide_legend(nrow = 1, byrow = TRUE)) +
    theme(
      legend.position = "bottom",
      axis.text.x = element_blank(),
      axis.title.x = element_blank(),
      strip.text.x = element_blank()
    ),
  
  # Blood pressure medication
  discontinuation_per_drug_not_combined %>%
    gather("key", "value", -drugclass, -analysis_var, -group, -variable) %>%
    mutate(
      key = factor(key, levels = c("perc_12m", "perc_6m", "perc_3m"))
    ) %>%
    filter(variable == "statins" & key == "perc_12m") %>%
    ggplot(aes(y = group, x = value, fill = drugclass)) +
    geom_col(position = "dodge") +
    scale_fill_manual(values = c("Pooled" = "black", "SGLT2" = "#E69F00", "GLP1" = "#56B4E9", "SU" = "#CC79A7", "DPP4" = "#0072B2", "TZD" = "#D55E00", "MFN" = "grey"), breaks = c("GLP1", "DPP4", "SGLT2", "TZD", "SU"), labels = c("GLP-1RA", "DPP4i", "SGLT2i", "TZD", "SU"), name = " ", guide = guide_legend(reverse = TRUE)) +
    scale_y_discrete(limits = rev, breaks = c(0, 1), labels = c("No", "Yes")) +
    ggtitle("Blood pressure medication") +
    scale_x_continuous("Proportion of discontinuation", labels = scales::percent, limits = c(0, 0.5), breaks = seq(0, 0.5, 0.1)) +
    facet_wrap(~drugclass, nrow = 1, labeller = as_labeller(c("GLP1" = "GLP-1RA", "DPP4" = "DPP4i", "SGLT2" = "SGLT2i", "TZD" = "TZD", "SU" = "SU"))) +
    theme_bw() +
    guides(fill = guide_legend(nrow = 1, byrow = TRUE)) +
    theme(
      legend.position = "bottom",
      axis.text.x = element_blank(),
      axis.title.x = element_blank(),
      strip.text.x = element_blank()
    ),
  
  # Frailty
  discontinuation_per_drug_not_combined %>%
    gather("key", "value", -drugclass, -analysis_var, -group, -variable) %>%
    mutate(
      key = factor(key, levels = c("perc_12m", "perc_6m", "perc_3m"))
    ) %>%
    filter(variable == "frailty" & key == "perc_12m") %>%
    ggplot(aes(y = group, x = value, fill = drugclass)) +
    geom_col(position = "dodge") +
    scale_fill_manual(values = c("Pooled" = "black", "SGLT2" = "#E69F00", "GLP1" = "#56B4E9", "SU" = "#CC79A7", "DPP4" = "#0072B2", "TZD" = "#D55E00", "MFN" = "grey"), breaks = c("GLP1", "DPP4", "SGLT2", "TZD", "SU"), labels = c("GLP-1RA", "DPP4i", "SGLT2i", "TZD", "SU"), name = " ", guide = guide_legend(reverse = TRUE)) +
    scale_y_discrete(limits = rev, breaks = c(0, 1), labels = c("No", "Yes")) +
    ggtitle("Frailty") +
    scale_x_continuous("Proportion of discontinuation", labels = scales::percent, limits = c(0, 0.5), breaks = seq(0, 0.5, 0.1)) +
    facet_wrap(~drugclass, nrow = 1, labeller = as_labeller(c("GLP1" = "GLP-1RA", "DPP4" = "DPP4i", "SGLT2" = "SGLT2i", "TZD" = "TZD", "SU" = "SU"))) +
    theme_bw() +
    guides(fill = guide_legend(nrow = 1, byrow = TRUE)) +
    theme(
      legend.position = "bottom",
      axis.text.x = element_blank(),
      axis.title.x = element_blank(),
      strip.text.x = element_blank()
    ),
  
  # stopping metformin
  discontinuation_per_drug_not_combined %>%
    gather("key", "value", -drugclass, -analysis_var, -group, -variable) %>%
    mutate(
      key = factor(key, levels = c("perc_12m", "perc_6m", "perc_3m"))
    ) %>%
    filter(variable == "stopdrug_MFN" & key == "perc_12m") %>%
    ggplot(aes(y = group, x = value, fill = drugclass)) +
    geom_col() +
    scale_fill_manual(values = c("Pooled" = "black", "SGLT2" = "#E69F00", "GLP1" = "#56B4E9", "SU" = "#CC79A7", "DPP4" = "#0072B2", "TZD" = "#D55E00", "MFN" = "grey"), breaks = c("GLP1", "DPP4", "SGLT2", "TZD", "SU"), labels = c("GLP-1RA", "DPP4i", "SGLT2i", "TZD", "SU"), name = " ", guide = guide_legend(reverse = TRUE)) +
    scale_y_discrete(limits = rev, breaks = c(0, 1), labels = c("No", "Yes")) +
    ggtitle("History of Metformin discontinuation within 3-months") +
    scale_x_continuous("Proportion of discontinuation", labels = scales::percent, limits = c(0, 0.5), breaks = seq(0, 0.5, 0.1)) +
    facet_wrap(~drugclass, nrow = 1, labeller = as_labeller(c("GLP1" = "GLP-1RA", "DPP4" = "DPP4i", "SGLT2" = "SGLT2i", "TZD" = "TZD", "SU" = "SU"))) +
    theme_bw() +
    guides(fill = guide_legend(nrow = 1, byrow = TRUE)) +
    theme(
      legend.position = "bottom",
      axis.text.x = element_blank(),
      axis.title.x = element_blank(),
      strip.text.x = element_blank()
    ),
  
  # numdrugs
  discontinuation_per_drug_not_combined %>%
    gather("key", "value", -drugclass, -analysis_var, -group, -variable) %>%
    mutate(
      key = factor(key, levels = c("perc_12m", "perc_6m", "perc_3m"))
    ) %>%
    filter(variable == "numdrugs" & key == "perc_12m") %>%
    ggplot(aes(y = group, x = value, fill = drugclass)) +
    geom_col() +
    scale_fill_manual(values = c("Pooled" = "black", "SGLT2" = "#E69F00", "GLP1" = "#56B4E9", "SU" = "#CC79A7", "DPP4" = "#0072B2", "TZD" = "#D55E00", "MFN" = "grey"), breaks = c("GLP1", "DPP4", "SGLT2", "TZD", "SU"), labels = c("GLP-1RA", "DPP4i", "SGLT2i", "TZD", "SU"), name = " ", guide = guide_legend(reverse = TRUE)) +
    scale_y_discrete(limits = rev, labels = c("2+", "1", "0")) +
    ggtitle("Number of other current glucose-lowering therapies") +
    scale_x_continuous("Proportion of discontinuation", labels = scales::percent, limits = c(0, 0.5), breaks = seq(0, 0.5, 0.1)) +
    facet_wrap(~drugclass, nrow = 1, labeller = as_labeller(c("GLP1" = "GLP-1RA", "DPP4" = "DPP4i", "SGLT2" = "SGLT2i", "TZD" = "TZD", "SU" = "SU"))) +
    theme_bw() +
    guides(fill = guide_legend(nrow = 1, byrow = TRUE)) +
    theme(
      legend.position = "bottom",
      axis.text.x = element_blank(),
      axis.title.x = element_blank(),
      strip.text.x = element_blank()
    ),
  
  # drugline
  discontinuation_per_drug_not_combined %>%
    gather("key", "value", -drugclass, -analysis_var, -group, -variable) %>%
    mutate(
      key = factor(key, levels = c("perc_12m", "perc_6m", "perc_3m"))
    ) %>%
    filter(variable == "drugline" & key == "perc_12m") %>%
    ggplot(aes(y = group, x = value, fill = drugclass)) +
    geom_col(position = "dodge") +
    scale_fill_manual(values = c("Pooled" = "black", "SGLT2" = "#E69F00", "GLP1" = "#56B4E9", "SU" = "#CC79A7", "DPP4" = "#0072B2", "TZD" = "#D55E00", "MFN" = "grey"), breaks = c("GLP1", "DPP4", "SGLT2", "TZD", "SU"), labels = c("GLP-1RA", "DPP4i", "SGLT2i", "TZD", "SU"), name = " ", guide = guide_legend(reverse = TRUE)) +
    scale_y_discrete(limits = rev) +
    ggtitle("Number of therapy classes ever prescribed") +
    scale_x_continuous("Proportion of discontinuation", labels = scales::percent, limits = c(0, 0.5), breaks = seq(0, 0.5, 0.1)) +
    facet_wrap(~drugclass, nrow = 1, labeller = as_labeller(c("GLP1" = "GLP-1RA", "DPP4" = "DPP4i", "SGLT2" = "SGLT2i", "TZD" = "TZD", "SU" = "SU"))) +
    theme_bw() +
    guides(fill = guide_legend(nrow = 1, byrow = TRUE)) +
    theme(
      legend.position = "bottom",
      axis.text.x = element_text(size = 10),
      axis.title.x = element_text(size = 11),
      strip.text.x = element_blank()
    )
  
  , ncol = 1
  
) + 
  plot_layout(
    guides = "collect",
    height = c(4, 5, 2, 3, 4, 5, 4, 4, 3, 2, 2, 2, 2, 3, 3),
    # axis_titles = "collect" # not working, manually done in the code
  ) +
  plot_layout(
    axis_titles = "collect"
  ) &
  theme(
    legend.position = "none",
    plot.title.position = "plot",
    axis.title.y = element_blank(),
    axis.text.y = element_text(size = 10),
    panel.spacing = unit(0.3, "cm", data = NULL)
  )


dev.off()


###############################################################################
###############################################################################
############################### Read Data In ##################################
###############################################################################
###############################################################################

# load dataset
cprd_dataset <- set_up_data(
  raw_data = "20240814_t2d_1stinstance",
  diagnosis = FALSE,
  therapies = c("DPP4", "GLP1", "SGLT2", "SU", "TZD"),
  dataset = "3m.disc.dataset",
  follow_up = "6-months",
  full_prescribing_history = TRUE
) %>%
  drop_na(-prehdl, -stopdrug_6m_6mFU, -stopdrug_12m_6mFU)


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


drugs = c("DPP4", "GLP1", "SGLT2", "SU", "TZD")


# imd2015_10
cprd_dataset <- cprd_dataset %>%
  mutate(imd2015_10 = as.numeric(imd2015_10))

imd2015_sum <- univariate_analysis(
  cprd_dataset, drugs = drugs, outcome = "stopdrug_3m_6mFU", 
  variable = "imd2015_10", variable_name = "Index of multiple deprivation", type = "Clinical features & biomarkers"
)


# drugline
cprd_dataset <- cprd_dataset %>%
  mutate(drugline = as.numeric(drugline))

drugline_sum <- univariate_analysis(
  cprd_dataset, drugs = drugs, outcome = "stopdrug_3m_6mFU", 
  variable = "drugline", variable_name = "Number of therapy classes ever prescribed", type = "Clinical features & biomarkers"
)

# numdrugs
cprd_dataset <- cprd_dataset %>%
  mutate(numdrugs = as.numeric(numdrugs))

numdrugs_sum <- univariate_analysis(
  cprd_dataset, drugs = drugs, outcome = "stopdrug_3m_6mFU", 
  variable = "numdrugs", variable_name = "Number of other current glucose-lowering therapies", type = "Clinical features & biomarkers"
)

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


# prehdl

cprd <- cprd_dataset %>%
  mutate_at(c("prehdl"), ~(scale(.) %>% as.vector))

prehdl_sum <- univariate_analysis(
  cprd_dataset, drugs = drugs, outcome = "stopdrug_3m_6mFU", 
  variable = "prehdl", variable_name = "HDL (per SD)", type = "Clinical features & biomarkers"
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
  variable = "stopdrug_3m_3mFU_MFN_hist", variable_name = "History of Metformin discontinuation within 3-months", type = "Behavioural"
)


# ethnicity

cprd_dataset_white_southasian <- cprd_dataset %>%
  filter(ethnicity_5cat %in% c(0, 1)) %>%
  mutate(ethnicity_5cat = factor(ethnicity_5cat, levels = c(0, 1), labels = c("White", "South Asian")))

ethnicity_5cat_sum_white_southasian <- univariate_analysis(
  cprd_dataset_white_southasian, drugs = drugs, outcome = "stopdrug_3m_6mFU", 
  variable = "ethnicity_5cat", variable_name = "Ethnicity: South Asian (ref White)", type = "Behavioural"
)

cprd_dataset_white_black <- cprd_dataset %>%
  filter(ethnicity_5cat %in% c(0, 2)) %>%
  mutate(ethnicity_5cat = factor(ethnicity_5cat, levels = c(0, 2), labels = c("White", "Black")))

ethnicity_5cat_sum_white_black <- univariate_analysis(
  cprd_dataset_white_black, drugs = drugs, outcome = "stopdrug_3m_6mFU", 
  variable = "ethnicity_5cat", variable_name = "Ethnicity: Black (ref White)", type = "Behavioural"
)

cprd_dataset_white_other <- cprd_dataset %>%
  filter(ethnicity_5cat %in% c(0, 3)) %>%
  mutate(ethnicity_5cat = factor(ethnicity_5cat, levels = c(0, 3), labels = c("White", "Other")))

ethnicity_5cat_sum_white_other <- univariate_analysis(
  cprd_dataset_white_other, drugs = drugs, outcome = "stopdrug_3m_6mFU", 
  variable = "ethnicity_5cat", variable_name = "Ethnicity: Other (ref White)", type = "Behavioural"
)

cprd_dataset_white_mixed <- cprd_dataset %>%
  filter(ethnicity_5cat %in% c(0, 4)) %>%
  mutate(ethnicity_5cat = factor(ethnicity_5cat, levels = c(0, 4), labels = c("White", "Mixed")))

ethnicity_5cat_sum_white_mixed <- univariate_analysis(
  cprd_dataset_white_mixed, drugs = drugs, outcome = "stopdrug_3m_6mFU", 
  variable = "ethnicity_5cat", variable_name = "Ethnicity: Mixed (ref White)", type = "Behavioural"
)

cprd_dataset_white_othermixed <- cprd_dataset %>%
  filter(ethnicity_5cat %in% c(0, 3, 4)) %>%
  mutate(ethnicity_5cat = factor(ethnicity_5cat, levels = c(0, 3, 4), labels = c("White", "Mixed & Other", "Mixed & Other")))

ethnicity_5cat_sum_white_othermixed <- univariate_analysis(
  cprd_dataset_white_othermixed, drugs = drugs, outcome = "stopdrug_3m_6mFU", 
  variable = "ethnicity_5cat", variable_name = "Ethnicity: Mixed & Other (ref White)", type = "Behavioural"
)

cprd_dataset_test <- cprd_dataset %>%
  mutate(ethnicity_5cat = factor(ethnicity_5cat, levels = c(0, 1, 2, 3, 4), labels = c("White", "Other", "Other" , "Other" , "Other")))

ethnicity_5cat_sum <- univariate_analysis(
  cprd_dataset_test, drugs = drugs, outcome = "stopdrug_3m_6mFU",
  variable = "ethnicity_5cat", variable_name = "Ethnicity (ref White)", type = "Behavioural"
)

# smoking_cat

cprd_dataset <- cprd_dataset %>%
  mutate(smoking_cat = factor(smoking_cat, levels = c("Active smoker", "Ex-smoker", "Non-smoker"), labels = c("Active smoker", "Other", "Other")))

smoking_cat_sum <- univariate_analysis(
  cprd_dataset, drugs = drugs, outcome = "stopdrug_3m_6mFU", 
  variable = "smoking_cat", variable_name = "Smoking (ref Active Smoker)", type = "Behavioural"
)


# preckdstage

cprd_dataset <- cprd_dataset %>%
  mutate(preckdstage = factor(preckdstage, levels = c("stage_0", "stage_1", "stage_2", "stage_3a", "stage_3b", "stage_4", "stage_5"), labels = c("stage 0/1/2", "stage 0/1/2", "stage 0/1/2", "stage 3/4/5", "stage 3/4/5", "stage 3/4/5", "stage 3/4/5")))

preckdstage_sum <- univariate_analysis(
  cprd_dataset, drugs = drugs, outcome = "stopdrug_3m_6mFU", 
  variable = "preckdstage", variable_name = "CKD stage 3/4/5 (ref stage 0/1/2)", type = "Comorbidity event"
)


# predrug_frailty_proxy

predrug_frailty_proxy_sum <- univariate_analysis(
  cprd_dataset, drugs = drugs, outcome = "stopdrug_3m_6mFU", 
  variable = "predrug_frailty_proxy", variable_name = "Frailty", type = "Comorbidity event"
)



coefficients_3m <- rbind(
  imd2015_sum,
  drugline_sum,
  dstartdate_dm_dur_sum,
  numdrugs_sum,
  dstartdate_age_sum,
  prehba1c_sum,
  preegfr_sum,
  prebmi_sum,
  gender_sum,
  predrug_bloodmed_sum,
  predrug_statins_sum,
  stopdrug_3m_3mFU_MFN_hist_sum,
  ethnicity_5cat_sum_white_southasian,
  ethnicity_5cat_sum_white_black,
  ethnicity_5cat_sum_white_othermixed,
  ethnicity_5cat_sum,
  smoking_cat_sum,
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

saveRDS(coefficients_3m, "results/Models/Predictions/univariate_coefficients_3m.rds")


pdf("results/figures/03.univariate_analysis_3m_variable_layout.pdf", width = 11, heigh = 12)

coefficients_3m %>%
  filter(variable %in% c("Age (per SD)", "Diabetes duration (per SD)", "Sex (ref Male)", "Smoking (ref Active Smoker)", "Ethnicity: South Asian (ref White)", "Ethnicity: Black (ref White)", "Ethnicity: Mixed & Other (ref White)", "Index of multiple deprivation", "BMI (per SD)", "HbA1c (per SD)", "eGFR (per SD)", "Statins", "Blood pressure medication", "Frailty", "History of Metformin discontinuation within 3-months", "Number of other current glucose-lowering therapies", "Number of therapy classes ever prescribed")) %>%
  mutate(
    variable = factor(variable, levels = c("Age (per SD)", "Diabetes duration (per SD)", "Sex (ref Male)", "Smoking (ref Active Smoker)", "Ethnicity: South Asian (ref White)", "Ethnicity: Black (ref White)", "Ethnicity: Mixed & Other (ref White)", "Index of multiple deprivation", "BMI (per SD)", "HbA1c (per SD)", "eGFR (per SD)", "Statins", "Blood pressure medication", "Frailty", "History of Metformin discontinuation within 3-months", "Number of other current glucose-lowering therapies", "Number of therapy classes ever prescribed")),
    drug = factor(drug, levels = rev(c("Pooled", "GLP1", "DPP4", "SGLT2", "TZD", "SU")))
  ) %>%
  ggplot() + 
  geom_vline(aes(xintercept = 1), colour = "black") +
  geom_errorbar(aes(x = coef, xmin = LCI, xmax = UCI, y = variable, colour = drug), position = position_dodge(width = 1)) +
  geom_pointrange(aes(x = coef, xmin = LCI, xmax = UCI, y = variable, colour = drug), position = position_dodge(width = 1)) +
  xlab("Odds Ratio") +
  scale_x_log10(breaks = c(0.25, 0.5, 1, 2, 4), limits = c(0.4, 4.1)) +
  facet_grid(variable~., scales = "free_y") +
  scale_colour_manual(values = c("Pooled" = "black", "SGLT2" = "#E69F00", "GLP1" = "#56B4E9", "SU" = "#CC79A7", "DPP4" = "#0072B2", "TZD" = "#D55E00", "MFN" = "grey"), 
                      breaks = rev(c("Pooled", "GLP1", "DPP4", "SGLT2", "TZD", "SU")), labels = rev(c("Pooled", "GLP-1RA", "DPP4i", "SGLT2i", "TZD", "SU")),
                      name = "Therapy", guide = guide_legend(reverse = TRUE)) +
  theme_bw() +
  theme(
    axis.title.y = element_blank(),
    legend.position = "right",
    strip.text = element_blank(),
    legend.direction = "vertical",
    axis.text.y = element_text(size = 15),
    axis.text.x = element_text(size = 15),
    axis.title.x = element_text(size = 18)
  )

dev.off()




pdf("results/figures/03.univariate_analysis_3m.pdf", width = 12, height = 8)

patchwork::wrap_plots(
  
  coefficients_3m %>%
    filter(type == "Behavioural") %>%
    filter(variable %in% c("Blood pressure medication", "Statins", "History of Metformin discontinuation within 3-months", "Ethnicity (ref White)", "Smoking (ref Active Smoker)")) %>%
    ggplot() + 
    geom_vline(aes(xintercept = 1), colour = "black") +
    geom_errorbar(aes(x = coef, xmin = LCI, xmax = UCI, y = variable, colour = drug),position = position_dodge(width = 1)) +
    geom_pointrange(aes(x = coef, xmin = LCI, xmax = UCI, y = variable, colour = drug), position = position_dodge(width = 1)) +
    xlab("Odds Ratio") +
    scale_x_log10(breaks = c(0.25, 0.5, 1, 2, 4), limits = c(0.4, 4.1)) +
    facet_grid(variable~., scales = "free_y") +
    scale_colour_manual(values = c("Pooled" = "black", "SGLT2" = "#E69F00", "GLP1" = "#56B4E9", "SU" = "#CC79A7", "DPP4" = "#0072B2", "TZD" = "#D55E00", "MFN" = "grey"), 
                        breaks = rev(c("Pooled", "GLP1", "DPP4", "SGLT2", "TZD", "SU")), labels = rev(c("Pooled", "GLP-1RA", "DPP4i", "SGLT2i", "TZD", "SU")),
                        name = "Therapy", guide = guide_legend(reverse = TRUE)) +
    theme_bw() +
    theme(
      axis.title.y = element_blank(),
      legend.position = "bottom",
      strip.text = element_blank()
    ),
  
  
  coefficients_3m %>%
    filter(type == "Clinical features & biomarkers") %>%
    filter(variable %in% c("Diabetes duration (per SD)", "Age (per SD)", "HbA1c (per SD)", "BMI (per SD)", "Sex (ref Male)")) %>%
    mutate(variable = factor(variable, levels = c("Age (per SD)", "Diabetes duration (per SD)", "BMI (per SD)", "HbA1c (per SD)", "Sex (ref Male)"))) %>%
    ggplot() + 
    geom_vline(aes(xintercept = 1), colour = "black") +
    geom_errorbar(aes(x = coef, xmin = LCI, xmax = UCI, y = variable, colour = drug),position = position_dodge(width = 1)) +
    geom_pointrange(aes(x = coef, xmin = LCI, xmax = UCI, y = variable, colour = drug), position = position_dodge(width = 1)) +
    xlab("Odds Ratio") +
    scale_x_log10(breaks = c(0.25, 0.5, 1, 2, 4), limits = c(0.4, 4.1)) +
    facet_grid(variable~., scales = "free_y") +
    scale_colour_manual(values = c("Pooled" = "black", "SGLT2" = "#E69F00", "GLP1" = "#56B4E9", "SU" = "#CC79A7", "DPP4" = "#0072B2", "TZD" = "#D55E00", "MFN" = "grey"), 
                        breaks = rev(c("Pooled", "GLP1", "DPP4", "SGLT2", "TZD", "SU")), labels = rev(c("Pooled", "GLP-1RA", "DPP4i", "SGLT2i", "TZD", "SU")),
                        name = "Therapy", guide = guide_legend(reverse = TRUE)) +
    theme_bw() +
    theme(
      axis.title.y = element_blank(),
      legend.position = "bottom",
      strip.text = element_blank()
    ),
  
  coefficients_3m %>%
    filter(type == "Comorbidity event") %>%
    filter(variable %in% c("Frailty")) %>%
    ggplot() + 
    geom_vline(aes(xintercept = 1), colour = "black") +
    geom_errorbar(aes(x = coef, xmin = LCI, xmax = UCI, y = variable, colour = drug),position = position_dodge(width = 1)) +
    geom_pointrange(aes(x = coef, xmin = LCI, xmax = UCI, y = variable, colour = drug), position = position_dodge(width = 1)) +
    xlab("Odds Ratio") +
    scale_x_log10(breaks = c(0.25, 0.5, 1, 2, 4), limits = c(0.4, 4.1)) +
    facet_grid(variable~., scales = "free_y") +
    scale_colour_manual(values = c("Pooled" = "black", "SGLT2" = "#E69F00", "GLP1" = "#56B4E9", "SU" = "#CC79A7", "DPP4" = "#0072B2", "TZD" = "#D55E00", "MFN" = "grey"), 
                        breaks = rev(c("Pooled", "GLP1", "DPP4", "SGLT2", "TZD", "SU")), labels = rev(c("Pooled", "GLP-1RA", "DPP4i", "SGLT2i", "TZD", "SU")),
                        name = "Therapy", guide = guide_legend(reverse = TRUE)) +
    theme_bw() +
    theme(
      axis.title.y = element_blank(),
      legend.position = "bottom",
      strip.text = element_blank()
    ),
  
  ncol = 1, nrow = 3
  
) + 
  plot_layout(guides = "collect", height = c(5, 5, 1), axis_titles = "collect", axes = "collect") +
  plot_annotation(tag_levels = list(c("Behavioural", "Clinical features & biomarkers", "Comorbidity event"))) & 
  theme(
    legend.direction = "vertical",
    axis.text.y = element_text(size = 15),
    axis.text.x = element_text(size = 15),
    axis.title.x = element_text(size = 18)
  )

dev.off()



## Table of coefficients
coefficients_3m %>%
  mutate(numbers = paste0(format(round(coef, 2), nsmall = 2), " (", format(round(LCI, 2), nsmall = 2), ",", format(round(UCI, 2), nsmall = 2),")")) %>%
  select(variable, drug, numbers) %>%
  pivot_wider(names_from = drug, values_from = numbers) %>%
  write.csv("results/tables/odds_ratios_3m.csv", na="") 


###############################################################################
###############################################################################
################# Univariate analysis - 6m discontinuation ####################
###############################################################################
###############################################################################


# load dataset
cprd_dataset <- cprd_dataset %>%
  drop_na(stopdrug_6m_6mFU)


drugs = c("DPP4", "GLP1", "SGLT2", "SU", "TZD")



# imd2015_10
cprd_dataset <- cprd_dataset %>%
  mutate(imd2015_10 = as.numeric(imd2015_10))

imd2015_sum <- univariate_analysis(
  cprd_dataset, drugs = drugs, outcome = "stopdrug_6m_6mFU", 
  variable = "imd2015_10", variable_name = "Index of multiple deprivation", type = "Clinical features & biomarkers"
)


# drugline
cprd_dataset <- cprd_dataset %>%
  mutate(drugline = as.numeric(drugline))

drugline_sum <- univariate_analysis(
  cprd_dataset, drugs = drugs, outcome = "stopdrug_6m_6mFU", 
  variable = "drugline", variable_name = "Number of therapy classes ever prescribed", type = "Clinical features & biomarkers"
)

# numdrugs
cprd_dataset <- cprd_dataset %>%
  mutate(numdrugs = as.numeric(numdrugs))

numdrugs_sum <- univariate_analysis(
  cprd_dataset, drugs = drugs, outcome = "stopdrug_3m_6mFU", 
  variable = "numdrugs", variable_name = "Number of other current glucose-lowering therapies", type = "Clinical features & biomarkers"
)



# dstartdate_dm_dur

dstartdate_dm_dur_sum <- univariate_analysis(
  cprd_dataset, drugs = drugs, outcome = "stopdrug_6m_6mFU", 
  variable = "dstartdate_dm_dur", variable_name = "Diabetes duration (per SD)", type = "Clinical features & biomarkers"
)


# dstartdate_age

dstartdate_age_sum <- univariate_analysis(
  cprd_dataset, drugs = drugs, outcome = "stopdrug_6m_6mFU", 
  variable = "dstartdate_age", variable_name = "Age (per SD)", type = "Clinical features & biomarkers"
)


# prehba1c

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

prebmi_sum <- univariate_analysis(
  cprd_dataset, drugs = drugs, outcome = "stopdrug_6m_6mFU", 
  variable = "prebmi", variable_name = "BMI (per SD)", type = "Clinical features & biomarkers"
)

# prehdl

prehdl_sum <- univariate_analysis(
  cprd_dataset, drugs = drugs, outcome = "stopdrug_6m_6mFU", 
  variable = "prehdl", variable_name = "HDL (per SD)", type = "Clinical features & biomarkers"
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
  variable = "stopdrug_3m_3mFU_MFN_hist", variable_name = "History of Metformin discontinuation within 3-months", type = "Behavioural"
)


# ethnicity

ethnicity_5cat_sum <- univariate_analysis(
  cprd_dataset, drugs = drugs, outcome = "stopdrug_6m_6mFU", 
  variable = "ethnicity_5cat", variable_name = "Ethnicity (ref White)", type = "Behavioural"
)


# smoking_cat

smoking_cat_sum <- univariate_analysis(
  cprd_dataset, drugs = drugs, outcome = "stopdrug_6m_6mFU", 
  variable = "smoking_cat", variable_name = "Smoking (ref Active Smoker)", type = "Behavioural"
)


# preckdstage

preckdstage_sum <- univariate_analysis(
  cprd_dataset, drugs = drugs, outcome = "stopdrug_6m_6mFU", 
  variable = "preckdstage", variable_name = "CKD stage 3/4/5 (ref stage 0/1/2)", type = "Comorbidity event"
)


# predrug_frailty_proxy

predrug_frailty_proxy_sum <- univariate_analysis(
  cprd_dataset, drugs = drugs, outcome = "stopdrug_6m_6mFU", 
  variable = "predrug_frailty_proxy", variable_name = "Frailty", type = "Comorbidity event"
)



coefficients_6m <- rbind(
  imd2015_sum,
  drugline_sum,
  numdrugs_sum,
  dstartdate_dm_dur_sum,
  dstartdate_age_sum,
  prehba1c_sum,
  preegfr_sum,
  prebmi_sum,
  gender_sum,
  predrug_bloodmed_sum,
  predrug_statins_sum,
  stopdrug_3m_3mFU_MFN_hist_sum,
  ethnicity_5cat_sum,
  smoking_cat_sum,
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


saveRDS(coefficients_6m, "results/Models/Predictions/univariate_coefficients_6m.rds")



pdf("results/figures/03.univariate_analysis_6m.pdf", width = 12, height = 8)

wrap_plots(
  
  coefficients_6m %>%
    filter(type == "Behavioural") %>%
    filter(variable %in% c("Blood pressure medication", "Statins", "History of Metformin discontinuation within 3-months", "Ethnicity (ref White)", "Smoking (ref Active Smoker)")) %>%
    ggplot() + 
    geom_vline(aes(xintercept = 1), colour = "black") +
    geom_errorbar(aes(x = coef, xmin = LCI, xmax = UCI, y = variable, colour = drug),position = position_dodge(width = 1)) +
    geom_pointrange(aes(x = coef, xmin = LCI, xmax = UCI, y = variable, colour = drug), position = position_dodge(width = 1)) +
    xlab("Odds Ratio") +
    scale_x_log10(breaks = c(0.25, 0.5, 1, 2, 4), limits = c(0.4, 4.1)) +
    facet_grid(variable~., scales = "free_y") +
    scale_colour_manual(values = c("Pooled" = "black", "SGLT2" = "#E69F00", "GLP1" = "#56B4E9", "SU" = "#CC79A7", "DPP4" = "#0072B2", "TZD" = "#D55E00", "MFN" = "grey"), 
                        breaks = rev(c("Pooled", "GLP1", "DPP4", "SGLT2", "TZD", "SU")), labels = rev(c("Pooled", "GLP-1RA", "DPP4i", "SGLT2i", "TZD", "SU")),
                        name = "Therapy", guide = guide_legend(reverse = TRUE)) +
    theme_bw() +
    theme(
      axis.title.y = element_blank(),
      legend.position = "bottom",
      strip.text = element_blank()
    ),
  
  
  coefficients_6m %>%
    filter(type == "Clinical features & biomarkers") %>%
    filter(variable %in% c("Diabetes duration (per SD)", "Age (per SD)", "HbA1c (per SD)", "BMI (per SD)", "Sex (ref Male)")) %>%
    mutate(variable = factor(variable, levels = c("Age (per SD)", "Diabetes duration (per SD)", "BMI (per SD)", "HbA1c (per SD)", "Sex (ref Male)"))) %>%
    ggplot() + 
    geom_vline(aes(xintercept = 1), colour = "black") +
    geom_errorbar(aes(x = coef, xmin = LCI, xmax = UCI, y = variable, colour = drug),position = position_dodge(width = 1)) +
    geom_pointrange(aes(x = coef, xmin = LCI, xmax = UCI, y = variable, colour = drug), position = position_dodge(width = 1)) +
    xlab("Odds Ratio") +
    scale_x_log10(breaks = c(0.25, 0.5, 1, 2, 4), limits = c(0.4, 4.1)) +
    facet_grid(variable~., scales = "free_y") +
    scale_colour_manual(values = c("Pooled" = "black", "SGLT2" = "#E69F00", "GLP1" = "#56B4E9", "SU" = "#CC79A7", "DPP4" = "#0072B2", "TZD" = "#D55E00", "MFN" = "grey"), 
                        breaks = rev(c("Pooled", "GLP1", "DPP4", "SGLT2", "TZD", "SU")), labels = rev(c("Pooled", "GLP-1RA", "DPP4i", "SGLT2i", "TZD", "SU")),
                        name = "Therapy", guide = guide_legend(reverse = TRUE)) +
    theme_bw() +
    theme(
      axis.title.y = element_blank(),
      legend.position = "bottom",
      strip.text = element_blank()
    ),
  
  coefficients_6m %>%
    filter(type == "Comorbidity event") %>%
    filter(variable %in% c("Frailty")) %>%
    ggplot() + 
    geom_vline(aes(xintercept = 1), colour = "black") +
    geom_errorbar(aes(x = coef, xmin = LCI, xmax = UCI, y = variable, colour = drug),position = position_dodge(width = 1)) +
    geom_pointrange(aes(x = coef, xmin = LCI, xmax = UCI, y = variable, colour = drug), position = position_dodge(width = 1)) +
    xlab("Odds Ratio") +
    scale_x_log10(breaks = c(0.25, 0.5, 1, 2, 4), limits = c(0.4, 4.1)) +
    facet_grid(variable~., scales = "free_y") +
    scale_colour_manual(values = c("Pooled" = "black", "SGLT2" = "#E69F00", "GLP1" = "#56B4E9", "SU" = "#CC79A7", "DPP4" = "#0072B2", "TZD" = "#D55E00", "MFN" = "grey"), 
                        breaks = rev(c("Pooled", "GLP1", "DPP4", "SGLT2", "TZD", "SU")), labels = rev(c("Pooled", "GLP-1RA", "DPP4i", "SGLT2i", "TZD", "SU")),
                        name = "Therapy", guide = guide_legend(reverse = TRUE)) +
    theme_bw() +
    theme(
      axis.title.y = element_blank(),
      legend.position = "bottom",
      strip.text = element_blank()
    ),
  
  ncol = 1, nrow = 3
  
) + 
  plot_layout(guides = "collect", height = c(5, 5, 1), axis_titles = "collect", axes = "collect") +
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


drugs = c("DPP4", "GLP1", "SGLT2", "SU", "TZD")


# imd2015_10
cprd_dataset <- cprd_dataset %>%
  mutate(imd2015_10 = as.numeric(imd2015_10))

imd2015_sum <- univariate_analysis(
  cprd_dataset, drugs = drugs, outcome = "stopdrug_12m_6mFU", 
  variable = "imd2015_10", variable_name = "Index of multiple deprivation", type = "Clinical features & biomarkers"
)


# drugline
cprd_dataset <- cprd_dataset %>%
  mutate(drugline = as.numeric(drugline))

drugline_sum <- univariate_analysis(
  cprd_dataset, drugs = drugs, outcome = "stopdrug_12m_6mFU", 
  variable = "drugline", variable_name = "Number of therapy classes ever prescribed", type = "Clinical features & biomarkers"
)


# numdrugs
cprd_dataset <- cprd_dataset %>%
  mutate(numdrugs = as.numeric(numdrugs))

numdrugs_sum <- univariate_analysis(
  cprd_dataset, drugs = drugs, outcome = "stopdrug_3m_6mFU", 
  variable = "numdrugs", variable_name = "Number of other current glucose-lowering therapies", type = "Clinical features & biomarkers"
)



# dstartdate_dm_dur

dstartdate_dm_dur_sum <- univariate_analysis(
  cprd_dataset, drugs = drugs, outcome = "stopdrug_12m_6mFU", 
  variable = "dstartdate_dm_dur", variable_name = "Diabetes duration (per SD)", type = "Clinical features & biomarkers"
)


# dstartdate_age

dstartdate_age_sum <- univariate_analysis(
  cprd_dataset, drugs = drugs, outcome = "stopdrug_12m_6mFU", 
  variable = "dstartdate_age", variable_name = "Age (per SD)", type = "Clinical features & biomarkers"
)


# prehba1c

prehba1c_sum <- univariate_analysis(
  cprd_dataset, drugs = drugs, outcome = "stopdrug_12m_6mFU", 
  variable = "prehba1c", variable_name = "HbA1c (per SD)", type = "Clinical features & biomarkers"
)


# preegfr

preegfr_sum <- univariate_analysis(
  cprd_dataset, drugs = drugs, outcome = "stopdrug_12m_6mFU", 
  variable = "preegfr", variable_name = "eGFR (per SD)", type = "Clinical features & biomarkers"
)


# prebmi

prebmi_sum <- univariate_analysis(
  cprd_dataset, drugs = drugs, outcome = "stopdrug_12m_6mFU", 
  variable = "prebmi", variable_name = "BMI (per SD)", type = "Clinical features & biomarkers"
)

# prehdl

prehdl_sum <- univariate_analysis(
  cprd_dataset, drugs = drugs, outcome = "stopdrug_12m_6mFU", 
  variable = "prehdl", variable_name = "HDL (per SD)", type = "Clinical features & biomarkers"
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
  variable = "stopdrug_3m_3mFU_MFN_hist", variable_name = "History of Metformin discontinuation within 3-months", type = "Behavioural"
)


# ethnicity

ethnicity_5cat_sum <- univariate_analysis(
  cprd_dataset, drugs = drugs, outcome = "stopdrug_3m_6mFU", 
  variable = "ethnicity_5cat", variable_name = "Ethnicity (ref White)", type = "Behavioural"
)


# smoking_cat

smoking_cat_sum <- univariate_analysis(
  cprd_dataset, drugs = drugs, outcome = "stopdrug_12m_6mFU", 
  variable = "smoking_cat", variable_name = "Smoking (ref Active Smoker)", type = "Behavioural"
)


# preckdstage

preckdstage_sum <- univariate_analysis(
  cprd_dataset, drugs = drugs, outcome = "stopdrug_12m_6mFU", 
  variable = "preckdstage", variable_name = "CKD stage 3/4/5 (ref stage 0/1/2)", type = "Comorbidity event"
)


# predrug_frailty_proxy

predrug_frailty_proxy_sum <- univariate_analysis(
  cprd_dataset, drugs = drugs, outcome = "stopdrug_12m_6mFU", 
  variable = "predrug_frailty_proxy", variable_name = "Frailty", type = "Comorbidity event"
)



coefficients_12m <- rbind(
  imd2015_sum,
  drugline_sum,
  numdrugs_sum,
  dstartdate_dm_dur_sum,
  dstartdate_age_sum,
  prehba1c_sum,
  preegfr_sum,
  prebmi_sum,
  gender_sum,
  predrug_bloodmed_sum,
  predrug_statins_sum,
  stopdrug_3m_3mFU_MFN_hist_sum,
  ethnicity_5cat_sum,
  smoking_cat_sum,
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


saveRDS(coefficients_12m, "results/Models/Predictions/univariate_coefficients_12m.rds")



pdf("results/figures/03.univariate_analysis_12m.pdf", width = 12, height = 8)



wrap_plots(
  
  coefficients_12m %>%
    filter(type == "Behavioural") %>%
    filter(variable %in% c("Blood pressure medication", "Statins", "History of Metformin discontinuation within 3-months", "Ethnicity (ref White)", "Smoking (ref Active Smoker)")) %>%
    ggplot() + 
    geom_vline(aes(xintercept = 1), colour = "black") +
    geom_errorbar(aes(x = coef, xmin = LCI, xmax = UCI, y = variable, colour = drug),position = position_dodge(width = 1)) +
    geom_pointrange(aes(x = coef, xmin = LCI, xmax = UCI, y = variable, colour = drug), position = position_dodge(width = 1)) +
    xlab("Odds Ratio") +
    scale_x_log10(breaks = c(0.25, 0.5, 1, 2, 4), limits = c(0.4, 4.1)) +
    facet_grid(variable~., scales = "free_y") +
    scale_colour_manual(values = c("Pooled" = "black", "SGLT2" = "#E69F00", "GLP1" = "#56B4E9", "SU" = "#CC79A7", "DPP4" = "#0072B2", "TZD" = "#D55E00", "MFN" = "grey"), 
                        breaks = rev(c("Pooled", "GLP1", "DPP4", "SGLT2", "TZD", "SU")), labels = rev(c("Pooled", "GLP-1RA", "DPP4i", "SGLT2i", "TZD", "SU")),
                        name = "Therapy", guide = guide_legend(reverse = TRUE)) +
    theme_bw() +
    theme(
      axis.title.y = element_blank(),
      legend.position = "bottom",
      strip.text = element_blank()
    ),
  
  
  coefficients_12m %>%
    filter(type == "Clinical features & biomarkers") %>%
    filter(variable %in% c("Diabetes duration (per SD)", "Age (per SD)", "HbA1c (per SD)", "BMI (per SD)", "Sex (ref Male)")) %>%
    mutate(variable = factor(variable, levels = c("Age (per SD)", "Diabetes duration (per SD)", "BMI (per SD)", "HbA1c (per SD)", "Sex (ref Male)"))) %>%
    ggplot() + 
    geom_vline(aes(xintercept = 1), colour = "black") +
    geom_errorbar(aes(x = coef, xmin = LCI, xmax = UCI, y = variable, colour = drug),position = position_dodge(width = 1)) +
    geom_pointrange(aes(x = coef, xmin = LCI, xmax = UCI, y = variable, colour = drug), position = position_dodge(width = 1)) +
    xlab("Odds Ratio") +
    scale_x_log10(breaks = c(0.25, 0.5, 1, 2, 4), limits = c(0.4, 4.1)) +
    facet_grid(variable~., scales = "free_y") +
    scale_colour_manual(values = c("Pooled" = "black", "SGLT2" = "#E69F00", "GLP1" = "#56B4E9", "SU" = "#CC79A7", "DPP4" = "#0072B2", "TZD" = "#D55E00", "MFN" = "grey"), 
                        breaks = rev(c("Pooled", "GLP1", "DPP4", "SGLT2", "TZD", "SU")), labels = rev(c("Pooled", "GLP-1RA", "DPP4i", "SGLT2i", "TZD", "SU")),
                        name = "Therapy", guide = guide_legend(reverse = TRUE)) +
    theme_bw() +
    theme(
      axis.title.y = element_blank(),
      legend.position = "bottom",
      strip.text = element_blank()
    ),
  
  coefficients_12m %>%
    filter(type == "Comorbidity event") %>%
    filter(variable %in% c("Frailty")) %>%
    ggplot() + 
    geom_vline(aes(xintercept = 1), colour = "black") +
    geom_errorbar(aes(x = coef, xmin = LCI, xmax = UCI, y = variable, colour = drug),position = position_dodge(width = 1)) +
    geom_pointrange(aes(x = coef, xmin = LCI, xmax = UCI, y = variable, colour = drug), position = position_dodge(width = 1)) +
    xlab("Odds Ratio") +
    scale_x_log10(breaks = c(0.25, 0.5, 1, 2, 4), limits = c(0.4, 4.1)) +
    facet_grid(variable~., scales = "free_y") +
    scale_colour_manual(values = c("Pooled" = "black", "SGLT2" = "#E69F00", "GLP1" = "#56B4E9", "SU" = "#CC79A7", "DPP4" = "#0072B2", "TZD" = "#D55E00", "MFN" = "grey"), 
                        breaks = rev(c("Pooled", "GLP1", "DPP4", "SGLT2", "TZD", "SU")), labels = rev(c("Pooled", "GLP-1RA", "DPP4i", "SGLT2i", "TZD", "SU")),
                        name = "Therapy", guide = guide_legend(reverse = TRUE)) +
    theme_bw() +
    theme(
      axis.title.y = element_blank(),
      legend.position = "bottom",
      strip.text = element_blank()
    ),
  
  ncol = 1, nrow = 3
  
) + 
  plot_layout(guides = "collect", height = c(5, 5, 1), axis_titles = "collect", axes = "collect") +
  plot_annotation(tag_levels = list(c("Behavioural", "Clinical features & biomarkers", "Comorbidity event"))) & 
  theme(
    legend.direction = "vertical",
    axis.text.y = element_text(size = 15),
    axis.text.x = element_text(size = 15),
    axis.title.x = element_text(size = 18)
  )

dev.off()
