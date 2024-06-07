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
  raw_data = "20240308_t2d_1stinstance",
  diagnosis = FALSE,
  therapies = c("DPP4", "GLP1", "MFN", "SGLT2", "SU", "TZD"),
  dataset = "3m.disc.dataset",
  full_prescribing_history = TRUE
) %>%
  # drop_na(-stopdrug_6m_6mFU, -stopdrug_12m_6mFU)
  drop_na(stopdrug_3m_6mFU)




###############################################################################
###############################################################################

# plot for the therapy discontinuation

pdf("results/figures/03.discontinuation_per_drug_combined.pdf", width = 6, height = 5)

data <- cprd_dataset %>%
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
    drugclass = factor(drugclass, levels = c("MFN", "GLP1", "DPP4", "SGLT2", "TZD", "SU")),
    key = factor(key, levels = c("perc_12m", "perc_6m", "perc_3m"))
  ) %>%
  ggplot(aes(x = value, y = drugclass, fill = drugclass, alpha = key)) +
  geom_bar(stat = "identity", colour = "black") +
  scale_fill_manual(values = c("Pooled" = "black", "SGLT2" = "#E69F00", "GLP1" = "#56B4E9", "SU" = "#CC79A7", "DPP4" = "#0072B2", "TZD" = "#D55E00", "MFN" = "grey"), 
                    breaks = rev(c("Pooled", "MFN", "GLP1", "DPP4", "SGLT2", "TZD", "SU")), labels = rev(c("Pooled", "Metformin", "GLP-1RA", "DPP4i", "SGLT2i", "TZD", "SU")),
                    name = "Therapy", guide = guide_legend(reverse = TRUE)) +
  scale_y_discrete(limits = rev, breaks = c("MFN", "GLP1", "DPP4", "SGLT2", "TZD", "SU"), labels = c("Metformin", "GLP-1RA", "DPP4i", "SGLT2i", "TZD", "SU")) +
  scale_alpha_discrete(limits = rev, range = c(1, 0.4), breaks = c("perc_12m", "perc_6m", "perc_3m"), labels = c("12-month", "6-month", "3-month"), name = "Discontinuation", guide = guide_legend(reverse = TRUE, row = 1)) +
  scale_x_continuous("Proportion of discontinuation", labels = scales::percent, limits = c(0, 0.35), breaks = seq(0, 0.35, 0.05)) +
  ggtitle("Therapies (2014-onwards)") +
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



cprd_dataset %>%
  filter(dstartdate > "2017-01-01") %>%
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
      filter(dstartdate > "2017-01-01") %>%
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
      filter(dstartdate > "2017-01-01") %>%
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
    drugclass = factor(drugclass, levels = c("MFN", "GLP1", "DPP4", "SGLT2", "TZD", "SU")),
    key = factor(key, levels = c("perc_12m", "perc_6m", "perc_3m"))
  ) %>%
  ggplot(aes(x = value, y = drugclass, fill = drugclass, alpha = key)) +
  geom_bar(stat = "identity", colour = "black") +
  scale_fill_manual(values = c("Pooled" = "black", "SGLT2" = "#E69F00", "GLP1" = "#56B4E9", "SU" = "#CC79A7", "DPP4" = "#0072B2", "TZD" = "#D55E00", "MFN" = "grey"), 
                    breaks = rev(c("Pooled", "MFN", "GLP1", "DPP4", "SGLT2", "TZD", "SU")), labels = rev(c("Pooled", "Metformin", "GLP-1RA", "DPP4i", "SGLT2i", "TZD", "SU")),
                    name = "Therapy", guide = guide_legend(reverse = TRUE)) +
  scale_y_discrete(limits = rev, breaks = c("MFN", "GLP1", "DPP4", "SGLT2", "TZD", "SU"), labels = c("Metformin", "GLP-1RA", "DPP4i", "SGLT2i", "TZD", "SU")) +
  scale_alpha_discrete(limits = rev, range = c(1, 0.4), breaks = c("perc_12m", "perc_6m", "perc_3m"), labels = c("12-month", "6-month", "3-month"), name = "Discontinuation", guide = guide_legend(reverse = TRUE, row = 1)) +
  scale_x_continuous("Proportion of discontinuation", labels = scales::percent, limits = c(0, 0.35), breaks = seq(0, 0.35, 0.05)) +
  ggtitle("Therapies (2017-onwards)") +
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


dev.off()

###############################################################################
###############################################################################


# Start plot by variables


pdf("results/figures/03.3m_discontinuation_by_variable.pdf", width = 7, height = 13)

patchwork::wrap_plots(
  
  # Overall
  cprd_dataset %>%
    select(stopdrug_3m_6mFU, drugclass) %>%
    mutate(drugclass = factor(drugclass, levels = c("MFN", "GLP1", "DPP4", "SGLT2", "TZD", "SU"), labels = c("MFN", "other", "other", "other", "other", "other"))) %>%
    rename("analysis_var" = "stopdrug_3m_6mFU") %>%
    table() %>%
    as.data.frame() %>%
    group_by(drugclass) %>%
    mutate(total = sum(Freq), Freq = Freq/total) %>%
    filter(analysis_var == 1) %>%
    ggplot(aes(y = analysis_var, x = Freq, fill = drugclass)) +
    geom_col() +
    scale_fill_manual(values = c("other" = "black", "SGLT2" = "#E69F00", "GLP1" = "#56B4E9", "SU" = "#CC79A7", "DPP4" = "#0072B2", "TZD" = "#D55E00", "MFN" = "grey"), breaks = c("MFN", "other"), labels = c("Metformin", "Other glucose-lowering therapies"), name = " ", guide = guide_legend(reverse = FALSE)) +
    scale_y_discrete(breaks = c(1), labels = c("Overall")) +
    scale_x_continuous("Proportion of discontinuation", labels = scales::percent, limits = c(0, 0.3)) +
    facet_wrap(~drugclass, nrow = 1, labeller = as_labeller(c("MFN" = "Metformin", "other" = "Other glucose-lowering therapies"))) +
    theme_bw() +
    guides(fill = guide_legend(nrow = 1, byrow = TRUE)) +
    theme(
      legend.position = "bottom",
      axis.text.x = element_blank(),
      axis.title.x = element_blank()
    ),
  
  # Age
  cprd_dataset %>%
    select(dstartdate_age, stopdrug_3m_6mFU, drugclass) %>%
    mutate(drugclass = factor(drugclass, levels = c("MFN", "GLP1", "DPP4", "SGLT2", "TZD", "SU"), labels = c("MFN", "other", "other", "other", "other", "other"))) %>%
    mutate(group = cut(dstartdate_age, breaks = c(18, 50, 60, 70, 150))) %>%
    rename("analysis_var" = "stopdrug_3m_6mFU") %>%
    select(-dstartdate_age) %>%
    table() %>%
    as.data.frame() %>%
    group_by(group, drugclass) %>%
    mutate(total = sum(Freq), Freq = Freq/total) %>%
    filter(analysis_var == 1) %>%
    ggplot(aes(y = group, x = Freq, fill = drugclass)) +
    geom_col() +
    scale_fill_manual(values = c("other" = "black", "SGLT2" = "#E69F00", "GLP1" = "#56B4E9", "SU" = "#CC79A7", "DPP4" = "#0072B2", "TZD" = "#D55E00", "MFN" = "grey"), breaks = c("MFN", "other"), labels = c("Metformin", "Other glucose-lowering therapies"), name = " ", guide = guide_legend(reverse = FALSE)) +
    scale_y_discrete(limits = rev, labels = c("70+", "60-69", "50-59", "<50")) +
    scale_x_continuous("Proportion of discontinuation", labels = scales::percent, limits = c(0, 0.3)) +
    ggtitle("Age group, years") +
    facet_wrap(~drugclass, nrow = 1, labeller = as_labeller(c("MFN" = "Metformin", "other" = "Other glucose-lowering therapies"))) +
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
  cprd_dataset %>%
    select(gender, stopdrug_3m_6mFU, drugclass) %>%
    mutate(drugclass = factor(drugclass, levels = c("MFN", "GLP1", "DPP4", "SGLT2", "TZD", "SU"), labels = c("MFN", "other", "other", "other", "other", "other"))) %>%
    rename("analysis_var" = "stopdrug_3m_6mFU") %>%
    table() %>%
    as.data.frame() %>%
    group_by(gender, drugclass) %>%
    mutate(total = sum(Freq), Freq = Freq/total) %>%
    filter(analysis_var == 1) %>%
    ggplot(aes(y = gender, x = Freq, fill = drugclass)) +
    geom_col() +
    scale_fill_manual(values = c("other" = "black", "SGLT2" = "#E69F00", "GLP1" = "#56B4E9", "SU" = "#CC79A7", "DPP4" = "#0072B2", "TZD" = "#D55E00", "MFN" = "grey"), breaks = c("MFN", "other"), labels = c("Metformin", "Other glucose-lowering therapies"), name = " ", guide = guide_legend(reverse = FALSE)) +
    scale_y_discrete(limits = rev, labels = c("Female", "Male")) +
    scale_x_continuous("Proportion of discontinuation", labels = scales::percent, limits = c(0, 0.3)) +
    ggtitle("Sex") +
    facet_wrap(~drugclass, nrow = 1, labeller = as_labeller(c("MFN" = "Metformin", "other" = "Other glucose-lowering therapies"))) +
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
  cprd_dataset %>%
    select(imd2015_10, stopdrug_3m_6mFU, drugclass) %>%
    mutate(drugclass = factor(drugclass, levels = c("MFN", "GLP1", "DPP4", "SGLT2", "TZD", "SU"), labels = c("MFN", "other", "other", "other", "other", "other"))) %>%
    rename("analysis_var" = "stopdrug_3m_6mFU") %>%
    table() %>%
    as.data.frame() %>%
    group_by(imd2015_10, drugclass) %>%
    mutate(total = sum(Freq), Freq = Freq/total) %>%
    filter(analysis_var == 1) %>%
    ggplot(aes(y = imd2015_10, x = Freq, fill = drugclass)) +
    geom_col() +
    scale_fill_manual(values = c("other" = "black", "SGLT2" = "#E69F00", "GLP1" = "#56B4E9", "SU" = "#CC79A7", "DPP4" = "#0072B2", "TZD" = "#D55E00", "MFN" = "grey"), breaks = c("MFN", "other"), labels = c("Metformin", "Other glucose-lowering therapies"), name = " ", guide = guide_legend(reverse = FALSE)) +
    scale_y_discrete(limits = rev, labels = c("5 (most deprived)", "4", "3", "2", "1 (least deprived)")) +
    scale_x_continuous("Proportion of discontinuation", labels = scales::percent, limits = c(0, 0.3)) +
    ggtitle("Index of multiple deprivation quintile") +
    facet_wrap(~drugclass, nrow = 1, labeller = as_labeller(c("MFN" = "Metformin", "other" = "Other glucose-lowering therapies"))) +
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
  cprd_dataset %>%
    select(prebmi, stopdrug_3m_6mFU, drugclass) %>%
    mutate(drugclass = factor(drugclass, levels = c("MFN", "GLP1", "DPP4", "SGLT2", "TZD", "SU"), labels = c("MFN", "other", "other", "other", "other", "other"))) %>%
    mutate(group = cut(prebmi, breaks = c(1, 25, 30, 35, 500))) %>%
    rename("analysis_var" = "stopdrug_3m_6mFU") %>%
    select(-prebmi) %>%
    table() %>%
    as.data.frame() %>%
    group_by(group, drugclass) %>%
    mutate(total = sum(Freq), Freq = Freq/total) %>%
    filter(analysis_var == 1) %>%
    ggplot(aes(y = group, x = Freq, fill = drugclass)) +
    geom_col() +
    scale_fill_manual(values = c("other" = "black", "SGLT2" = "#E69F00", "GLP1" = "#56B4E9", "SU" = "#CC79A7", "DPP4" = "#0072B2", "TZD" = "#D55E00", "MFN" = "grey"), breaks = c("MFN", "other"), labels = c("Metformin", "Other glucose-lowering therapies"), name = " ", guide = guide_legend(reverse = FALSE)) +
    scale_y_discrete(limits = rev, labels = c("35+", "30-34.9", "25-29.9", "<25")) +
    ggtitle("BMI, kg/m2") +
    scale_x_continuous("Proportion of discontinuation", labels = scales::percent, limits = c(0, 0.3)) +
    facet_wrap(~drugclass, nrow = 1, labeller = as_labeller(c("MFN" = "Metformin", "other" = "Other glucose-lowering therapies"))) +
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
  cprd_dataset %>%
    select(dstartdate_dm_dur, stopdrug_3m_6mFU, drugclass) %>%
    mutate(drugclass = factor(drugclass, levels = c("MFN", "GLP1", "DPP4", "SGLT2", "TZD", "SU"), labels = c("MFN", "other", "other", "other", "other", "other"))) %>%
    mutate(group = cut(dstartdate_dm_dur, breaks = c(-1, 1, 3, 6, 10, 150))) %>%
    rename("analysis_var" = "stopdrug_3m_6mFU") %>%
    select(-dstartdate_dm_dur) %>%
    table() %>%
    as.data.frame() %>%
    group_by(group, drugclass) %>%
    mutate(total = sum(Freq), Freq = Freq/total) %>%
    filter(analysis_var == 1) %>%
    ggplot(aes(y = group, x = Freq, fill = drugclass)) +
    geom_col() +
    scale_fill_manual(values = c("other" = "black", "SGLT2" = "#E69F00", "GLP1" = "#56B4E9", "SU" = "#CC79A7", "DPP4" = "#0072B2", "TZD" = "#D55E00", "MFN" = "grey"), breaks = c("MFN", "other"), labels = c("Metformin", "Other glucose-lowering therapies"), name = " ", guide = guide_legend(reverse = FALSE)) +
    scale_y_discrete(limits = rev, labels = c("10+", "6-9", "3-5", "1-2", "<1")) +
    ggtitle("Duration of diabetes, years") +
    scale_x_continuous("Proportion of discontinuation", labels = scales::percent, limits = c(0, 0.3)) +
    facet_wrap(~drugclass, nrow = 1, labeller = as_labeller(c("MFN" = "Metformin", "other" = "Other glucose-lowering therapies"))) +
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
  cprd_dataset %>%
    select(prehba1c, stopdrug_3m_6mFU, drugclass) %>%
    mutate(drugclass = factor(drugclass, levels = c("MFN", "GLP1", "DPP4", "SGLT2", "TZD", "SU"), labels = c("MFN", "other", "other", "other", "other", "other"))) %>%
    mutate(group = cut(prehba1c, breaks = c(52, 64, 75, 86, 250))) %>%
    rename("analysis_var" = "stopdrug_3m_6mFU") %>%
    select(-prehba1c) %>%
    table() %>%
    as.data.frame() %>%
    group_by(group, drugclass) %>%
    mutate(total = sum(Freq), Freq = Freq/total) %>%
    filter(analysis_var == 1) %>%
    ggplot(aes(y = group, x = Freq, fill = drugclass)) +
    geom_col() +
    scale_fill_manual(values = c("other" = "black", "SGLT2" = "#E69F00", "GLP1" = "#56B4E9", "SU" = "#CC79A7", "DPP4" = "#0072B2", "TZD" = "#D55E00", "MFN" = "grey"), breaks = c("MFN", "other"), labels = c("Metformin", "Other glucose-lowering therapies"), name = " ", guide = guide_legend(reverse = FALSE)) +
    scale_y_discrete(limits = rev, labels = c("86+", "75-86", "64-75", "53-64")) +
    ggtitle("HbA1c, mmol/mol") +
    scale_x_continuous("Proportion of discontinuation", labels = scales::percent, limits = c(0, 0.3)) +
    facet_wrap(~drugclass, nrow = 1, labeller = as_labeller(c("MFN" = "Metformin", "other" = "Other glucose-lowering therapies"))) +
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
  cprd_dataset %>%
    select(drugline, stopdrug_3m_6mFU, drugclass) %>%
    mutate(drugline = factor(drugline, levels = c("1", "2", "3", "4", "5+"), labels = c("1", "2", "3+", "3+", "3+"))) %>%
    mutate(drugclass = factor(drugclass, levels = c("MFN", "GLP1", "DPP4", "SGLT2", "TZD", "SU"), labels = c("MFN", "other", "other", "other", "other", "other"))) %>%
    rename("analysis_var" = "stopdrug_3m_6mFU") %>%
    table() %>%
    as.data.frame() %>%
    group_by(drugline, drugclass) %>%
    mutate(total = sum(Freq), Freq = Freq/total) %>%
    filter(analysis_var == 1) %>%
    ggplot(aes(y = drugline, x = Freq, fill = drugclass)) +
    geom_col() +
    scale_fill_manual(values = c("other" = "black", "SGLT2" = "#E69F00", "GLP1" = "#56B4E9", "SU" = "#CC79A7", "DPP4" = "#0072B2", "TZD" = "#D55E00", "MFN" = "grey"), breaks = c("MFN", "other"), labels = c("Metformin", "Other glucose-lowering therapies"), name = " ", guide = guide_legend(reverse = FALSE)) +
    scale_y_discrete(limits = rev) +
    ggtitle("Number of therapy classes ever prescribed") +
    scale_x_continuous("Proportion of discontinuation", labels = scales::percent, limits = c(0, 0.3)) +
    facet_wrap(~drugclass, nrow = 1, labeller = as_labeller(c("MFN" = "Metformin", "other" = "Other glucose-lowering therapies"))) +
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
  cprd_dataset %>%
    select(predrug_frailty_proxy, stopdrug_3m_6mFU, drugclass) %>%
    mutate(drugclass = factor(drugclass, levels = c("MFN", "GLP1", "DPP4", "SGLT2", "TZD", "SU"), labels = c("MFN", "other", "other", "other", "other", "other"))) %>%
    rename("analysis_var" = "stopdrug_3m_6mFU") %>%
    table() %>%
    as.data.frame() %>%
    group_by(predrug_frailty_proxy, drugclass) %>%
    mutate(total = sum(Freq), Freq = Freq/total) %>%
    filter(analysis_var == 1) %>%
    ggplot(aes(y = predrug_frailty_proxy, x = Freq, fill = drugclass)) +
    geom_col() +
    scale_fill_manual(values = c("other" = "black", "SGLT2" = "#E69F00", "GLP1" = "#56B4E9", "SU" = "#CC79A7", "DPP4" = "#0072B2", "TZD" = "#D55E00", "MFN" = "grey"), breaks = c("MFN", "other"), labels = c("Metformin", "Other glucose-lowering therapies"), name = " ", guide = guide_legend(reverse = FALSE)) +
    scale_y_discrete(limits = rev, breaks = c(0, 1), labels = c("No", "Yes")) +
    ggtitle("Frailty") +
    scale_x_continuous("Proportion of discontinuation", labels = scales::percent, limits = c(0, 0.3)) +
    facet_wrap(~drugclass, nrow = 1, labeller = as_labeller(c("MFN" = "Metformin", "other" = "Other glucose-lowering therapies"))) +
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
  cprd_dataset %>%
    select(ethnicity_5cat, stopdrug_3m_6mFU, drugclass) %>%
    mutate(ethnicity_5cat = factor(ethnicity_5cat, levels = c("4", "3", "2", "1", "0"), labels = c("Mixed & Other", "Mixed & Other", "Black", "South Asian", "White"))) %>%
    mutate(drugclass = factor(drugclass, levels = c("MFN", "GLP1", "DPP4", "SGLT2", "TZD", "SU"), labels = c("MFN", "other", "other", "other", "other", "other"))) %>%
    rename("analysis_var" = "stopdrug_3m_6mFU") %>%
    table() %>%
    as.data.frame() %>%
    group_by(ethnicity_5cat, drugclass) %>%
    mutate(total = sum(Freq), Freq = Freq/total) %>%
    filter(analysis_var == 1) %>%
    ggplot(aes(y = ethnicity_5cat, x = Freq, fill = drugclass)) +
    geom_col() +
    scale_fill_manual(values = c("other" = "black", "SGLT2" = "#E69F00", "GLP1" = "#56B4E9", "SU" = "#CC79A7", "DPP4" = "#0072B2", "TZD" = "#D55E00", "MFN" = "grey"), breaks = c("MFN", "other"), labels = c("Metformin", "Other glucose-lowering therapies"), name = " ", guide = guide_legend(reverse = FALSE)) +
    ggtitle("Ethnicity") +
    scale_x_continuous("Proportion of discontinuation", labels = scales::percent, limits = c(0, 0.3)) +
    facet_wrap(~drugclass, nrow = 1, labeller = as_labeller(c("MFN" = "Metformin", "other" = "Other glucose-lowering therapies"))) +
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
    height = c(1, 4, 2, 5, 4, 5, 4, 3, 2, 4),
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


pdf("results/figures/03.3m_discontinuation_by_variable_by_drug.pdf", width = 14, height = 13)

patchwork::wrap_plots(
  
  # Age - Metformin
  cprd_dataset %>%
    select(dstartdate_age, stopdrug_3m_6mFU, drugclass) %>%
    mutate(drugclass = factor(drugclass, levels = c("MFN", "GLP1", "DPP4", "SGLT2", "TZD", "SU"))) %>%
    mutate(group = cut(dstartdate_age, breaks = c(18, 50, 60, 70, 150))) %>%
    rename("analysis_var" = "stopdrug_3m_6mFU") %>%
    select(-dstartdate_age) %>%
    table() %>%
    as.data.frame() %>%
    group_by(group, drugclass) %>%
    mutate(total = sum(Freq), Freq = Freq/total) %>%
    filter(analysis_var == 1) %>%
    ggplot(aes(y = group, x = Freq, fill = drugclass)) +
    geom_col(position = "dodge") +
    scale_fill_manual(values = c("Pooled" = "black", "SGLT2" = "#E69F00", "GLP1" = "#56B4E9", "SU" = "#CC79A7", "DPP4" = "#0072B2", "TZD" = "#D55E00", "MFN" = "grey"), breaks = c("MFN", "GLP1", "DPP4", "SGLT2", "TZD", "SU"), labels = c("Metformin", "GLP-1RA", "DPP4i", "SGLT2i", "TZD", "SU"), name = " ", guide = guide_legend(reverse = TRUE)) +
    scale_y_discrete(limits = rev, labels = c("70+", "60-69", "50-59", "<50")) +
    ggtitle("Age group, years") +
    scale_x_continuous("Proportion of discontinuation", labels = scales::percent, limits = c(0, 0.5), breaks = seq(0, 0.5, 0.1)) +
    facet_wrap(~drugclass, nrow = 1, labeller = as_labeller(c("MFN" = "Metformin", "GLP1" = "GLP-1RA", "DPP4" = "DPP4i", "SGLT2" = "SGLT2i", "TZD" = "TZD", "SU" = "SU"))) +
    theme_bw() +
    guides(fill = guide_legend(nrow = 1, byrow = TRUE)) +
    theme(
      legend.position = "bottom",
      axis.text.x = element_blank(),
      axis.title.x = element_blank()
    ),
  
  # Sex 
  cprd_dataset %>%
    select(gender, stopdrug_3m_6mFU, drugclass) %>%
    mutate(drugclass = factor(drugclass, levels = c("MFN", "GLP1", "DPP4", "SGLT2", "TZD", "SU"))) %>%
    rename("analysis_var" = "stopdrug_3m_6mFU") %>%
    table() %>%
    as.data.frame() %>%
    group_by(gender, drugclass) %>%
    mutate(total = sum(Freq), Freq = Freq/total) %>%
    filter(analysis_var == 1) %>%
    ggplot(aes(y = gender, x = Freq, fill = drugclass)) +
    geom_col(position = "dodge") +
    scale_fill_manual(values = c("Pooled" = "black", "SGLT2" = "#E69F00", "GLP1" = "#56B4E9", "SU" = "#CC79A7", "DPP4" = "#0072B2", "TZD" = "#D55E00", "MFN" = "grey"), breaks = c("MFN", "GLP1", "DPP4", "SGLT2", "TZD", "SU"), labels = c("Metformin", "GLP-1RA", "DPP4i", "SGLT2i", "TZD", "SU"), name = " ", guide = guide_legend(reverse = TRUE)) +
    scale_y_discrete(limits = rev, labels = c("Female", "Male")) +
    ggtitle("Sex") +
    scale_x_continuous("Proportion of discontinuation", labels = scales::percent, limits = c(0, 0.5), breaks = seq(0, 0.5, 0.1)) +
    facet_wrap(~drugclass, nrow = 1, labeller = as_labeller(c("MFN" = "Metformin", "GLP1" = "GLP-1RA", "DPP4" = "DPP4i", "SGLT2" = "SGLT2i", "TZD" = "TZD", "SU" = "SU"))) +
    theme_bw() +
    guides(fill = guide_legend(nrow = 1, byrow = TRUE)) +
    theme(
      legend.position = "bottom",
      axis.text.x = element_blank(),
      axis.title.x = element_blank(),
      strip.text.x = element_blank()
    ),
  
  # Quintile
  cprd_dataset %>%
    select(imd2015_10, stopdrug_3m_6mFU, drugclass) %>%
    mutate(drugclass = factor(drugclass, levels = c("MFN", "GLP1", "DPP4", "SGLT2", "TZD", "SU"))) %>%
    rename("analysis_var" = "stopdrug_3m_6mFU") %>%
    table() %>%
    as.data.frame() %>%
    group_by(imd2015_10, drugclass) %>%
    mutate(total = sum(Freq), Freq = Freq/total) %>%
    filter(analysis_var == 1) %>%
    ggplot(aes(y = imd2015_10, x = Freq, fill = drugclass)) +
    geom_col(position = "dodge") +
    scale_fill_manual(values = c("Pooled" = "black", "SGLT2" = "#E69F00", "GLP1" = "#56B4E9", "SU" = "#CC79A7", "DPP4" = "#0072B2", "TZD" = "#D55E00", "MFN" = "grey"), breaks = c("MFN", "GLP1", "DPP4", "SGLT2", "TZD", "SU"), labels = c("Metformin", "GLP-1RA", "DPP4i", "SGLT2i", "TZD", "SU"), name = " ", guide = guide_legend(reverse = TRUE)) +
    scale_y_discrete(limits = rev, labels = c("5 (most deprived)", "4", "3", "2", "1 (least deprived)")) +
    ggtitle("Index of multiple deprivation quintile") +
    scale_x_continuous("Proportion of discontinuation", labels = scales::percent, limits = c(0, 0.5), breaks = seq(0, 0.5, 0.1)) +
    facet_wrap(~drugclass, nrow = 1, labeller = as_labeller(c("MFN" = "Metformin", "GLP1" = "GLP-1RA", "DPP4" = "DPP4i", "SGLT2" = "SGLT2i", "TZD" = "TZD", "SU" = "SU"))) +
    theme_bw() +
    guides(fill = guide_legend(nrow = 1, byrow = TRUE)) +
    theme(
      legend.position = "bottom",
      axis.text.x = element_blank(),
      axis.title.x = element_blank(),
      strip.text.x = element_blank()
    ),
  
  # BMI
  cprd_dataset %>%
    select(prebmi, stopdrug_3m_6mFU, drugclass) %>%
    mutate(drugclass = factor(drugclass, levels = c("MFN", "GLP1", "DPP4", "SGLT2", "TZD", "SU"))) %>%
    mutate(group = cut(prebmi, breaks = c(1, 25, 30, 35, 500))) %>%
    rename("analysis_var" = "stopdrug_3m_6mFU") %>%
    select(-prebmi) %>%
    table() %>%
    as.data.frame() %>%
    group_by(group, drugclass) %>%
    mutate(total = sum(Freq), Freq = Freq/total) %>%
    filter(analysis_var == 1) %>%
    ggplot(aes(y = group, x = Freq, fill = drugclass)) +
    geom_col(position = "dodge") +
    scale_fill_manual(values = c("Pooled" = "black", "SGLT2" = "#E69F00", "GLP1" = "#56B4E9", "SU" = "#CC79A7", "DPP4" = "#0072B2", "TZD" = "#D55E00", "MFN" = "grey"), breaks = c("MFN", "GLP1", "DPP4", "SGLT2", "TZD", "SU"), labels = c("Metformin", "GLP-1RA", "DPP4i", "SGLT2i", "TZD", "SU"), name = " ", guide = guide_legend(reverse = TRUE)) +
    scale_y_discrete(limits = rev, labels = c("35+", "30-34.9", "25-29.9", "<25")) +
    ggtitle("BMI, kg/m2") +
    scale_x_continuous("Proportion of discontinuation", labels = scales::percent, limits = c(0, 0.5), breaks = seq(0, 0.5, 0.1)) +
    facet_wrap(~drugclass, nrow = 1, labeller = as_labeller(c("MFN" = "Metformin", "GLP1" = "GLP-1RA", "DPP4" = "DPP4i", "SGLT2" = "SGLT2i", "TZD" = "TZD", "SU" = "SU"))) +
    theme_bw() +
    guides(fill = guide_legend(nrow = 1, byrow = TRUE)) +
    theme(
      legend.position = "bottom",
      axis.text.x = element_blank(),
      axis.title.x = element_blank(),
      strip.text.x = element_blank()
    ),
  
  # Diabetes duration
  cprd_dataset %>%
    select(dstartdate_dm_dur, stopdrug_3m_6mFU, drugclass) %>%
    mutate(drugclass = factor(drugclass, levels = c("MFN", "GLP1", "DPP4", "SGLT2", "TZD", "SU"))) %>%
    mutate(group = cut(dstartdate_dm_dur, breaks = c(-1, 1, 3, 6, 10, 150))) %>%
    rename("analysis_var" = "stopdrug_3m_6mFU") %>%
    select(-dstartdate_dm_dur) %>%
    table() %>%
    as.data.frame() %>%
    group_by(group, drugclass) %>%
    mutate(total = sum(Freq), Freq = Freq/total) %>%
    filter(analysis_var == 1) %>%
    ggplot(aes(y = group, x = Freq, fill = drugclass)) +
    geom_col(position = "dodge") +
    scale_fill_manual(values = c("Pooled" = "black", "SGLT2" = "#E69F00", "GLP1" = "#56B4E9", "SU" = "#CC79A7", "DPP4" = "#0072B2", "TZD" = "#D55E00", "MFN" = "grey"), breaks = c("MFN", "GLP1", "DPP4", "SGLT2", "TZD", "SU"), labels = c("Metformin", "GLP-1RA", "DPP4i", "SGLT2i", "TZD", "SU"), name = " ", guide = guide_legend(reverse = TRUE)) +
    scale_y_discrete(limits = rev, labels = c("10+", "6-9", "3-5", "1-2", "<1")) +
    ggtitle("Duration of diabetes, years") +
    scale_x_continuous("Proportion of discontinuation", labels = scales::percent, limits = c(0, 0.5), breaks = seq(0, 0.5, 0.1)) +
    facet_wrap(~drugclass, nrow = 1, labeller = as_labeller(c("MFN" = "Metformin", "GLP1" = "GLP-1RA", "DPP4" = "DPP4i", "SGLT2" = "SGLT2i", "TZD" = "TZD", "SU" = "SU"))) +
    theme_bw() +
    guides(fill = guide_legend(nrow = 1, byrow = TRUE)) +
    theme(
      legend.position = "bottom",
      axis.text.x = element_blank(),
      axis.title.x = element_blank(),
      strip.text.x = element_blank()
    ),
  
  # HbA1c
  cprd_dataset %>%
    select(prehba1c, stopdrug_3m_6mFU, drugclass) %>%
    mutate(drugclass = factor(drugclass, levels = c("MFN", "GLP1", "DPP4", "SGLT2", "TZD", "SU"))) %>%
    mutate(group = cut(prehba1c, breaks = c(52, 64, 75, 86, 250))) %>%
    rename("analysis_var" = "stopdrug_3m_6mFU") %>%
    select(-prehba1c) %>%
    table() %>%
    as.data.frame() %>%
    group_by(group, drugclass) %>%
    mutate(total = sum(Freq), Freq = Freq/total) %>%
    filter(analysis_var == 1) %>%
    ggplot(aes(y = group, x = Freq, fill = drugclass)) +
    geom_col(position = "dodge") +
    scale_fill_manual(values = c("Pooled" = "black", "SGLT2" = "#E69F00", "GLP1" = "#56B4E9", "SU" = "#CC79A7", "DPP4" = "#0072B2", "TZD" = "#D55E00", "MFN" = "grey"), breaks = c("MFN", "GLP1", "DPP4", "SGLT2", "TZD", "SU"), labels = c("Metformin", "GLP-1RA", "DPP4i", "SGLT2i", "TZD", "SU"), name = " ", guide = guide_legend(reverse = TRUE)) +
    scale_y_discrete(limits = rev, labels = c("86+", "75-86", "64-75", "53-64")) +
    ggtitle("HbA1c, mmol/mol") +
    scale_x_continuous("Proportion of discontinuation", labels = scales::percent, limits = c(0, 0.5), breaks = seq(0, 0.5, 0.1)) +
    facet_wrap(~drugclass, nrow = 1, labeller = as_labeller(c("MFN" = "Metformin", "GLP1" = "GLP-1RA", "DPP4" = "DPP4i", "SGLT2" = "SGLT2i", "TZD" = "TZD", "SU" = "SU"))) +
    theme_bw() +
    guides(fill = guide_legend(nrow = 1, byrow = TRUE)) +
    theme(
      legend.position = "bottom",
      axis.text.x = element_blank(),
      axis.title.x = element_blank(),
      strip.text.x = element_blank()
    ),
  
  # drugline
  cprd_dataset %>%
    select(drugline, stopdrug_3m_6mFU, drugclass) %>%
    mutate(drugclass = factor(drugclass, levels = c("MFN", "GLP1", "DPP4", "SGLT2", "TZD", "SU"))) %>%
    mutate(drugline = factor(drugline, levels = c("1", "2", "3", "4", "5+"), labels = c("1", "2", "3+", "3+", "3+"))) %>%
    rename("analysis_var" = "stopdrug_3m_6mFU") %>%
    table() %>%
    as.data.frame() %>%
    group_by(drugline, drugclass) %>%
    mutate(total = sum(Freq), Freq = Freq/total) %>%
    filter(analysis_var == 1) %>%
    ggplot(aes(y = drugline, x = Freq, fill = drugclass)) +
    geom_col(position = "dodge") +
    scale_fill_manual(values = c("Pooled" = "black", "SGLT2" = "#E69F00", "GLP1" = "#56B4E9", "SU" = "#CC79A7", "DPP4" = "#0072B2", "TZD" = "#D55E00", "MFN" = "grey"), breaks = c("MFN", "GLP1", "DPP4", "SGLT2", "TZD", "SU"), labels = c("Metformin", "GLP-1RA", "DPP4i", "SGLT2i", "TZD", "SU"), name = " ", guide = guide_legend(reverse = TRUE)) +
    scale_y_discrete(limits = rev) +
    ggtitle("Number of therapy classes ever prescribed") +
    scale_x_continuous("Proportion of discontinuation", labels = scales::percent, limits = c(0, 0.5), breaks = seq(0, 0.5, 0.1)) +
    facet_wrap(~drugclass, nrow = 1, labeller = as_labeller(c("MFN" = "Metformin", "GLP1" = "GLP-1RA", "DPP4" = "DPP4i", "SGLT2" = "SGLT2i", "TZD" = "TZD", "SU" = "SU"))) +
    theme_bw() +
    guides(fill = guide_legend(nrow = 1, byrow = TRUE)) +
    theme(
      legend.position = "bottom",
      axis.text.x = element_blank(),
      axis.title.x = element_blank(),
      strip.text.x = element_blank()
    ),
  
  # Frailty
  cprd_dataset %>%
    select(predrug_frailty_proxy, stopdrug_3m_6mFU, drugclass) %>%
    mutate(drugclass = factor(drugclass, levels = c("MFN", "GLP1", "DPP4", "SGLT2", "TZD", "SU"))) %>%
    rename("analysis_var" = "stopdrug_3m_6mFU") %>%
    table() %>%
    as.data.frame() %>%
    group_by(predrug_frailty_proxy, drugclass) %>%
    mutate(total = sum(Freq), Freq = Freq/total) %>%
    filter(analysis_var == 1) %>%
    ggplot(aes(y = predrug_frailty_proxy, x = Freq, fill = drugclass)) +
    geom_col(position = "dodge") +
    scale_fill_manual(values = c("Pooled" = "black", "SGLT2" = "#E69F00", "GLP1" = "#56B4E9", "SU" = "#CC79A7", "DPP4" = "#0072B2", "TZD" = "#D55E00", "MFN" = "grey"), breaks = c("MFN", "GLP1", "DPP4", "SGLT2", "TZD", "SU"), labels = c("Metformin", "GLP-1RA", "DPP4i", "SGLT2i", "TZD", "SU"), name = " ", guide = guide_legend(reverse = TRUE)) +
    scale_y_discrete(limits = rev, breaks = c(0, 1), labels = c("No", "Yes")) +
    ggtitle("Frailty") +
    scale_x_continuous("Proportion of discontinuation", labels = scales::percent, limits = c(0, 0.5), breaks = seq(0, 0.5, 0.1)) +
    facet_wrap(~drugclass, nrow = 1, labeller = as_labeller(c("MFN" = "Metformin", "GLP1" = "GLP-1RA", "DPP4" = "DPP4i", "SGLT2" = "SGLT2i", "TZD" = "TZD", "SU" = "SU"))) +
    theme_bw() +
    guides(fill = guide_legend(nrow = 1, byrow = TRUE)) +
    theme(
      legend.position = "bottom",
      axis.text.x = element_blank(),
      axis.title.x = element_blank(),
      strip.text.x = element_blank()
    ),
  
  # Ethnicity
  cprd_dataset %>%
    select(ethnicity_5cat, stopdrug_3m_6mFU, drugclass) %>%
    mutate(drugclass = factor(drugclass, levels = c("MFN", "GLP1", "DPP4", "SGLT2", "TZD", "SU"))) %>%
    mutate(ethnicity_5cat = factor(ethnicity_5cat, levels = c("4", "3", "2", "1", "0"), labels = c("Mixed & Other", "Mixed & Other", "Black", "South Asian", "White"))) %>%
    rename("analysis_var" = "stopdrug_3m_6mFU") %>%
    table() %>%
    as.data.frame() %>%
    group_by(ethnicity_5cat, drugclass) %>%
    mutate(total = sum(Freq), Freq = Freq/total) %>%
    filter(analysis_var == 1) %>%
    ggplot(aes(y = ethnicity_5cat, x = Freq, fill = drugclass)) +
    geom_col(position = "dodge") +
    scale_fill_manual(values = c("Pooled" = "black", "SGLT2" = "#E69F00", "GLP1" = "#56B4E9", "SU" = "#CC79A7", "DPP4" = "#0072B2", "TZD" = "#D55E00", "MFN" = "grey"), breaks = c("MFN", "GLP1", "DPP4", "SGLT2", "TZD", "SU"), labels = c("Metformin", "GLP-1RA", "DPP4i", "SGLT2i", "TZD", "SU"), name = " ", guide = guide_legend(reverse = TRUE)) +
    ggtitle("Ethnicity") +
    scale_x_continuous("Proportion of discontinuation", labels = scales::percent, limits = c(0, 0.5), breaks = seq(0, 0.5, 0.1)) +
    facet_wrap(~drugclass, nrow = 1, labeller = as_labeller(c("MFN" = "Metformin", "GLP1" = "GLP-1RA", "DPP4" = "DPP4i", "SGLT2" = "SGLT2i", "TZD" = "TZD", "SU" = "SU"))) +
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
    height = c(4, 2, 5, 4, 5, 4, 3, 2, 4),
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





# load dataset
cprd_dataset <- set_up_data(
  raw_data = "20240308_t2d_1stinstance",
  diagnosis = FALSE,
  therapies = c("DPP4", "GLP1", "MFN", "SGLT2", "SU", "TZD"),
  dataset = "6m.disc.dataset",
  full_prescribing_history = TRUE
) %>%
  drop_na(-stopdrug_12m_6mFU)



# Start plot by variables


pdf("results/figures/03.6m_discontinuation_by_variable.pdf", width = 7, height = 13)

patchwork::wrap_plots(
  
  # Overall
  cprd_dataset %>%
    select(stopdrug_6m_6mFU, drugclass) %>%
    mutate(drugclass = factor(drugclass, levels = c("MFN", "GLP1", "DPP4", "SGLT2", "TZD", "SU"), labels = c("MFN", "other", "other", "other", "other", "other"))) %>%
    rename("analysis_var" = "stopdrug_6m_6mFU") %>%
    table() %>%
    as.data.frame() %>%
    group_by(drugclass) %>%
    mutate(total = sum(Freq), Freq = Freq/total) %>%
    filter(analysis_var == 1) %>%
    ggplot(aes(y = analysis_var, x = Freq, fill = drugclass)) +
    geom_col() +
    scale_fill_manual(values = c("other" = "black", "SGLT2" = "#E69F00", "GLP1" = "#56B4E9", "SU" = "#CC79A7", "DPP4" = "#0072B2", "TZD" = "#D55E00", "MFN" = "grey"), breaks = c("MFN", "other"), labels = c("Metformin", "Other glucose-lowering therapies"), name = " ", guide = guide_legend(reverse = FALSE)) +
    scale_y_discrete(breaks = c(1), labels = c("Overall")) +
    scale_x_continuous("Proportion of discontinuation", labels = scales::percent, limits = c(0, 0.3)) +
    facet_wrap(~drugclass, nrow = 1, labeller = as_labeller(c("MFN" = "Metformin", "other" = "Other glucose-lowering therapies"))) +
    theme_bw() +
    guides(fill = guide_legend(nrow = 1, byrow = TRUE)) +
    theme(
      legend.position = "bottom",
      axis.text.x = element_blank(),
      axis.title.x = element_blank()
    ),
  
  # Age
  cprd_dataset %>%
    select(dstartdate_age, stopdrug_6m_6mFU, drugclass) %>%
    mutate(drugclass = factor(drugclass, levels = c("MFN", "GLP1", "DPP4", "SGLT2", "TZD", "SU"), labels = c("MFN", "other", "other", "other", "other", "other"))) %>%
    mutate(group = cut(dstartdate_age, breaks = c(18, 50, 60, 70, 150))) %>%
    rename("analysis_var" = "stopdrug_6m_6mFU") %>%
    select(-dstartdate_age) %>%
    table() %>%
    as.data.frame() %>%
    group_by(group, drugclass) %>%
    mutate(total = sum(Freq), Freq = Freq/total) %>%
    filter(analysis_var == 1) %>%
    ggplot(aes(y = group, x = Freq, fill = drugclass)) +
    geom_col() +
    scale_fill_manual(values = c("other" = "black", "SGLT2" = "#E69F00", "GLP1" = "#56B4E9", "SU" = "#CC79A7", "DPP4" = "#0072B2", "TZD" = "#D55E00", "MFN" = "grey"), breaks = c("MFN", "other"), labels = c("Metformin", "Other glucose-lowering therapies"), name = " ", guide = guide_legend(reverse = FALSE)) +
    scale_y_discrete(limits = rev, labels = c("70+", "60-69", "50-59", "<50")) +
    scale_x_continuous("Proportion of discontinuation", labels = scales::percent, limits = c(0, 0.3)) +
    ggtitle("Age group, years") +
    facet_wrap(~drugclass, nrow = 1, labeller = as_labeller(c("MFN" = "Metformin", "other" = "Other glucose-lowering therapies"))) +
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
  cprd_dataset %>%
    select(gender, stopdrug_6m_6mFU, drugclass) %>%
    mutate(drugclass = factor(drugclass, levels = c("MFN", "GLP1", "DPP4", "SGLT2", "TZD", "SU"), labels = c("MFN", "other", "other", "other", "other", "other"))) %>%
    rename("analysis_var" = "stopdrug_6m_6mFU") %>%
    table() %>%
    as.data.frame() %>%
    group_by(gender, drugclass) %>%
    mutate(total = sum(Freq), Freq = Freq/total) %>%
    filter(analysis_var == 1) %>%
    ggplot(aes(y = gender, x = Freq, fill = drugclass)) +
    geom_col() +
    scale_fill_manual(values = c("other" = "black", "SGLT2" = "#E69F00", "GLP1" = "#56B4E9", "SU" = "#CC79A7", "DPP4" = "#0072B2", "TZD" = "#D55E00", "MFN" = "grey"), breaks = c("MFN", "other"), labels = c("Metformin", "Other glucose-lowering therapies"), name = " ", guide = guide_legend(reverse = FALSE)) +
    scale_y_discrete(limits = rev, labels = c("Female", "Male")) +
    scale_x_continuous("Proportion of discontinuation", labels = scales::percent, limits = c(0, 0.3)) +
    ggtitle("Sex") +
    facet_wrap(~drugclass, nrow = 1, labeller = as_labeller(c("MFN" = "Metformin", "other" = "Other glucose-lowering therapies"))) +
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
  cprd_dataset %>%
    select(imd2015_10, stopdrug_6m_6mFU, drugclass) %>%
    mutate(drugclass = factor(drugclass, levels = c("MFN", "GLP1", "DPP4", "SGLT2", "TZD", "SU"), labels = c("MFN", "other", "other", "other", "other", "other"))) %>%
    rename("analysis_var" = "stopdrug_6m_6mFU") %>%
    table() %>%
    as.data.frame() %>%
    group_by(imd2015_10, drugclass) %>%
    mutate(total = sum(Freq), Freq = Freq/total) %>%
    filter(analysis_var == 1) %>%
    ggplot(aes(y = imd2015_10, x = Freq, fill = drugclass)) +
    geom_col() +
    scale_fill_manual(values = c("other" = "black", "SGLT2" = "#E69F00", "GLP1" = "#56B4E9", "SU" = "#CC79A7", "DPP4" = "#0072B2", "TZD" = "#D55E00", "MFN" = "grey"), breaks = c("MFN", "other"), labels = c("Metformin", "Other glucose-lowering therapies"), name = " ", guide = guide_legend(reverse = FALSE)) +
    scale_y_discrete(limits = rev, labels = c("5 (most deprived)", "4", "3", "2", "1 (least deprived)")) +
    scale_x_continuous("Proportion of discontinuation", labels = scales::percent, limits = c(0, 0.3)) +
    ggtitle("Index of multiple deprivation quintile") +
    facet_wrap(~drugclass, nrow = 1, labeller = as_labeller(c("MFN" = "Metformin", "other" = "Other glucose-lowering therapies"))) +
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
  cprd_dataset %>%
    select(prebmi, stopdrug_6m_6mFU, drugclass) %>%
    mutate(drugclass = factor(drugclass, levels = c("MFN", "GLP1", "DPP4", "SGLT2", "TZD", "SU"), labels = c("MFN", "other", "other", "other", "other", "other"))) %>%
    mutate(group = cut(prebmi, breaks = c(1, 25, 30, 35, 500))) %>%
    rename("analysis_var" = "stopdrug_6m_6mFU") %>%
    select(-prebmi) %>%
    table() %>%
    as.data.frame() %>%
    group_by(group, drugclass) %>%
    mutate(total = sum(Freq), Freq = Freq/total) %>%
    filter(analysis_var == 1) %>%
    ggplot(aes(y = group, x = Freq, fill = drugclass)) +
    geom_col() +
    scale_fill_manual(values = c("other" = "black", "SGLT2" = "#E69F00", "GLP1" = "#56B4E9", "SU" = "#CC79A7", "DPP4" = "#0072B2", "TZD" = "#D55E00", "MFN" = "grey"), breaks = c("MFN", "other"), labels = c("Metformin", "Other glucose-lowering therapies"), name = " ", guide = guide_legend(reverse = FALSE)) +
    scale_y_discrete(limits = rev, labels = c("35+", "30-34.9", "25-29.9", "<25")) +
    ggtitle("BMI, kg/m2") +
    scale_x_continuous("Proportion of discontinuation", labels = scales::percent, limits = c(0, 0.3)) +
    facet_wrap(~drugclass, nrow = 1, labeller = as_labeller(c("MFN" = "Metformin", "other" = "Other glucose-lowering therapies"))) +
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
  cprd_dataset %>%
    select(dstartdate_dm_dur, stopdrug_6m_6mFU, drugclass) %>%
    mutate(drugclass = factor(drugclass, levels = c("MFN", "GLP1", "DPP4", "SGLT2", "TZD", "SU"), labels = c("MFN", "other", "other", "other", "other", "other"))) %>%
    mutate(group = cut(dstartdate_dm_dur, breaks = c(-1, 1, 3, 6, 10, 150))) %>%
    rename("analysis_var" = "stopdrug_6m_6mFU") %>%
    select(-dstartdate_dm_dur) %>%
    table() %>%
    as.data.frame() %>%
    group_by(group, drugclass) %>%
    mutate(total = sum(Freq), Freq = Freq/total) %>%
    filter(analysis_var == 1) %>%
    ggplot(aes(y = group, x = Freq, fill = drugclass)) +
    geom_col() +
    scale_fill_manual(values = c("other" = "black", "SGLT2" = "#E69F00", "GLP1" = "#56B4E9", "SU" = "#CC79A7", "DPP4" = "#0072B2", "TZD" = "#D55E00", "MFN" = "grey"), breaks = c("MFN", "other"), labels = c("Metformin", "Other glucose-lowering therapies"), name = " ", guide = guide_legend(reverse = FALSE)) +
    scale_y_discrete(limits = rev, labels = c("10+", "6-9", "3-5", "1-2", "<1")) +
    ggtitle("Duration of diabetes, years") +
    scale_x_continuous("Proportion of discontinuation", labels = scales::percent, limits = c(0, 0.3)) +
    facet_wrap(~drugclass, nrow = 1, labeller = as_labeller(c("MFN" = "Metformin", "other" = "Other glucose-lowering therapies"))) +
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
  cprd_dataset %>%
    select(prehba1c, stopdrug_6m_6mFU, drugclass) %>%
    mutate(drugclass = factor(drugclass, levels = c("MFN", "GLP1", "DPP4", "SGLT2", "TZD", "SU"), labels = c("MFN", "other", "other", "other", "other", "other"))) %>%
    mutate(group = cut(prehba1c, breaks = c(52, 64, 75, 86, 250))) %>%
    rename("analysis_var" = "stopdrug_6m_6mFU") %>%
    select(-prehba1c) %>%
    table() %>%
    as.data.frame() %>%
    group_by(group, drugclass) %>%
    mutate(total = sum(Freq), Freq = Freq/total) %>%
    filter(analysis_var == 1) %>%
    ggplot(aes(y = group, x = Freq, fill = drugclass)) +
    geom_col() +
    scale_fill_manual(values = c("other" = "black", "SGLT2" = "#E69F00", "GLP1" = "#56B4E9", "SU" = "#CC79A7", "DPP4" = "#0072B2", "TZD" = "#D55E00", "MFN" = "grey"), breaks = c("MFN", "other"), labels = c("Metformin", "Other glucose-lowering therapies"), name = " ", guide = guide_legend(reverse = FALSE)) +
    scale_y_discrete(limits = rev, labels = c("86+", "75-86", "64-75", "53-64")) +
    ggtitle("HbA1c, mmol/mol") +
    scale_x_continuous("Proportion of discontinuation", labels = scales::percent, limits = c(0, 0.3)) +
    facet_wrap(~drugclass, nrow = 1, labeller = as_labeller(c("MFN" = "Metformin", "other" = "Other glucose-lowering therapies"))) +
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
  cprd_dataset %>%
    select(drugline, stopdrug_6m_6mFU, drugclass) %>%
    mutate(drugline = factor(drugline, levels = c("1", "2", "3", "4", "5+"), labels = c("1", "2", "3+", "3+", "3+"))) %>%
    mutate(drugclass = factor(drugclass, levels = c("MFN", "GLP1", "DPP4", "SGLT2", "TZD", "SU"), labels = c("MFN", "other", "other", "other", "other", "other"))) %>%
    rename("analysis_var" = "stopdrug_6m_6mFU") %>%
    table() %>%
    as.data.frame() %>%
    group_by(drugline, drugclass) %>%
    mutate(total = sum(Freq), Freq = Freq/total) %>%
    filter(analysis_var == 1) %>%
    ggplot(aes(y = drugline, x = Freq, fill = drugclass)) +
    geom_col() +
    scale_fill_manual(values = c("other" = "black", "SGLT2" = "#E69F00", "GLP1" = "#56B4E9", "SU" = "#CC79A7", "DPP4" = "#0072B2", "TZD" = "#D55E00", "MFN" = "grey"), breaks = c("MFN", "other"), labels = c("Metformin", "Other glucose-lowering therapies"), name = " ", guide = guide_legend(reverse = FALSE)) +
    scale_y_discrete(limits = rev) +
    ggtitle("Number of therapy classes ever prescribed") +
    scale_x_continuous("Proportion of discontinuation", labels = scales::percent, limits = c(0, 0.3)) +
    facet_wrap(~drugclass, nrow = 1, labeller = as_labeller(c("MFN" = "Metformin", "other" = "Other glucose-lowering therapies"))) +
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
  cprd_dataset %>%
    select(predrug_frailty_proxy, stopdrug_6m_6mFU, drugclass) %>%
    mutate(drugclass = factor(drugclass, levels = c("MFN", "GLP1", "DPP4", "SGLT2", "TZD", "SU"), labels = c("MFN", "other", "other", "other", "other", "other"))) %>%
    rename("analysis_var" = "stopdrug_6m_6mFU") %>%
    table() %>%
    as.data.frame() %>%
    group_by(predrug_frailty_proxy, drugclass) %>%
    mutate(total = sum(Freq), Freq = Freq/total) %>%
    filter(analysis_var == 1) %>%
    ggplot(aes(y = predrug_frailty_proxy, x = Freq, fill = drugclass)) +
    geom_col() +
    scale_fill_manual(values = c("other" = "black", "SGLT2" = "#E69F00", "GLP1" = "#56B4E9", "SU" = "#CC79A7", "DPP4" = "#0072B2", "TZD" = "#D55E00", "MFN" = "grey"), breaks = c("MFN", "other"), labels = c("Metformin", "Other glucose-lowering therapies"), name = " ", guide = guide_legend(reverse = FALSE)) +
    scale_y_discrete(limits = rev, breaks = c(0, 1), labels = c("No", "Yes")) +
    ggtitle("Frailty") +
    scale_x_continuous("Proportion of discontinuation", labels = scales::percent, limits = c(0, 0.3)) +
    facet_wrap(~drugclass, nrow = 1, labeller = as_labeller(c("MFN" = "Metformin", "other" = "Other glucose-lowering therapies"))) +
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
  cprd_dataset %>%
    select(ethnicity_5cat, stopdrug_6m_6mFU, drugclass) %>%
    mutate(ethnicity_5cat = factor(ethnicity_5cat, levels = c("4", "3", "2", "1", "0"), labels = c("Mixed & Other", "Mixed & Other", "Black", "South Asian", "White"))) %>%
    mutate(drugclass = factor(drugclass, levels = c("MFN", "GLP1", "DPP4", "SGLT2", "TZD", "SU"), labels = c("MFN", "other", "other", "other", "other", "other"))) %>%
    rename("analysis_var" = "stopdrug_6m_6mFU") %>%
    table() %>%
    as.data.frame() %>%
    group_by(ethnicity_5cat, drugclass) %>%
    mutate(total = sum(Freq), Freq = Freq/total) %>%
    filter(analysis_var == 1) %>%
    ggplot(aes(y = ethnicity_5cat, x = Freq, fill = drugclass)) +
    geom_col() +
    scale_fill_manual(values = c("other" = "black", "SGLT2" = "#E69F00", "GLP1" = "#56B4E9", "SU" = "#CC79A7", "DPP4" = "#0072B2", "TZD" = "#D55E00", "MFN" = "grey"), breaks = c("MFN", "other"), labels = c("Metformin", "Other glucose-lowering therapies"), name = " ", guide = guide_legend(reverse = FALSE)) +
    ggtitle("Ethnicity") +
    scale_x_continuous("Proportion of discontinuation", labels = scales::percent, limits = c(0, 0.3)) +
    facet_wrap(~drugclass, nrow = 1, labeller = as_labeller(c("MFN" = "Metformin", "other" = "Other glucose-lowering therapies"))) +
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
    height = c(1, 4, 2, 5, 4, 5, 4, 3, 2, 4),
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



pdf("results/figures/03.6m_discontinuation_by_variable_by_drug.pdf", width = 14, height = 13)

patchwork::wrap_plots(
  
  # Age - Metformin
  cprd_dataset %>%
    select(dstartdate_age, stopdrug_6m_6mFU, drugclass) %>%
    mutate(drugclass = factor(drugclass, levels = c("MFN", "GLP1", "DPP4", "SGLT2", "TZD", "SU"))) %>%
    mutate(group = cut(dstartdate_age, breaks = c(18, 50, 60, 70, 150))) %>%
    rename("analysis_var" = "stopdrug_6m_6mFU") %>%
    select(-dstartdate_age) %>%
    table() %>%
    as.data.frame() %>%
    group_by(group, drugclass) %>%
    mutate(total = sum(Freq), Freq = Freq/total) %>%
    filter(analysis_var == 1) %>%
    ggplot(aes(y = group, x = Freq, fill = drugclass)) +
    geom_col(position = "dodge") +
    scale_fill_manual(values = c("Pooled" = "black", "SGLT2" = "#E69F00", "GLP1" = "#56B4E9", "SU" = "#CC79A7", "DPP4" = "#0072B2", "TZD" = "#D55E00", "MFN" = "grey"), breaks = c("MFN", "GLP1", "DPP4", "SGLT2", "TZD", "SU"), labels = c("Metformin", "GLP-1RA", "DPP4i", "SGLT2i", "TZD", "SU"), name = " ", guide = guide_legend(reverse = TRUE)) +
    scale_y_discrete(limits = rev, labels = c("70+", "60-69", "50-59", "<50")) +
    ggtitle("Age group, years") +
    scale_x_continuous("Proportion of discontinuation", labels = scales::percent, limits = c(0, 0.5), breaks = seq(0, 0.5, 0.1)) +
    facet_wrap(~drugclass, nrow = 1, labeller = as_labeller(c("MFN" = "Metformin", "GLP1" = "GLP-1RA", "DPP4" = "DPP4i", "SGLT2" = "SGLT2i", "TZD" = "TZD", "SU" = "SU"))) +
    theme_bw() +
    guides(fill = guide_legend(nrow = 1, byrow = TRUE)) +
    theme(
      legend.position = "bottom",
      axis.text.x = element_blank(),
      axis.title.x = element_blank()
    ),
  
  # Sex 
  cprd_dataset %>%
    select(gender, stopdrug_6m_6mFU, drugclass) %>%
    mutate(drugclass = factor(drugclass, levels = c("MFN", "GLP1", "DPP4", "SGLT2", "TZD", "SU"))) %>%
    rename("analysis_var" = "stopdrug_6m_6mFU") %>%
    table() %>%
    as.data.frame() %>%
    group_by(gender, drugclass) %>%
    mutate(total = sum(Freq), Freq = Freq/total) %>%
    filter(analysis_var == 1) %>%
    ggplot(aes(y = gender, x = Freq, fill = drugclass)) +
    geom_col(position = "dodge") +
    scale_fill_manual(values = c("Pooled" = "black", "SGLT2" = "#E69F00", "GLP1" = "#56B4E9", "SU" = "#CC79A7", "DPP4" = "#0072B2", "TZD" = "#D55E00", "MFN" = "grey"), breaks = c("MFN", "GLP1", "DPP4", "SGLT2", "TZD", "SU"), labels = c("Metformin", "GLP-1RA", "DPP4i", "SGLT2i", "TZD", "SU"), name = " ", guide = guide_legend(reverse = TRUE)) +
    scale_y_discrete(limits = rev, labels = c("Female", "Male")) +
    ggtitle("Sex") +
    scale_x_continuous("Proportion of discontinuation", labels = scales::percent, limits = c(0, 0.5), breaks = seq(0, 0.5, 0.1)) +
    facet_wrap(~drugclass, nrow = 1, labeller = as_labeller(c("MFN" = "Metformin", "GLP1" = "GLP-1RA", "DPP4" = "DPP4i", "SGLT2" = "SGLT2i", "TZD" = "TZD", "SU" = "SU"))) +
    theme_bw() +
    guides(fill = guide_legend(nrow = 1, byrow = TRUE)) +
    theme(
      legend.position = "bottom",
      axis.text.x = element_blank(),
      axis.title.x = element_blank(),
      strip.text.x = element_blank()
    ),
  
  # Quintile
  cprd_dataset %>%
    select(imd2015_10, stopdrug_6m_6mFU, drugclass) %>%
    mutate(drugclass = factor(drugclass, levels = c("MFN", "GLP1", "DPP4", "SGLT2", "TZD", "SU"))) %>%
    rename("analysis_var" = "stopdrug_6m_6mFU") %>%
    table() %>%
    as.data.frame() %>%
    group_by(imd2015_10, drugclass) %>%
    mutate(total = sum(Freq), Freq = Freq/total) %>%
    filter(analysis_var == 1) %>%
    ggplot(aes(y = imd2015_10, x = Freq, fill = drugclass)) +
    geom_col(position = "dodge") +
    scale_fill_manual(values = c("Pooled" = "black", "SGLT2" = "#E69F00", "GLP1" = "#56B4E9", "SU" = "#CC79A7", "DPP4" = "#0072B2", "TZD" = "#D55E00", "MFN" = "grey"), breaks = c("MFN", "GLP1", "DPP4", "SGLT2", "TZD", "SU"), labels = c("Metformin", "GLP-1RA", "DPP4i", "SGLT2i", "TZD", "SU"), name = " ", guide = guide_legend(reverse = TRUE)) +
    scale_y_discrete(limits = rev, labels = c("5 (most deprived)", "4", "3", "2", "1 (least deprived)")) +
    ggtitle("Index of multiple deprivation quintile") +
    scale_x_continuous("Proportion of discontinuation", labels = scales::percent, limits = c(0, 0.5), breaks = seq(0, 0.5, 0.1)) +
    facet_wrap(~drugclass, nrow = 1, labeller = as_labeller(c("MFN" = "Metformin", "GLP1" = "GLP-1RA", "DPP4" = "DPP4i", "SGLT2" = "SGLT2i", "TZD" = "TZD", "SU" = "SU"))) +
    theme_bw() +
    guides(fill = guide_legend(nrow = 1, byrow = TRUE)) +
    theme(
      legend.position = "bottom",
      axis.text.x = element_blank(),
      axis.title.x = element_blank(),
      strip.text.x = element_blank()
    ),
  
  # BMI
  cprd_dataset %>%
    select(prebmi, stopdrug_6m_6mFU, drugclass) %>%
    mutate(drugclass = factor(drugclass, levels = c("MFN", "GLP1", "DPP4", "SGLT2", "TZD", "SU"))) %>%
    mutate(group = cut(prebmi, breaks = c(1, 25, 30, 35, 500))) %>%
    rename("analysis_var" = "stopdrug_6m_6mFU") %>%
    select(-prebmi) %>%
    table() %>%
    as.data.frame() %>%
    group_by(group, drugclass) %>%
    mutate(total = sum(Freq), Freq = Freq/total) %>%
    filter(analysis_var == 1) %>%
    ggplot(aes(y = group, x = Freq, fill = drugclass)) +
    geom_col(position = "dodge") +
    scale_fill_manual(values = c("Pooled" = "black", "SGLT2" = "#E69F00", "GLP1" = "#56B4E9", "SU" = "#CC79A7", "DPP4" = "#0072B2", "TZD" = "#D55E00", "MFN" = "grey"), breaks = c("MFN", "GLP1", "DPP4", "SGLT2", "TZD", "SU"), labels = c("Metformin", "GLP-1RA", "DPP4i", "SGLT2i", "TZD", "SU"), name = " ", guide = guide_legend(reverse = TRUE)) +
    scale_y_discrete(limits = rev, labels = c("35+", "30-34.9", "25-29.9", "<25")) +
    ggtitle("BMI, kg/m2") +
    scale_x_continuous("Proportion of discontinuation", labels = scales::percent, limits = c(0, 0.5), breaks = seq(0, 0.5, 0.1)) +
    facet_wrap(~drugclass, nrow = 1, labeller = as_labeller(c("MFN" = "Metformin", "GLP1" = "GLP-1RA", "DPP4" = "DPP4i", "SGLT2" = "SGLT2i", "TZD" = "TZD", "SU" = "SU"))) +
    theme_bw() +
    guides(fill = guide_legend(nrow = 1, byrow = TRUE)) +
    theme(
      legend.position = "bottom",
      axis.text.x = element_blank(),
      axis.title.x = element_blank(),
      strip.text.x = element_blank()
    ),
  
  # Diabetes duration
  cprd_dataset %>%
    select(dstartdate_dm_dur, stopdrug_6m_6mFU, drugclass) %>%
    mutate(drugclass = factor(drugclass, levels = c("MFN", "GLP1", "DPP4", "SGLT2", "TZD", "SU"))) %>%
    mutate(group = cut(dstartdate_dm_dur, breaks = c(-1, 1, 3, 6, 10, 150))) %>%
    rename("analysis_var" = "stopdrug_6m_6mFU") %>%
    select(-dstartdate_dm_dur) %>%
    table() %>%
    as.data.frame() %>%
    group_by(group, drugclass) %>%
    mutate(total = sum(Freq), Freq = Freq/total) %>%
    filter(analysis_var == 1) %>%
    ggplot(aes(y = group, x = Freq, fill = drugclass)) +
    geom_col(position = "dodge") +
    scale_fill_manual(values = c("Pooled" = "black", "SGLT2" = "#E69F00", "GLP1" = "#56B4E9", "SU" = "#CC79A7", "DPP4" = "#0072B2", "TZD" = "#D55E00", "MFN" = "grey"), breaks = c("MFN", "GLP1", "DPP4", "SGLT2", "TZD", "SU"), labels = c("Metformin", "GLP-1RA", "DPP4i", "SGLT2i", "TZD", "SU"), name = " ", guide = guide_legend(reverse = TRUE)) +
    scale_y_discrete(limits = rev, labels = c("10+", "6-9", "3-5", "1-2", "<1")) +
    ggtitle("Duration of diabetes, years") +
    scale_x_continuous("Proportion of discontinuation", labels = scales::percent, limits = c(0, 0.5), breaks = seq(0, 0.5, 0.1)) +
    facet_wrap(~drugclass, nrow = 1, labeller = as_labeller(c("MFN" = "Metformin", "GLP1" = "GLP-1RA", "DPP4" = "DPP4i", "SGLT2" = "SGLT2i", "TZD" = "TZD", "SU" = "SU"))) +
    theme_bw() +
    guides(fill = guide_legend(nrow = 1, byrow = TRUE)) +
    theme(
      legend.position = "bottom",
      axis.text.x = element_blank(),
      axis.title.x = element_blank(),
      strip.text.x = element_blank()
    ),
  
  # HbA1c
  cprd_dataset %>%
    select(prehba1c, stopdrug_6m_6mFU, drugclass) %>%
    mutate(drugclass = factor(drugclass, levels = c("MFN", "GLP1", "DPP4", "SGLT2", "TZD", "SU"))) %>%
    mutate(group = cut(prehba1c, breaks = c(52, 64, 75, 86, 250))) %>%
    rename("analysis_var" = "stopdrug_6m_6mFU") %>%
    select(-prehba1c) %>%
    table() %>%
    as.data.frame() %>%
    group_by(group, drugclass) %>%
    mutate(total = sum(Freq), Freq = Freq/total) %>%
    filter(analysis_var == 1) %>%
    ggplot(aes(y = group, x = Freq, fill = drugclass)) +
    geom_col(position = "dodge") +
    scale_fill_manual(values = c("Pooled" = "black", "SGLT2" = "#E69F00", "GLP1" = "#56B4E9", "SU" = "#CC79A7", "DPP4" = "#0072B2", "TZD" = "#D55E00", "MFN" = "grey"), breaks = c("MFN", "GLP1", "DPP4", "SGLT2", "TZD", "SU"), labels = c("Metformin", "GLP-1RA", "DPP4i", "SGLT2i", "TZD", "SU"), name = " ", guide = guide_legend(reverse = TRUE)) +
    scale_y_discrete(limits = rev, labels = c("86+", "75-86", "64-75", "53-64")) +
    ggtitle("HbA1c, mmol/mol") +
    scale_x_continuous("Proportion of discontinuation", labels = scales::percent, limits = c(0, 0.5), breaks = seq(0, 0.5, 0.1)) +
    facet_wrap(~drugclass, nrow = 1, labeller = as_labeller(c("MFN" = "Metformin", "GLP1" = "GLP-1RA", "DPP4" = "DPP4i", "SGLT2" = "SGLT2i", "TZD" = "TZD", "SU" = "SU"))) +
    theme_bw() +
    guides(fill = guide_legend(nrow = 1, byrow = TRUE)) +
    theme(
      legend.position = "bottom",
      axis.text.x = element_blank(),
      axis.title.x = element_blank(),
      strip.text.x = element_blank()
    ),
  
  # drugline
  cprd_dataset %>%
    select(drugline, stopdrug_6m_6mFU, drugclass) %>%
    mutate(drugclass = factor(drugclass, levels = c("MFN", "GLP1", "DPP4", "SGLT2", "TZD", "SU"))) %>%
    mutate(drugline = factor(drugline, levels = c("1", "2", "3", "4", "5+"), labels = c("1", "2", "3+", "3+", "3+"))) %>%
    rename("analysis_var" = "stopdrug_6m_6mFU") %>%
    table() %>%
    as.data.frame() %>%
    group_by(drugline, drugclass) %>%
    mutate(total = sum(Freq), Freq = Freq/total) %>%
    filter(analysis_var == 1) %>%
    ggplot(aes(y = drugline, x = Freq, fill = drugclass)) +
    geom_col(position = "dodge") +
    scale_fill_manual(values = c("Pooled" = "black", "SGLT2" = "#E69F00", "GLP1" = "#56B4E9", "SU" = "#CC79A7", "DPP4" = "#0072B2", "TZD" = "#D55E00", "MFN" = "grey"), breaks = c("MFN", "GLP1", "DPP4", "SGLT2", "TZD", "SU"), labels = c("Metformin", "GLP-1RA", "DPP4i", "SGLT2i", "TZD", "SU"), name = " ", guide = guide_legend(reverse = TRUE)) +
    scale_y_discrete(limits = rev) +
    ggtitle("Number of therapy classes ever prescribed") +
    scale_x_continuous("Proportion of discontinuation", labels = scales::percent, limits = c(0, 0.5), breaks = seq(0, 0.5, 0.1)) +
    facet_wrap(~drugclass, nrow = 1, labeller = as_labeller(c("MFN" = "Metformin", "GLP1" = "GLP-1RA", "DPP4" = "DPP4i", "SGLT2" = "SGLT2i", "TZD" = "TZD", "SU" = "SU"))) +
    theme_bw() +
    guides(fill = guide_legend(nrow = 1, byrow = TRUE)) +
    theme(
      legend.position = "bottom",
      axis.text.x = element_blank(),
      axis.title.x = element_blank(),
      strip.text.x = element_blank()
    ),
  
  # Frailty
  cprd_dataset %>%
    select(predrug_frailty_proxy, stopdrug_6m_6mFU, drugclass) %>%
    mutate(drugclass = factor(drugclass, levels = c("MFN", "GLP1", "DPP4", "SGLT2", "TZD", "SU"))) %>%
    rename("analysis_var" = "stopdrug_6m_6mFU") %>%
    table() %>%
    as.data.frame() %>%
    group_by(predrug_frailty_proxy, drugclass) %>%
    mutate(total = sum(Freq), Freq = Freq/total) %>%
    filter(analysis_var == 1) %>%
    ggplot(aes(y = predrug_frailty_proxy, x = Freq, fill = drugclass)) +
    geom_col(position = "dodge") +
    scale_fill_manual(values = c("Pooled" = "black", "SGLT2" = "#E69F00", "GLP1" = "#56B4E9", "SU" = "#CC79A7", "DPP4" = "#0072B2", "TZD" = "#D55E00", "MFN" = "grey"), breaks = c("MFN", "GLP1", "DPP4", "SGLT2", "TZD", "SU"), labels = c("Metformin", "GLP-1RA", "DPP4i", "SGLT2i", "TZD", "SU"), name = " ", guide = guide_legend(reverse = TRUE)) +
    scale_y_discrete(limits = rev, breaks = c(0, 1), labels = c("No", "Yes")) +
    ggtitle("Frailty") +
    scale_x_continuous("Proportion of discontinuation", labels = scales::percent, limits = c(0, 0.5), breaks = seq(0, 0.5, 0.1)) +
    facet_wrap(~drugclass, nrow = 1, labeller = as_labeller(c("MFN" = "Metformin", "GLP1" = "GLP-1RA", "DPP4" = "DPP4i", "SGLT2" = "SGLT2i", "TZD" = "TZD", "SU" = "SU"))) +
    theme_bw() +
    guides(fill = guide_legend(nrow = 1, byrow = TRUE)) +
    theme(
      legend.position = "bottom",
      axis.text.x = element_blank(),
      axis.title.x = element_blank(),
      strip.text.x = element_blank()
    ),
  
  # Ethnicity
  cprd_dataset %>%
    select(ethnicity_5cat, stopdrug_6m_6mFU, drugclass) %>%
    mutate(drugclass = factor(drugclass, levels = c("MFN", "GLP1", "DPP4", "SGLT2", "TZD", "SU"))) %>%
    mutate(ethnicity_5cat = factor(ethnicity_5cat, levels = c("4", "3", "2", "1", "0"), labels = c("Mixed & Other", "Mixed & Other", "Black", "South Asian", "White"))) %>%
    rename("analysis_var" = "stopdrug_6m_6mFU") %>%
    table() %>%
    as.data.frame() %>%
    group_by(ethnicity_5cat, drugclass) %>%
    mutate(total = sum(Freq), Freq = Freq/total) %>%
    filter(analysis_var == 1) %>%
    ggplot(aes(y = ethnicity_5cat, x = Freq, fill = drugclass)) +
    geom_col(position = "dodge") +
    scale_fill_manual(values = c("Pooled" = "black", "SGLT2" = "#E69F00", "GLP1" = "#56B4E9", "SU" = "#CC79A7", "DPP4" = "#0072B2", "TZD" = "#D55E00", "MFN" = "grey"), breaks = c("MFN", "GLP1", "DPP4", "SGLT2", "TZD", "SU"), labels = c("Metformin", "GLP-1RA", "DPP4i", "SGLT2i", "TZD", "SU"), name = " ", guide = guide_legend(reverse = TRUE)) +
    ggtitle("Ethnicity") +
    scale_x_continuous("Proportion of discontinuation", labels = scales::percent, limits = c(0, 0.5), breaks = seq(0, 0.5, 0.1)) +
    facet_wrap(~drugclass, nrow = 1, labeller = as_labeller(c("MFN" = "Metformin", "GLP1" = "GLP-1RA", "DPP4" = "DPP4i", "SGLT2" = "SGLT2i", "TZD" = "TZD", "SU" = "SU"))) +
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
    height = c(4, 2, 5, 4, 5, 4, 3, 2, 4),
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


# load dataset
cprd_dataset <- set_up_data(
  raw_data = "20240308_t2d_1stinstance",
  diagnosis = FALSE,
  therapies = c("DPP4", "GLP1", "MFN", "SGLT2", "SU", "TZD"),
  dataset = "6m.disc.dataset",
  full_prescribing_history = TRUE
) %>%
  drop_na()




# Start plot by variables


pdf("results/figures/03.12m_discontinuation_by_variable.pdf", width = 7, height = 13)

patchwork::wrap_plots(
  
  # Overall
  cprd_dataset %>%
    select(stopdrug_12m_6mFU, drugclass) %>%
    mutate(drugclass = factor(drugclass, levels = c("MFN", "GLP1", "DPP4", "SGLT2", "TZD", "SU"), labels = c("MFN", "other", "other", "other", "other", "other"))) %>%
    rename("analysis_var" = "stopdrug_12m_6mFU") %>%
    table() %>%
    as.data.frame() %>%
    group_by(drugclass) %>%
    mutate(total = sum(Freq), Freq = Freq/total) %>%
    filter(analysis_var == 1) %>%
    ggplot(aes(y = analysis_var, x = Freq, fill = drugclass)) +
    geom_col() +
    scale_fill_manual(values = c("other" = "black", "SGLT2" = "#E69F00", "GLP1" = "#56B4E9", "SU" = "#CC79A7", "DPP4" = "#0072B2", "TZD" = "#D55E00", "MFN" = "grey"), breaks = c("MFN", "other"), labels = c("Metformin", "Other glucose-lowering therapies"), name = " ", guide = guide_legend(reverse = FALSE)) +
    scale_y_discrete(breaks = c(1), labels = c("Overall")) +
    scale_x_continuous("Proportion of discontinuation", labels = scales::percent, limits = c(0, 0.3)) +
    facet_wrap(~drugclass, nrow = 1, labeller = as_labeller(c("MFN" = "Metformin", "other" = "Other glucose-lowering therapies"))) +
    theme_bw() +
    guides(fill = guide_legend(nrow = 1, byrow = TRUE)) +
    theme(
      legend.position = "bottom",
      axis.text.x = element_blank(),
      axis.title.x = element_blank()
    ),
  
  # Age
  cprd_dataset %>%
    select(dstartdate_age, stopdrug_12m_6mFU, drugclass) %>%
    mutate(drugclass = factor(drugclass, levels = c("MFN", "GLP1", "DPP4", "SGLT2", "TZD", "SU"), labels = c("MFN", "other", "other", "other", "other", "other"))) %>%
    mutate(group = cut(dstartdate_age, breaks = c(18, 50, 60, 70, 150))) %>%
    rename("analysis_var" = "stopdrug_12m_6mFU") %>%
    select(-dstartdate_age) %>%
    table() %>%
    as.data.frame() %>%
    group_by(group, drugclass) %>%
    mutate(total = sum(Freq), Freq = Freq/total) %>%
    filter(analysis_var == 1) %>%
    ggplot(aes(y = group, x = Freq, fill = drugclass)) +
    geom_col() +
    scale_fill_manual(values = c("other" = "black", "SGLT2" = "#E69F00", "GLP1" = "#56B4E9", "SU" = "#CC79A7", "DPP4" = "#0072B2", "TZD" = "#D55E00", "MFN" = "grey"), breaks = c("MFN", "other"), labels = c("Metformin", "Other glucose-lowering therapies"), name = " ", guide = guide_legend(reverse = FALSE)) +
    scale_y_discrete(limits = rev, labels = c("70+", "60-69", "50-59", "<50")) +
    scale_x_continuous("Proportion of discontinuation", labels = scales::percent, limits = c(0, 0.3)) +
    ggtitle("Age group, years") +
    facet_wrap(~drugclass, nrow = 1, labeller = as_labeller(c("MFN" = "Metformin", "other" = "Other glucose-lowering therapies"))) +
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
  cprd_dataset %>%
    select(gender, stopdrug_12m_6mFU, drugclass) %>%
    mutate(drugclass = factor(drugclass, levels = c("MFN", "GLP1", "DPP4", "SGLT2", "TZD", "SU"), labels = c("MFN", "other", "other", "other", "other", "other"))) %>%
    rename("analysis_var" = "stopdrug_12m_6mFU") %>%
    table() %>%
    as.data.frame() %>%
    group_by(gender, drugclass) %>%
    mutate(total = sum(Freq), Freq = Freq/total) %>%
    filter(analysis_var == 1) %>%
    ggplot(aes(y = gender, x = Freq, fill = drugclass)) +
    geom_col() +
    scale_fill_manual(values = c("other" = "black", "SGLT2" = "#E69F00", "GLP1" = "#56B4E9", "SU" = "#CC79A7", "DPP4" = "#0072B2", "TZD" = "#D55E00", "MFN" = "grey"), breaks = c("MFN", "other"), labels = c("Metformin", "Other glucose-lowering therapies"), name = " ", guide = guide_legend(reverse = FALSE)) +
    scale_y_discrete(limits = rev, labels = c("Female", "Male")) +
    scale_x_continuous("Proportion of discontinuation", labels = scales::percent, limits = c(0, 0.3)) +
    ggtitle("Sex") +
    facet_wrap(~drugclass, nrow = 1, labeller = as_labeller(c("MFN" = "Metformin", "other" = "Other glucose-lowering therapies"))) +
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
  cprd_dataset %>%
    select(imd2015_10, stopdrug_12m_6mFU, drugclass) %>%
    mutate(drugclass = factor(drugclass, levels = c("MFN", "GLP1", "DPP4", "SGLT2", "TZD", "SU"), labels = c("MFN", "other", "other", "other", "other", "other"))) %>%
    rename("analysis_var" = "stopdrug_12m_6mFU") %>%
    table() %>%
    as.data.frame() %>%
    group_by(imd2015_10, drugclass) %>%
    mutate(total = sum(Freq), Freq = Freq/total) %>%
    filter(analysis_var == 1) %>%
    ggplot(aes(y = imd2015_10, x = Freq, fill = drugclass)) +
    geom_col() +
    scale_fill_manual(values = c("other" = "black", "SGLT2" = "#E69F00", "GLP1" = "#56B4E9", "SU" = "#CC79A7", "DPP4" = "#0072B2", "TZD" = "#D55E00", "MFN" = "grey"), breaks = c("MFN", "other"), labels = c("Metformin", "Other glucose-lowering therapies"), name = " ", guide = guide_legend(reverse = FALSE)) +
    scale_y_discrete(limits = rev, labels = c("5 (most deprived)", "4", "3", "2", "1 (least deprived)")) +
    scale_x_continuous("Proportion of discontinuation", labels = scales::percent, limits = c(0, 0.3)) +
    ggtitle("Index of multiple deprivation quintile") +
    facet_wrap(~drugclass, nrow = 1, labeller = as_labeller(c("MFN" = "Metformin", "other" = "Other glucose-lowering therapies"))) +
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
  cprd_dataset %>%
    select(prebmi, stopdrug_12m_6mFU, drugclass) %>%
    mutate(drugclass = factor(drugclass, levels = c("MFN", "GLP1", "DPP4", "SGLT2", "TZD", "SU"), labels = c("MFN", "other", "other", "other", "other", "other"))) %>%
    mutate(group = cut(prebmi, breaks = c(1, 25, 30, 35, 500))) %>%
    rename("analysis_var" = "stopdrug_12m_6mFU") %>%
    select(-prebmi) %>%
    table() %>%
    as.data.frame() %>%
    group_by(group, drugclass) %>%
    mutate(total = sum(Freq), Freq = Freq/total) %>%
    filter(analysis_var == 1) %>%
    ggplot(aes(y = group, x = Freq, fill = drugclass)) +
    geom_col() +
    scale_fill_manual(values = c("other" = "black", "SGLT2" = "#E69F00", "GLP1" = "#56B4E9", "SU" = "#CC79A7", "DPP4" = "#0072B2", "TZD" = "#D55E00", "MFN" = "grey"), breaks = c("MFN", "other"), labels = c("Metformin", "Other glucose-lowering therapies"), name = " ", guide = guide_legend(reverse = FALSE)) +
    scale_y_discrete(limits = rev, labels = c("35+", "30-34.9", "25-29.9", "<25")) +
    ggtitle("BMI, kg/m2") +
    scale_x_continuous("Proportion of discontinuation", labels = scales::percent, limits = c(0, 0.3)) +
    facet_wrap(~drugclass, nrow = 1, labeller = as_labeller(c("MFN" = "Metformin", "other" = "Other glucose-lowering therapies"))) +
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
  cprd_dataset %>%
    select(dstartdate_dm_dur, stopdrug_12m_6mFU, drugclass) %>%
    mutate(drugclass = factor(drugclass, levels = c("MFN", "GLP1", "DPP4", "SGLT2", "TZD", "SU"), labels = c("MFN", "other", "other", "other", "other", "other"))) %>%
    mutate(group = cut(dstartdate_dm_dur, breaks = c(-1, 1, 3, 6, 10, 150))) %>%
    rename("analysis_var" = "stopdrug_12m_6mFU") %>%
    select(-dstartdate_dm_dur) %>%
    table() %>%
    as.data.frame() %>%
    group_by(group, drugclass) %>%
    mutate(total = sum(Freq), Freq = Freq/total) %>%
    filter(analysis_var == 1) %>%
    ggplot(aes(y = group, x = Freq, fill = drugclass)) +
    geom_col() +
    scale_fill_manual(values = c("other" = "black", "SGLT2" = "#E69F00", "GLP1" = "#56B4E9", "SU" = "#CC79A7", "DPP4" = "#0072B2", "TZD" = "#D55E00", "MFN" = "grey"), breaks = c("MFN", "other"), labels = c("Metformin", "Other glucose-lowering therapies"), name = " ", guide = guide_legend(reverse = FALSE)) +
    scale_y_discrete(limits = rev, labels = c("10+", "6-9", "3-5", "1-2", "<1")) +
    ggtitle("Duration of diabetes, years") +
    scale_x_continuous("Proportion of discontinuation", labels = scales::percent, limits = c(0, 0.3)) +
    facet_wrap(~drugclass, nrow = 1, labeller = as_labeller(c("MFN" = "Metformin", "other" = "Other glucose-lowering therapies"))) +
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
  cprd_dataset %>%
    select(prehba1c, stopdrug_12m_6mFU, drugclass) %>%
    mutate(drugclass = factor(drugclass, levels = c("MFN", "GLP1", "DPP4", "SGLT2", "TZD", "SU"), labels = c("MFN", "other", "other", "other", "other", "other"))) %>%
    mutate(group = cut(prehba1c, breaks = c(52, 64, 75, 86, 250))) %>%
    rename("analysis_var" = "stopdrug_12m_6mFU") %>%
    select(-prehba1c) %>%
    table() %>%
    as.data.frame() %>%
    group_by(group, drugclass) %>%
    mutate(total = sum(Freq), Freq = Freq/total) %>%
    filter(analysis_var == 1) %>%
    ggplot(aes(y = group, x = Freq, fill = drugclass)) +
    geom_col() +
    scale_fill_manual(values = c("other" = "black", "SGLT2" = "#E69F00", "GLP1" = "#56B4E9", "SU" = "#CC79A7", "DPP4" = "#0072B2", "TZD" = "#D55E00", "MFN" = "grey"), breaks = c("MFN", "other"), labels = c("Metformin", "Other glucose-lowering therapies"), name = " ", guide = guide_legend(reverse = FALSE)) +
    scale_y_discrete(limits = rev, labels = c("86+", "75-86", "64-75", "53-64")) +
    ggtitle("HbA1c, mmol/mol") +
    scale_x_continuous("Proportion of discontinuation", labels = scales::percent, limits = c(0, 0.3)) +
    facet_wrap(~drugclass, nrow = 1, labeller = as_labeller(c("MFN" = "Metformin", "other" = "Other glucose-lowering therapies"))) +
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
  cprd_dataset %>%
    select(drugline, stopdrug_12m_6mFU, drugclass) %>%
    mutate(drugline = factor(drugline, levels = c("1", "2", "3", "4", "5+"), labels = c("1", "2", "3+", "3+", "3+"))) %>%
    mutate(drugclass = factor(drugclass, levels = c("MFN", "GLP1", "DPP4", "SGLT2", "TZD", "SU"), labels = c("MFN", "other", "other", "other", "other", "other"))) %>%
    rename("analysis_var" = "stopdrug_12m_6mFU") %>%
    table() %>%
    as.data.frame() %>%
    group_by(drugline, drugclass) %>%
    mutate(total = sum(Freq), Freq = Freq/total) %>%
    filter(analysis_var == 1) %>%
    ggplot(aes(y = drugline, x = Freq, fill = drugclass)) +
    geom_col() +
    scale_fill_manual(values = c("other" = "black", "SGLT2" = "#E69F00", "GLP1" = "#56B4E9", "SU" = "#CC79A7", "DPP4" = "#0072B2", "TZD" = "#D55E00", "MFN" = "grey"), breaks = c("MFN", "other"), labels = c("Metformin", "Other glucose-lowering therapies"), name = " ", guide = guide_legend(reverse = FALSE)) +
    scale_y_discrete(limits = rev) +
    ggtitle("Number of therapy classes ever prescribed") +
    scale_x_continuous("Proportion of discontinuation", labels = scales::percent, limits = c(0, 0.3)) +
    facet_wrap(~drugclass, nrow = 1, labeller = as_labeller(c("MFN" = "Metformin", "other" = "Other glucose-lowering therapies"))) +
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
  cprd_dataset %>%
    select(predrug_frailty_proxy, stopdrug_12m_6mFU, drugclass) %>%
    mutate(drugclass = factor(drugclass, levels = c("MFN", "GLP1", "DPP4", "SGLT2", "TZD", "SU"), labels = c("MFN", "other", "other", "other", "other", "other"))) %>%
    rename("analysis_var" = "stopdrug_12m_6mFU") %>%
    table() %>%
    as.data.frame() %>%
    group_by(predrug_frailty_proxy, drugclass) %>%
    mutate(total = sum(Freq), Freq = Freq/total) %>%
    filter(analysis_var == 1) %>%
    ggplot(aes(y = predrug_frailty_proxy, x = Freq, fill = drugclass)) +
    geom_col() +
    scale_fill_manual(values = c("other" = "black", "SGLT2" = "#E69F00", "GLP1" = "#56B4E9", "SU" = "#CC79A7", "DPP4" = "#0072B2", "TZD" = "#D55E00", "MFN" = "grey"), breaks = c("MFN", "other"), labels = c("Metformin", "Other glucose-lowering therapies"), name = " ", guide = guide_legend(reverse = FALSE)) +
    scale_y_discrete(limits = rev, breaks = c(0, 1), labels = c("No", "Yes")) +
    ggtitle("Frailty") +
    scale_x_continuous("Proportion of discontinuation", labels = scales::percent, limits = c(0, 0.3)) +
    facet_wrap(~drugclass, nrow = 1, labeller = as_labeller(c("MFN" = "Metformin", "other" = "Other glucose-lowering therapies"))) +
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
  cprd_dataset %>%
    select(ethnicity_5cat, stopdrug_12m_6mFU, drugclass) %>%
    mutate(ethnicity_5cat = factor(ethnicity_5cat, levels = c("4", "3", "2", "1", "0"), labels = c("Mixed & Other", "Mixed & Other", "Black", "South Asian", "White"))) %>%
    mutate(drugclass = factor(drugclass, levels = c("MFN", "GLP1", "DPP4", "SGLT2", "TZD", "SU"), labels = c("MFN", "other", "other", "other", "other", "other"))) %>%
    rename("analysis_var" = "stopdrug_12m_6mFU") %>%
    table() %>%
    as.data.frame() %>%
    group_by(ethnicity_5cat, drugclass) %>%
    mutate(total = sum(Freq), Freq = Freq/total) %>%
    filter(analysis_var == 1) %>%
    ggplot(aes(y = ethnicity_5cat, x = Freq, fill = drugclass)) +
    geom_col() +
    scale_fill_manual(values = c("other" = "black", "SGLT2" = "#E69F00", "GLP1" = "#56B4E9", "SU" = "#CC79A7", "DPP4" = "#0072B2", "TZD" = "#D55E00", "MFN" = "grey"), breaks = c("MFN", "other"), labels = c("Metformin", "Other glucose-lowering therapies"), name = " ", guide = guide_legend(reverse = FALSE)) +
    ggtitle("Ethnicity") +
    scale_x_continuous("Proportion of discontinuation", labels = scales::percent, limits = c(0, 0.3)) +
    facet_wrap(~drugclass, nrow = 1, labeller = as_labeller(c("MFN" = "Metformin", "other" = "Other glucose-lowering therapies"))) +
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
    height = c(1, 4, 2, 5, 4, 5, 4, 3, 2, 4),
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


pdf("results/figures/03.12m_discontinuation_by_variable_by_drug.pdf", width = 14, height = 13)

patchwork::wrap_plots(
  
  # Age - Metformin
  cprd_dataset %>%
    select(dstartdate_age, stopdrug_12m_6mFU, drugclass) %>%
    mutate(drugclass = factor(drugclass, levels = c("MFN", "GLP1", "DPP4", "SGLT2", "TZD", "SU"))) %>%
    mutate(group = cut(dstartdate_age, breaks = c(18, 50, 60, 70, 150))) %>%
    rename("analysis_var" = "stopdrug_12m_6mFU") %>%
    select(-dstartdate_age) %>%
    table() %>%
    as.data.frame() %>%
    group_by(group, drugclass) %>%
    mutate(total = sum(Freq), Freq = Freq/total) %>%
    filter(analysis_var == 1) %>%
    ggplot(aes(y = group, x = Freq, fill = drugclass)) +
    geom_col(position = "dodge") +
    scale_fill_manual(values = c("Pooled" = "black", "SGLT2" = "#E69F00", "GLP1" = "#56B4E9", "SU" = "#CC79A7", "DPP4" = "#0072B2", "TZD" = "#D55E00", "MFN" = "grey"), breaks = c("MFN", "GLP1", "DPP4", "SGLT2", "TZD", "SU"), labels = c("Metformin", "GLP-1RA", "DPP4i", "SGLT2i", "TZD", "SU"), name = " ", guide = guide_legend(reverse = TRUE)) +
    scale_y_discrete(limits = rev, labels = c("70+", "60-69", "50-59", "<50")) +
    ggtitle("Age group, years") +
    scale_x_continuous("Proportion of discontinuation", labels = scales::percent, limits = c(0, 0.5), breaks = seq(0, 0.5, 0.1)) +
    facet_wrap(~drugclass, nrow = 1, labeller = as_labeller(c("MFN" = "Metformin", "GLP1" = "GLP-1RA", "DPP4" = "DPP4i", "SGLT2" = "SGLT2i", "TZD" = "TZD", "SU" = "SU"))) +
    theme_bw() +
    guides(fill = guide_legend(nrow = 1, byrow = TRUE)) +
    theme(
      legend.position = "bottom",
      axis.text.x = element_blank(),
      axis.title.x = element_blank()
    ),
  
  # Sex 
  cprd_dataset %>%
    select(gender, stopdrug_12m_6mFU, drugclass) %>%
    mutate(drugclass = factor(drugclass, levels = c("MFN", "GLP1", "DPP4", "SGLT2", "TZD", "SU"))) %>%
    rename("analysis_var" = "stopdrug_12m_6mFU") %>%
    table() %>%
    as.data.frame() %>%
    group_by(gender, drugclass) %>%
    mutate(total = sum(Freq), Freq = Freq/total) %>%
    filter(analysis_var == 1) %>%
    ggplot(aes(y = gender, x = Freq, fill = drugclass)) +
    geom_col(position = "dodge") +
    scale_fill_manual(values = c("Pooled" = "black", "SGLT2" = "#E69F00", "GLP1" = "#56B4E9", "SU" = "#CC79A7", "DPP4" = "#0072B2", "TZD" = "#D55E00", "MFN" = "grey"), breaks = c("MFN", "GLP1", "DPP4", "SGLT2", "TZD", "SU"), labels = c("Metformin", "GLP-1RA", "DPP4i", "SGLT2i", "TZD", "SU"), name = " ", guide = guide_legend(reverse = TRUE)) +
    scale_y_discrete(limits = rev, labels = c("Female", "Male")) +
    ggtitle("Sex") +
    scale_x_continuous("Proportion of discontinuation", labels = scales::percent, limits = c(0, 0.5), breaks = seq(0, 0.5, 0.1)) +
    facet_wrap(~drugclass, nrow = 1, labeller = as_labeller(c("MFN" = "Metformin", "GLP1" = "GLP-1RA", "DPP4" = "DPP4i", "SGLT2" = "SGLT2i", "TZD" = "TZD", "SU" = "SU"))) +
    theme_bw() +
    guides(fill = guide_legend(nrow = 1, byrow = TRUE)) +
    theme(
      legend.position = "bottom",
      axis.text.x = element_blank(),
      axis.title.x = element_blank(),
      strip.text.x = element_blank()
    ),
  
  # Quintile
  cprd_dataset %>%
    select(imd2015_10, stopdrug_12m_6mFU, drugclass) %>%
    mutate(drugclass = factor(drugclass, levels = c("MFN", "GLP1", "DPP4", "SGLT2", "TZD", "SU"))) %>%
    rename("analysis_var" = "stopdrug_12m_6mFU") %>%
    table() %>%
    as.data.frame() %>%
    group_by(imd2015_10, drugclass) %>%
    mutate(total = sum(Freq), Freq = Freq/total) %>%
    filter(analysis_var == 1) %>%
    ggplot(aes(y = imd2015_10, x = Freq, fill = drugclass)) +
    geom_col(position = "dodge") +
    scale_fill_manual(values = c("Pooled" = "black", "SGLT2" = "#E69F00", "GLP1" = "#56B4E9", "SU" = "#CC79A7", "DPP4" = "#0072B2", "TZD" = "#D55E00", "MFN" = "grey"), breaks = c("MFN", "GLP1", "DPP4", "SGLT2", "TZD", "SU"), labels = c("Metformin", "GLP-1RA", "DPP4i", "SGLT2i", "TZD", "SU"), name = " ", guide = guide_legend(reverse = TRUE)) +
    scale_y_discrete(limits = rev, labels = c("5 (most deprived)", "4", "3", "2", "1 (least deprived)")) +
    ggtitle("Index of multiple deprivation quintile") +
    scale_x_continuous("Proportion of discontinuation", labels = scales::percent, limits = c(0, 0.5), breaks = seq(0, 0.5, 0.1)) +
    facet_wrap(~drugclass, nrow = 1, labeller = as_labeller(c("MFN" = "Metformin", "GLP1" = "GLP-1RA", "DPP4" = "DPP4i", "SGLT2" = "SGLT2i", "TZD" = "TZD", "SU" = "SU"))) +
    theme_bw() +
    guides(fill = guide_legend(nrow = 1, byrow = TRUE)) +
    theme(
      legend.position = "bottom",
      axis.text.x = element_blank(),
      axis.title.x = element_blank(),
      strip.text.x = element_blank()
    ),
  
  # BMI
  cprd_dataset %>%
    select(prebmi, stopdrug_12m_6mFU, drugclass) %>%
    mutate(drugclass = factor(drugclass, levels = c("MFN", "GLP1", "DPP4", "SGLT2", "TZD", "SU"))) %>%
    mutate(group = cut(prebmi, breaks = c(1, 25, 30, 35, 500))) %>%
    rename("analysis_var" = "stopdrug_12m_6mFU") %>%
    select(-prebmi) %>%
    table() %>%
    as.data.frame() %>%
    group_by(group, drugclass) %>%
    mutate(total = sum(Freq), Freq = Freq/total) %>%
    filter(analysis_var == 1) %>%
    ggplot(aes(y = group, x = Freq, fill = drugclass)) +
    geom_col(position = "dodge") +
    scale_fill_manual(values = c("Pooled" = "black", "SGLT2" = "#E69F00", "GLP1" = "#56B4E9", "SU" = "#CC79A7", "DPP4" = "#0072B2", "TZD" = "#D55E00", "MFN" = "grey"), breaks = c("MFN", "GLP1", "DPP4", "SGLT2", "TZD", "SU"), labels = c("Metformin", "GLP-1RA", "DPP4i", "SGLT2i", "TZD", "SU"), name = " ", guide = guide_legend(reverse = TRUE)) +
    scale_y_discrete(limits = rev, labels = c("35+", "30-34.9", "25-29.9", "<25")) +
    ggtitle("BMI, kg/m2") +
    scale_x_continuous("Proportion of discontinuation", labels = scales::percent, limits = c(0, 0.5), breaks = seq(0, 0.5, 0.1)) +
    facet_wrap(~drugclass, nrow = 1, labeller = as_labeller(c("MFN" = "Metformin", "GLP1" = "GLP-1RA", "DPP4" = "DPP4i", "SGLT2" = "SGLT2i", "TZD" = "TZD", "SU" = "SU"))) +
    theme_bw() +
    guides(fill = guide_legend(nrow = 1, byrow = TRUE)) +
    theme(
      legend.position = "bottom",
      axis.text.x = element_blank(),
      axis.title.x = element_blank(),
      strip.text.x = element_blank()
    ),
  
  # Diabetes duration
  cprd_dataset %>%
    select(dstartdate_dm_dur, stopdrug_12m_6mFU, drugclass) %>%
    mutate(drugclass = factor(drugclass, levels = c("MFN", "GLP1", "DPP4", "SGLT2", "TZD", "SU"))) %>%
    mutate(group = cut(dstartdate_dm_dur, breaks = c(-1, 1, 3, 6, 10, 150))) %>%
    rename("analysis_var" = "stopdrug_12m_6mFU") %>%
    select(-dstartdate_dm_dur) %>%
    table() %>%
    as.data.frame() %>%
    group_by(group, drugclass) %>%
    mutate(total = sum(Freq), Freq = Freq/total) %>%
    filter(analysis_var == 1) %>%
    ggplot(aes(y = group, x = Freq, fill = drugclass)) +
    geom_col(position = "dodge") +
    scale_fill_manual(values = c("Pooled" = "black", "SGLT2" = "#E69F00", "GLP1" = "#56B4E9", "SU" = "#CC79A7", "DPP4" = "#0072B2", "TZD" = "#D55E00", "MFN" = "grey"), breaks = c("MFN", "GLP1", "DPP4", "SGLT2", "TZD", "SU"), labels = c("Metformin", "GLP-1RA", "DPP4i", "SGLT2i", "TZD", "SU"), name = " ", guide = guide_legend(reverse = TRUE)) +
    scale_y_discrete(limits = rev, labels = c("10+", "6-9", "3-5", "1-2", "<1")) +
    ggtitle("Duration of diabetes, years") +
    scale_x_continuous("Proportion of discontinuation", labels = scales::percent, limits = c(0, 0.5), breaks = seq(0, 0.5, 0.1)) +
    facet_wrap(~drugclass, nrow = 1, labeller = as_labeller(c("MFN" = "Metformin", "GLP1" = "GLP-1RA", "DPP4" = "DPP4i", "SGLT2" = "SGLT2i", "TZD" = "TZD", "SU" = "SU"))) +
    theme_bw() +
    guides(fill = guide_legend(nrow = 1, byrow = TRUE)) +
    theme(
      legend.position = "bottom",
      axis.text.x = element_blank(),
      axis.title.x = element_blank(),
      strip.text.x = element_blank()
    ),
  
  # HbA1c
  cprd_dataset %>%
    select(prehba1c, stopdrug_12m_6mFU, drugclass) %>%
    mutate(drugclass = factor(drugclass, levels = c("MFN", "GLP1", "DPP4", "SGLT2", "TZD", "SU"))) %>%
    mutate(group = cut(prehba1c, breaks = c(52, 64, 75, 86, 250))) %>%
    rename("analysis_var" = "stopdrug_12m_6mFU") %>%
    select(-prehba1c) %>%
    table() %>%
    as.data.frame() %>%
    group_by(group, drugclass) %>%
    mutate(total = sum(Freq), Freq = Freq/total) %>%
    filter(analysis_var == 1) %>%
    ggplot(aes(y = group, x = Freq, fill = drugclass)) +
    geom_col(position = "dodge") +
    scale_fill_manual(values = c("Pooled" = "black", "SGLT2" = "#E69F00", "GLP1" = "#56B4E9", "SU" = "#CC79A7", "DPP4" = "#0072B2", "TZD" = "#D55E00", "MFN" = "grey"), breaks = c("MFN", "GLP1", "DPP4", "SGLT2", "TZD", "SU"), labels = c("Metformin", "GLP-1RA", "DPP4i", "SGLT2i", "TZD", "SU"), name = " ", guide = guide_legend(reverse = TRUE)) +
    scale_y_discrete(limits = rev, labels = c("86+", "75-86", "64-75", "53-64")) +
    ggtitle("HbA1c, mmol/mol") +
    scale_x_continuous("Proportion of discontinuation", labels = scales::percent, limits = c(0, 0.5), breaks = seq(0, 0.5, 0.1)) +
    facet_wrap(~drugclass, nrow = 1, labeller = as_labeller(c("MFN" = "Metformin", "GLP1" = "GLP-1RA", "DPP4" = "DPP4i", "SGLT2" = "SGLT2i", "TZD" = "TZD", "SU" = "SU"))) +
    theme_bw() +
    guides(fill = guide_legend(nrow = 1, byrow = TRUE)) +
    theme(
      legend.position = "bottom",
      axis.text.x = element_blank(),
      axis.title.x = element_blank(),
      strip.text.x = element_blank()
    ),
  
  # drugline
  cprd_dataset %>%
    select(drugline, stopdrug_12m_6mFU, drugclass) %>%
    mutate(drugclass = factor(drugclass, levels = c("MFN", "GLP1", "DPP4", "SGLT2", "TZD", "SU"))) %>%
    mutate(drugline = factor(drugline, levels = c("1", "2", "3", "4", "5+"), labels = c("1", "2", "3+", "3+", "3+"))) %>%
    rename("analysis_var" = "stopdrug_12m_6mFU") %>%
    table() %>%
    as.data.frame() %>%
    group_by(drugline, drugclass) %>%
    mutate(total = sum(Freq), Freq = Freq/total) %>%
    filter(analysis_var == 1) %>%
    ggplot(aes(y = drugline, x = Freq, fill = drugclass)) +
    geom_col(position = "dodge") +
    scale_fill_manual(values = c("Pooled" = "black", "SGLT2" = "#E69F00", "GLP1" = "#56B4E9", "SU" = "#CC79A7", "DPP4" = "#0072B2", "TZD" = "#D55E00", "MFN" = "grey"), breaks = c("MFN", "GLP1", "DPP4", "SGLT2", "TZD", "SU"), labels = c("Metformin", "GLP-1RA", "DPP4i", "SGLT2i", "TZD", "SU"), name = " ", guide = guide_legend(reverse = TRUE)) +
    scale_y_discrete(limits = rev) +
    ggtitle("Number of therapy classes ever prescribed") +
    scale_x_continuous("Proportion of discontinuation", labels = scales::percent, limits = c(0, 0.5), breaks = seq(0, 0.5, 0.1)) +
    facet_wrap(~drugclass, nrow = 1, labeller = as_labeller(c("MFN" = "Metformin", "GLP1" = "GLP-1RA", "DPP4" = "DPP4i", "SGLT2" = "SGLT2i", "TZD" = "TZD", "SU" = "SU"))) +
    theme_bw() +
    guides(fill = guide_legend(nrow = 1, byrow = TRUE)) +
    theme(
      legend.position = "bottom",
      axis.text.x = element_blank(),
      axis.title.x = element_blank(),
      strip.text.x = element_blank()
    ),
  
  # Frailty
  cprd_dataset %>%
    select(predrug_frailty_proxy, stopdrug_12m_6mFU, drugclass) %>%
    mutate(drugclass = factor(drugclass, levels = c("MFN", "GLP1", "DPP4", "SGLT2", "TZD", "SU"))) %>%
    rename("analysis_var" = "stopdrug_12m_6mFU") %>%
    table() %>%
    as.data.frame() %>%
    group_by(predrug_frailty_proxy, drugclass) %>%
    mutate(total = sum(Freq), Freq = Freq/total) %>%
    filter(analysis_var == 1) %>%
    ggplot(aes(y = predrug_frailty_proxy, x = Freq, fill = drugclass)) +
    geom_col(position = "dodge") +
    scale_fill_manual(values = c("Pooled" = "black", "SGLT2" = "#E69F00", "GLP1" = "#56B4E9", "SU" = "#CC79A7", "DPP4" = "#0072B2", "TZD" = "#D55E00", "MFN" = "grey"), breaks = c("MFN", "GLP1", "DPP4", "SGLT2", "TZD", "SU"), labels = c("Metformin", "GLP-1RA", "DPP4i", "SGLT2i", "TZD", "SU"), name = " ", guide = guide_legend(reverse = TRUE)) +
    scale_y_discrete(limits = rev, breaks = c(0, 1), labels = c("No", "Yes")) +
    ggtitle("Frailty") +
    scale_x_continuous("Proportion of discontinuation", labels = scales::percent, limits = c(0, 0.5), breaks = seq(0, 0.5, 0.1)) +
    facet_wrap(~drugclass, nrow = 1, labeller = as_labeller(c("MFN" = "Metformin", "GLP1" = "GLP-1RA", "DPP4" = "DPP4i", "SGLT2" = "SGLT2i", "TZD" = "TZD", "SU" = "SU"))) +
    theme_bw() +
    guides(fill = guide_legend(nrow = 1, byrow = TRUE)) +
    theme(
      legend.position = "bottom",
      axis.text.x = element_blank(),
      axis.title.x = element_blank(),
      strip.text.x = element_blank()
    ),
  
  # Ethnicity
  cprd_dataset %>%
    select(ethnicity_5cat, stopdrug_12m_6mFU, drugclass) %>%
    mutate(drugclass = factor(drugclass, levels = c("MFN", "GLP1", "DPP4", "SGLT2", "TZD", "SU"))) %>%
    mutate(ethnicity_5cat = factor(ethnicity_5cat, levels = c("4", "3", "2", "1", "0"), labels = c("Mixed & Other", "Mixed & Other", "Black", "South Asian", "White"))) %>%
    rename("analysis_var" = "stopdrug_12m_6mFU") %>%
    table() %>%
    as.data.frame() %>%
    group_by(ethnicity_5cat, drugclass) %>%
    mutate(total = sum(Freq), Freq = Freq/total) %>%
    filter(analysis_var == 1) %>%
    ggplot(aes(y = ethnicity_5cat, x = Freq, fill = drugclass)) +
    geom_col(position = "dodge") +
    scale_fill_manual(values = c("Pooled" = "black", "SGLT2" = "#E69F00", "GLP1" = "#56B4E9", "SU" = "#CC79A7", "DPP4" = "#0072B2", "TZD" = "#D55E00", "MFN" = "grey"), breaks = c("MFN", "GLP1", "DPP4", "SGLT2", "TZD", "SU"), labels = c("Metformin", "GLP-1RA", "DPP4i", "SGLT2i", "TZD", "SU"), name = " ", guide = guide_legend(reverse = TRUE)) +
    ggtitle("Ethnicity") +
    scale_x_continuous("Proportion of discontinuation", labels = scales::percent, limits = c(0, 0.5), breaks = seq(0, 0.5, 0.1)) +
    facet_wrap(~drugclass, nrow = 1, labeller = as_labeller(c("MFN" = "Metformin", "GLP1" = "GLP-1RA", "DPP4" = "DPP4i", "SGLT2" = "SGLT2i", "TZD" = "TZD", "SU" = "SU"))) +
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
    height = c(4, 2, 5, 4, 5, 4, 3, 2, 4),
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
############################### Read Data In ##################################
###############################################################################
###############################################################################

# load dataset
cprd_dataset <- set_up_data(
  raw_data = "20240308_t2d_1stinstance",
  diagnosis = FALSE,
  therapies = c("DPP4", "GLP1", "MFN", "SGLT2", "SU", "TZD"),
  dataset = "3m.disc.dataset",
  full_prescribing_history = TRUE
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



coefficients_3m <- rbind(
  dstartdate_dm_dur_sum,
  dstartdate_age_sum,
  prehba1c_sum,
  preegfr_sum,
  prebmi_sum,
  prealt_sum,
  prehdl_sum,
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

saveRDS(coefficients_3m, "results/Models/Predictions/univariate_coefficients_3m.rds")


pdf("results/figures/03.univariate_analysis_3m.pdf", width = 12, height = 8)

wrap_plots(
  
  coefficients_3m %>%
    filter(type == "Behavioural") %>%
    filter(variable %in% c("Blood pressure medication", "Statins", "MFN discontinuation", "Ethnicity (ref White)", "Smoking (ref Active Smoker)")) %>%
    ggplot() + 
    geom_vline(aes(xintercept = 1), colour = "black") +
    geom_errorbar(aes(x = coef, xmin = LCI, xmax = UCI, y = variable, colour = drug),position = position_dodge(width = 1)) +
    geom_pointrange(aes(x = coef, xmin = LCI, xmax = UCI, y = variable, colour = drug), position = position_dodge(width = 1)) +
    xlab("Odds Ratio") +
    scale_x_log10(breaks = c(0.25, 0.5, 1, 2, 4), limits = c(0.4, 4.1)) +
    facet_grid(variable~., scales = "free_y") +
    scale_colour_manual(values = c("Pooled" = "black", "SGLT2" = "#E69F00", "GLP1" = "#56B4E9", "SU" = "#CC79A7", "DPP4" = "#0072B2", "TZD" = "#D55E00", "MFN" = "grey"), 
                        breaks = rev(c("Pooled", "MFN", "GLP1", "DPP4", "SGLT2", "TZD", "SU")), labels = rev(c("Pooled", "Metformin", "GLP-1RA", "DPP4i", "SGLT2i", "TZD", "SU")),
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
    ggplot() + 
    geom_vline(aes(xintercept = 1), colour = "black") +
    geom_errorbar(aes(x = coef, xmin = LCI, xmax = UCI, y = variable, colour = drug),position = position_dodge(width = 1)) +
    geom_pointrange(aes(x = coef, xmin = LCI, xmax = UCI, y = variable, colour = drug), position = position_dodge(width = 1)) +
    xlab("Odds Ratio") +
    scale_x_log10(breaks = c(0.25, 0.5, 1, 2, 4), limits = c(0.4, 4.1)) +
    facet_grid(variable~., scales = "free_y") +
    scale_colour_manual(values = c("Pooled" = "black", "SGLT2" = "#E69F00", "GLP1" = "#56B4E9", "SU" = "#CC79A7", "DPP4" = "#0072B2", "TZD" = "#D55E00", "MFN" = "grey"), 
                        breaks = rev(c("Pooled", "MFN", "GLP1", "DPP4", "SGLT2", "TZD", "SU")), labels = rev(c("Pooled", "Metformin", "GLP-1RA", "DPP4i", "SGLT2i", "TZD", "SU")),
                        name = "Therapy", guide = guide_legend(reverse = TRUE)) +
    theme_bw() +
    theme(
      axis.title.y = element_blank(),
      legend.position = "bottom",
      strip.text = element_blank()
    ),
  
  coefficients_3m %>%
    filter(type == "Comorbidity event") %>%
    filter(variable %in% c("Frailty proxy")) %>%
    ggplot() + 
    geom_vline(aes(xintercept = 1), colour = "black") +
    geom_errorbar(aes(x = coef, xmin = LCI, xmax = UCI, y = variable, colour = drug),position = position_dodge(width = 1)) +
    geom_pointrange(aes(x = coef, xmin = LCI, xmax = UCI, y = variable, colour = drug), position = position_dodge(width = 1)) +
    xlab("Odds Ratio") +
    scale_x_log10(breaks = c(0.25, 0.5, 1, 2, 4), limits = c(0.4, 4.1)) +
    facet_grid(variable~., scales = "free_y") +
    scale_colour_manual(values = c("Pooled" = "black", "SGLT2" = "#E69F00", "GLP1" = "#56B4E9", "SU" = "#CC79A7", "DPP4" = "#0072B2", "TZD" = "#D55E00", "MFN" = "grey"), 
                        breaks = rev(c("Pooled", "MFN", "GLP1", "DPP4", "SGLT2", "TZD", "SU")), labels = rev(c("Pooled", "Metformin", "GLP-1RA", "DPP4i", "SGLT2i", "TZD", "SU")),
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
################# Univariate analysis - 6m discontinuation ####################
###############################################################################
###############################################################################


# load dataset
cprd_dataset <- cprd_dataset %>%
  drop_na(stopdrug_6m_6mFU)


drugs = c("MFN", "DPP4", "GLP1", "SGLT2", "SU", "TZD")


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

# prealt

prealt_sum <- univariate_analysis(
  cprd_dataset, drugs = drugs, outcome = "stopdrug_6m_6mFU", 
  variable = "prealt", variable_name = "ALT (per SD)", type = "Clinical features & biomarkers"
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
  variable = "stopdrug_3m_3mFU_MFN_hist", variable_name = "MFN discontinuation", type = "Behavioural"
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


# alcohol_cat

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



coefficients_6m <- rbind(
  dstartdate_dm_dur_sum,
  dstartdate_age_sum,
  prehba1c_sum,
  preegfr_sum,
  prebmi_sum,
  prealt_sum,
  prehdl_sum,
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


saveRDS(coefficients_6m, "results/Models/Predictions/univariate_coefficients_6m.rds")



pdf("results/figures/03.univariate_analysis_6m.pdf", width = 12, height = 8)

wrap_plots(
  
  coefficients_6m %>%
    filter(type == "Behavioural") %>%
    filter(variable %in% c("Blood pressure medication", "Statins", "MFN discontinuation", "Ethnicity (ref White)", "Smoking (ref Active Smoker)")) %>%
    ggplot() + 
    geom_vline(aes(xintercept = 1), colour = "black") +
    geom_errorbar(aes(x = coef, xmin = LCI, xmax = UCI, y = variable, colour = drug),position = position_dodge(width = 1)) +
    geom_pointrange(aes(x = coef, xmin = LCI, xmax = UCI, y = variable, colour = drug), position = position_dodge(width = 1)) +
    xlab("Odds Ratio") +
    scale_x_log10(breaks = c(0.25, 0.5, 1, 2, 4), limits = c(0.4, 4.1)) +
    facet_grid(variable~., scales = "free_y") +
    scale_colour_manual(values = c("Pooled" = "black", "SGLT2" = "#E69F00", "GLP1" = "#56B4E9", "SU" = "#CC79A7", "DPP4" = "#0072B2", "TZD" = "#D55E00", "MFN" = "grey"), 
                        breaks = rev(c("Pooled", "MFN", "GLP1", "DPP4", "SGLT2", "TZD", "SU")), labels = rev(c("Pooled", "Metformin", "GLP-1RA", "DPP4i", "SGLT2i", "TZD", "SU")),
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
    ggplot() + 
    geom_vline(aes(xintercept = 1), colour = "black") +
    geom_errorbar(aes(x = coef, xmin = LCI, xmax = UCI, y = variable, colour = drug),position = position_dodge(width = 1)) +
    geom_pointrange(aes(x = coef, xmin = LCI, xmax = UCI, y = variable, colour = drug), position = position_dodge(width = 1)) +
    xlab("Odds Ratio") +
    scale_x_log10(breaks = c(0.25, 0.5, 1, 2, 4), limits = c(0.4, 4.1)) +
    facet_grid(variable~., scales = "free_y") +
    scale_colour_manual(values = c("Pooled" = "black", "SGLT2" = "#E69F00", "GLP1" = "#56B4E9", "SU" = "#CC79A7", "DPP4" = "#0072B2", "TZD" = "#D55E00", "MFN" = "grey"), 
                        breaks = rev(c("Pooled", "MFN", "GLP1", "DPP4", "SGLT2", "TZD", "SU")), labels = rev(c("Pooled", "Metformin", "GLP-1RA", "DPP4i", "SGLT2i", "TZD", "SU")),
                        name = "Therapy", guide = guide_legend(reverse = TRUE)) +
    theme_bw() +
    theme(
      axis.title.y = element_blank(),
      legend.position = "bottom",
      strip.text = element_blank()
    ),
  
  coefficients_6m %>%
    filter(type == "Comorbidity event") %>%
    filter(variable %in% c("Frailty proxy")) %>%
    ggplot() + 
    geom_vline(aes(xintercept = 1), colour = "black") +
    geom_errorbar(aes(x = coef, xmin = LCI, xmax = UCI, y = variable, colour = drug),position = position_dodge(width = 1)) +
    geom_pointrange(aes(x = coef, xmin = LCI, xmax = UCI, y = variable, colour = drug), position = position_dodge(width = 1)) +
    xlab("Odds Ratio") +
    scale_x_log10(breaks = c(0.25, 0.5, 1, 2, 4), limits = c(0.4, 4.1)) +
    facet_grid(variable~., scales = "free_y") +
    scale_colour_manual(values = c("Pooled" = "black", "SGLT2" = "#E69F00", "GLP1" = "#56B4E9", "SU" = "#CC79A7", "DPP4" = "#0072B2", "TZD" = "#D55E00", "MFN" = "grey"), 
                        breaks = rev(c("Pooled", "MFN", "GLP1", "DPP4", "SGLT2", "TZD", "SU")), labels = rev(c("Pooled", "Metformin", "GLP-1RA", "DPP4i", "SGLT2i", "TZD", "SU")),
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


drugs = c("MFN", "DPP4", "GLP1", "SGLT2", "SU", "TZD")


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

# prealt

prealt_sum <- univariate_analysis(
  cprd_dataset, drugs = drugs, outcome = "stopdrug_12m_6mFU", 
  variable = "prealt", variable_name = "ALT (per SD)", type = "Clinical features & biomarkers"
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
  variable = "stopdrug_3m_3mFU_MFN_hist", variable_name = "MFN discontinuation", type = "Behavioural"
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


# alcohol_cat

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



coefficients_12m <- rbind(
  dstartdate_dm_dur_sum,
  dstartdate_age_sum,
  prehba1c_sum,
  preegfr_sum,
  prebmi_sum,
  prealt_sum,
  prehdl_sum,
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


saveRDS(coefficients_12m, "results/Models/Predictions/univariate_coefficients_12m.rds")



pdf("results/figures/03.univariate_analysis_12m.pdf", width = 12, height = 8)



wrap_plots(
  
  coefficients_12m %>%
    filter(type == "Behavioural") %>%
    filter(variable %in% c("Blood pressure medication", "Statins", "MFN discontinuation", "Ethnicity (ref White)", "Smoking (ref Active Smoker)")) %>%
    ggplot() + 
    geom_vline(aes(xintercept = 1), colour = "black") +
    geom_errorbar(aes(x = coef, xmin = LCI, xmax = UCI, y = variable, colour = drug),position = position_dodge(width = 1)) +
    geom_pointrange(aes(x = coef, xmin = LCI, xmax = UCI, y = variable, colour = drug), position = position_dodge(width = 1)) +
    xlab("Odds Ratio") +
    scale_x_log10(breaks = c(0.25, 0.5, 1, 2, 4), limits = c(0.4, 4.1)) +
    facet_grid(variable~., scales = "free_y") +
    scale_colour_manual(values = c("Pooled" = "black", "SGLT2" = "#E69F00", "GLP1" = "#56B4E9", "SU" = "#CC79A7", "DPP4" = "#0072B2", "TZD" = "#D55E00", "MFN" = "grey"), 
                        breaks = rev(c("Pooled", "MFN", "GLP1", "DPP4", "SGLT2", "TZD", "SU")), labels = rev(c("Pooled", "Metformin", "GLP-1RA", "DPP4i", "SGLT2i", "TZD", "SU")),
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
    ggplot() + 
    geom_vline(aes(xintercept = 1), colour = "black") +
    geom_errorbar(aes(x = coef, xmin = LCI, xmax = UCI, y = variable, colour = drug),position = position_dodge(width = 1)) +
    geom_pointrange(aes(x = coef, xmin = LCI, xmax = UCI, y = variable, colour = drug), position = position_dodge(width = 1)) +
    xlab("Odds Ratio") +
    scale_x_log10(breaks = c(0.25, 0.5, 1, 2, 4), limits = c(0.4, 4.1)) +
    facet_grid(variable~., scales = "free_y") +
    scale_colour_manual(values = c("Pooled" = "black", "SGLT2" = "#E69F00", "GLP1" = "#56B4E9", "SU" = "#CC79A7", "DPP4" = "#0072B2", "TZD" = "#D55E00", "MFN" = "grey"), 
                        breaks = rev(c("Pooled", "MFN", "GLP1", "DPP4", "SGLT2", "TZD", "SU")), labels = rev(c("Pooled", "Metformin", "GLP-1RA", "DPP4i", "SGLT2i", "TZD", "SU")),
                        name = "Therapy", guide = guide_legend(reverse = TRUE)) +
    theme_bw() +
    theme(
      axis.title.y = element_blank(),
      legend.position = "bottom",
      strip.text = element_blank()
    ),
  
  coefficients_12m %>%
    filter(type == "Comorbidity event") %>%
    filter(variable %in% c("Frailty proxy")) %>%
    ggplot() + 
    geom_vline(aes(xintercept = 1), colour = "black") +
    geom_errorbar(aes(x = coef, xmin = LCI, xmax = UCI, y = variable, colour = drug),position = position_dodge(width = 1)) +
    geom_pointrange(aes(x = coef, xmin = LCI, xmax = UCI, y = variable, colour = drug), position = position_dodge(width = 1)) +
    xlab("Odds Ratio") +
    scale_x_log10(breaks = c(0.25, 0.5, 1, 2, 4), limits = c(0.4, 4.1)) +
    facet_grid(variable~., scales = "free_y") +
    scale_colour_manual(values = c("Pooled" = "black", "SGLT2" = "#E69F00", "GLP1" = "#56B4E9", "SU" = "#CC79A7", "DPP4" = "#0072B2", "TZD" = "#D55E00", "MFN" = "grey"), 
                        breaks = rev(c("Pooled", "MFN", "GLP1", "DPP4", "SGLT2", "TZD", "SU")), labels = rev(c("Pooled", "Metformin", "GLP-1RA", "DPP4i", "SGLT2i", "TZD", "SU")),
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
