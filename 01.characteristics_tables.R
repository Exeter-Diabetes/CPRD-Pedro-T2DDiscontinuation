######################################################################
##
##  In this file, we create characteristics tables for all therapies
##
######################################################################

# load working directory
setwd("Samples/T2D_Discontinuation")

# load functions
source("code/00.set_up_data_and_functions.R")

# load libraries
library(tidyverse)
library(ggplot2)
library(tableone)
library(naniar)
library(patchwork)

# load dataset
cprd_dataset <- set_up_data(
  raw_data = "20240814_t2d_1stinstance",
  diagnosis = TRUE,
  therapies = c("MFN", "GLP1", "DPP4", "SGLT2", "SU", "TZD"),
  dataset = "full.dataset",
  follow_up = "6-months",
  full_prescribing_history = TRUE
)

## Plot of initiations by year
plot_initiations_by_year <- patchwork::wrap_plots(
  cprd_dataset %>%
    select(drugclass, dstartdate) %>%
    mutate(
      dstartdate = format(as.Date(dstartdate, format="%d/%m/%Y"),"%Y"),
      dstartdate = as.numeric(dstartdate)
    ) %>%
    ggplot(aes(x = dstartdate, fill = drugclass)) +
    geom_bar(position="stack") +
    scale_y_continuous("Therapy Initiations") +
    scale_x_continuous("Calendar Year", breaks = 2014:2020) +
    scale_fill_manual(breaks = c("MFN", "GLP1", "DPP4", "SGLT2", "SU", "TZD"), values = c("#009E73", "#56B4E9", "#0072B2", "#E69F00", "#CC79A7", "#D55E00")) +
    theme_bw() +
    guides(
      fill = guide_legend("Therapies:", nrow = 1)
    ) +
    theme(
      legend.position = "bottom"
    ),
  cprd_dataset %>%
    select(drugclass, dstartdate) %>%
    mutate(
      dstartdate = format(as.Date(dstartdate, format="%d/%m/%Y"),"%Y"),
      dstartdate = as.numeric(dstartdate)
    ) %>%
    ggplot(aes(x = dstartdate, fill = drugclass)) +
    geom_bar(position="fill") +
    scale_y_continuous("Therapy Initiations (%)", labels = scales::percent) +
    scale_x_continuous("Calendar Year", breaks = 2014:2020) +
    scale_fill_manual(breaks = c("MFN", "GLP1", "DPP4", "SGLT2", "SU", "TZD"), values = c("#009E73", "#56B4E9", "#0072B2", "#E69F00", "#CC79A7", "#D55E00")) +
    theme_bw() +
    guides(
      fill = guide_legend("Therapies:", nrow = 1)
    ) +
    theme(
      legend.position = "bottom"
    )
) +
  plot_layout(guides = 'collect') &
  theme(
    legend.position = "bottom",
    axis.text = element_text(size = 15),
    axis.title = element_text(size = 20),
    legend.title = element_text(size = 20),
    legend.text = element_text(size = 20)
  )

pdf("results/figures/01.plot_initiations_by_year.pdf", width = 11, height = 6)
plot_initiations_by_year
dev.off()


## Plot of discontinuations by year
plot_discontinuation_by_year <- patchwork::wrap_plots(
  cprd_dataset %>%
    select(stopdrug_3m_6mFU, dstartdate) %>%
    mutate(
      dstartdate = format(as.Date(dstartdate, format="%d/%m/%Y"),"%Y"),
      dstartdate = as.numeric(dstartdate),
      stopdrug_3m_6mFU = ifelse(is.na(stopdrug_3m_6mFU), "Missing", stopdrug_3m_6mFU),
      stopdrug_3m_6mFU = factor(stopdrug_3m_6mFU, levels = c("1", "2", "Missing"), labels = c("No", "Yes", "Missing"))
    ) %>%
    ggplot(aes(x = dstartdate, fill = stopdrug_3m_6mFU)) +
    geom_bar(position="stack") +
    scale_y_continuous("Therapy Initiations") +
    scale_x_continuous("Calendar Year", breaks = 2014:2020) +
    scale_fill_manual(breaks = c("No", "Yes", "Missing"), values = c("#F8766D", "#00BFC4", "grey")) +
    theme_bw() +
    guides(
      fill = guide_legend("3-month discontinuation:")
    ) +
    theme(
      legend.position = "bottom"
    ),
  cprd_dataset %>%
    select(stopdrug_3m_6mFU, dstartdate, drugclass) %>%
    mutate(
      dstartdate = format(as.Date(dstartdate, format="%d/%m/%Y"),"%Y"),
      dstartdate = as.numeric(dstartdate),
      stopdrug_3m_6mFU = ifelse(is.na(stopdrug_3m_6mFU), "Missing", stopdrug_3m_6mFU),
      stopdrug_3m_6mFU = factor(stopdrug_3m_6mFU, levels = c("1", "2", "Missing"), labels = c("No", "Yes", "Missing"))
    ) %>%
    ggplot(aes(x = dstartdate, fill = stopdrug_3m_6mFU)) +
    geom_bar(position="fill") +
    scale_y_continuous("Therapy Discontinuations (%)", labels = scales::percent) +
    scale_x_continuous("Calendar Year", breaks = 2014:2020) +
    scale_fill_manual(breaks = c("No", "Yes", "Missing"), values = c("#F8766D", "#00BFC4", "grey")) +
    theme_bw() +
    guides(
      fill = guide_legend("3-month discontinuation:")
    ) +
    theme(
      legend.position = "bottom"
    )
) +
  plot_layout(guides = 'collect') &
  theme(
    legend.position = "bottom",
    axis.text = element_text(size = 15),
    axis.title = element_text(size = 20),
    legend.title = element_text(size = 20),
    legend.text = element_text(size = 20)
  )

pdf("results/figures/01.plot_discontinuation_by_year.pdf", width = 11, height = 6)
plot_discontinuation_by_year
dev.off()


## Plot of discontinuations by drugclass
plot_discontinuation_by_drugclass <- patchwork::wrap_plots(
  cprd_dataset %>%
    select(drugclass, stopdrug_3m_6mFU) %>%
    mutate(
      stopdrug_3m_6mFU = ifelse(is.na(stopdrug_3m_6mFU), "Missing", stopdrug_3m_6mFU),
      stopdrug_3m_6mFU = factor(stopdrug_3m_6mFU, levels = c("1", "2", "Missing"), labels = c("No", "Yes", "Missing"))
    ) %>%
    ggplot(aes(x = drugclass, fill = stopdrug_3m_6mFU)) +
    geom_bar(position="stack") +
    scale_y_continuous("Therapy Initiations") +
    scale_x_discrete("Therapy Initiated") +
    scale_fill_manual(breaks = c("No", "Yes", "Missing"), values = c("#F8766D", "#00BFC4", "grey")) +
    theme_bw() +
    guides(
      fill = guide_legend("3-month discontinuation:")
    ) +
    theme(
      legend.position = "bottom"
    ),
  cprd_dataset %>%
    select(drugclass, stopdrug_3m_6mFU) %>%
    mutate(
      stopdrug_3m_6mFU = ifelse(is.na(stopdrug_3m_6mFU), "Missing", stopdrug_3m_6mFU),
      stopdrug_3m_6mFU = factor(stopdrug_3m_6mFU, levels = c("1", "2", "Missing"), labels = c("No", "Yes", "Missing"))
    ) %>%
    ggplot(aes(x = drugclass, fill = stopdrug_3m_6mFU)) +
    geom_bar(position="fill") +
    scale_y_continuous("Therapy Discontinuations (%)", labels = scales::percent) +
    scale_x_discrete("Therapy Initiated") +
    scale_fill_manual(breaks = c("No", "Yes", "Missing"), values = c("#F8766D", "#00BFC4", "grey")) +
    theme_bw() +
    guides(
      fill = guide_legend("3-month discontinuation:")
    ) +
    theme(
      legend.position = "bottom"
    )
) +
  plot_layout(guides = 'collect') &
  theme(
    legend.position = "bottom",
    axis.text = element_text(size = 15),
    axis.title = element_text(size = 20),
    legend.title = element_text(size = 20),
    legend.text = element_text(size = 20)
  )

pdf("results/figures/01.plot_discontinuation_by_drugclass.pdf", width = 11, height = 6)
plot_discontinuation_by_drugclass
dev.off()





#########################################################################################################
#########################################################################################################
#
# Table of Characteristics
#
#########################################################################################################
#########################################################################################################


# load dataset
cprd_tables <- set_up_data(
  raw_data = "20240814_t2d_1stinstance",
  diagnosis = TRUE,
  therapies = c("GLP1", "DPP4", "SGLT2", "SU", "TZD"),
  dataset = "3m.disc.dataset",
  follow_up = "6-months",
  full_prescribing_history = TRUE
) %>%
  # drop_na(-stopdrug_12m_6mFU)
  # drop_na(-prehdl, -stopdrug_6m_6mFU, -stopdrug_12m_6mFU) %>%
  mutate(drugline = factor(drugline, levels = c("1", "2", "3", "4", "5+"), labels = c("1", "2", "3", "4+", "4+")))


# unique patients
unique(cprd_tables$patid) %>% length()

output_path <- "results/tables"


interim_table <- cprd_tables %>%
  mutate(
    dstartdate_dm_dur_NA = ifelse(is.na(dstartdate_dm_dur), 1, 0),
    dstartdate_age_NA = ifelse(is.na(dstartdate_age), 1, 0),
    prehba1c_NA = ifelse(is.na(prehba1c), 1, 0),
    preegfr_NA = ifelse(is.na(preegfr), 1, 0),
    prebmi_NA = ifelse(is.na(prebmi), 1, 0)
  )

vars <- c(
  # Outcome
  # "stopdrug_3m_6mFU",
  # "stopdrug_6m_6mFU",
  # "stopdrug_12m_6mFU",
  # Extra info
  "dstartdate_age", "dstartdate_age_NA", "dstartdate_dm_dur", "dstartdate_dm_dur_NA", "gender", "smoking_cat", 
  "ethnicity_5cat", "imd2015_10",
  # Biomarkers
  "prebmi", "prebmi_NA", "prehba1c", "prehba1c_NA", "preegfr", "preegfr_NA", 
  "predrug_statins", "predrug_bloodmed", "predrug_frailty_proxy", "stopdrug_3m_3mFU_MFN_hist",
  "numdrugs", "drugline"
  # Comorbidities
)



vars_cat <- c(
  # Outcome
  # "stopdrug_3m_6mFU",
  # "stopdrug_6m_6mFU",
  # "stopdrug_12m_6mFU",
  # Extra info
  "drugline", "numdrugs", "smoking_cat", "imd2015_10",
  "predrug_statins", "stopdrug_3m_3mFU_MFN_hist", "ethnicity_5cat", "gender", "predrug_bloodmed",
  "dstartdate_age_NA", "dstartdate_dm_dur_NA", "prebmi_NA", "prehba1c_NA", "preegfr_NA", 
  # Comorbidities
  ## Frailty proxy
  "predrug_frailty_proxy"
)

### Table by drugclass

table_characteristics_drugs <- CreateTableOne(
  vars = vars,
  factorVars = vars_cat,
  includeNA = TRUE,
  strata = c("drugclass"),
  data = interim_table,
  test = TRUE
)

table_characteristics_drugs_print <- print(table_characteristics_drugs, exact = "stage", quote = FALSE, noSpaces = TRUE, printToggle = FALSE, contDigits = 1)

write.csv(table_characteristics_drugs_print, file = paste0(output_path, "/table_characteristics_drugs.csv"))

### Discontinuation

vars <- c(
  # Outcome
  "stopdrug_3m_6mFU",
  "stopdrug_6m_6mFU",
  "stopdrug_12m_6mFU"
)

vars_cat <- c(
  # Outcome
  "stopdrug_3m_6mFU",
  "stopdrug_6m_6mFU",
  "stopdrug_12m_6mFU"
)

table_characteristics_outcome <- CreateTableOne(
  vars = vars,
  factorVars = vars_cat,
  includeNA = TRUE,
  strata = c("drugclass"),
  data = interim_table,
  test = TRUE
)

table_characteristics_outcome_print <- print(table_characteristics_outcome, exact = "stage", quote = FALSE, noSpaces = TRUE, printToggle = FALSE, contDigits = 1)

write.csv(table_characteristics_outcome_print, file = paste0(output_path, "/table_characteristics_outcome.csv"))

### Custom table


interim_table <- cprd_tables %>%
  mutate(
    dstartdate_age = ifelse(dstartdate_age < 50, "<50", ifelse(dstartdate_age >= 50 & dstartdate_age < 60, "50-59", ifelse(dstartdate_age >=60 & dstartdate_age < 70, "60-69", "70+"))),
    dstartdate_dm_dur = ifelse(dstartdate_dm_dur < 3, "<2", ifelse(dstartdate_dm_dur >= 3 & dstartdate_dm_dur < 6, "3-5", ifelse(dstartdate_dm_dur >= 6 & dstartdate_dm_dur < 10, "6-9", "10+"))),
    preegfr = ifelse(preegfr < 75, "<75", ifelse(preegfr >= 75 & preegfr < 90, "75-90", "90+")),
    prebmi = ifelse(prebmi < 30, "<30", ifelse(prebmi >= 30 & prebmi < 35, "30-35", "35+")),
    prehba1c = ifelse(prehba1c < 64, "53-64", ifelse(prehba1c >= 64 & prehba1c < 75, "64-75", ifelse(prehba1c >= 75 & prehba1c < 86, "75-86", "86+")))
  )


vars <- c(
  "dstartdate_age", "dstartdate_dm_dur",
  "prebmi", "prehba1c", "preegfr"
)

table_characteristics_drugs_v2 <- CreateTableOne(
  vars = vars,
  factorVars = vars,
  includeNA = TRUE,
  strata = c("drugclass"),
  data = interim_table,
  test = TRUE
)

table_characteristics_drugs_v2_print <- print(table_characteristics_drugs_v2, exact = "stage", quote = FALSE, noSpaces = TRUE, printToggle = FALSE, contDigits = 1)

write.csv(table_characteristics_drugs_v2_print, file = paste0(output_path, "/table_characteristics_drugs_v2.csv"))




### All drugs

table_characteristics_all <- CreateTableOne(
  vars = vars,
  factorVars = vars_cat,
  includeNA = TRUE,
  strata = c("stopdrug_3m_6mFU"),
  data = cprd_tables,
  test = TRUE
)

table_characteristics_all_print <- print(table_characteristics_all, exact = "stage", quote = FALSE, noSpaces = TRUE, printToggle = FALSE)

write.csv(table_characteristics_all_print, file = paste0(output_path, "/table_characteristics_all.csv"))


### Per discontinuation per drug

table_characteristics_disc_drug <- CreateTableOne(
  vars = vars,
  factorVars = vars_cat,
  includeNA = TRUE,
  strata = c("stopdrug_3m_6mFU", "drugclass"),
  data = cprd_tables,
  test = TRUE
)

table_characteristics_disc_drug_print <- print(table_characteristics_disc_drug, exact = "stage", quote = FALSE, noSpaces = TRUE, printToggle = FALSE)

write.csv(table_characteristics_disc_drug_print, file = paste0(output_path, "/table_characteristics_disc_drug.csv"))







#########################################################################################################
#########################################################################################################
#
# Summary plots
#
#########################################################################################################
#########################################################################################################

plot_discontinuation_drugclass <- cprd_dataset %>%
  group_by(drugclass) %>%
  summarise(
    `By 3-months` = sum(as.numeric(stopdrug_3m_6mFU)-1, na.rm = TRUE)/n(),
    `By 6-months` = sum(as.numeric(stopdrug_6m_6mFU)-1, na.rm = TRUE)/n()
  ) %>%
  ungroup() %>%
  gather("Discontinuation", "value", -drugclass) %>%
  mutate(
    category = paste(drugclass, Discontinuation, sep = "_"),
    category = factor(category, levels = rev(c("MFN_By 3-months", "MFN_By 6-months", 
                                               "DPP4_By 3-months", "DPP4_By 6-months", 
                                               "SGLT2_By 3-months", "SGLT2_By 6-months", 
                                               "GLP1_By 3-months", "GLP1_By 6-months", 
                                               "SU_By 3-months", "SU_By 6-months", 
                                               "TZD_By 3-months", "TZD_By 6-months"))
    )
  ) %>%
  ggplot() +
  geom_col(aes(y = category, x = value, fill = Discontinuation), position = "dodge") +
  theme_bw() +
  scale_y_discrete("Treatments", labels = rev(c("Metformin", "", 
                                                "DPP4i", "", 
                                                "SGLT2i", "", 
                                                "GLP-1RA", "", 
                                                "SU", "", 
                                                "TZD", ""))) +
  scale_x_continuous("Discontinuation", labels = scales::percent) +
  theme(
    legend.position = "bottom",
    axis.ticks.y = element_blank()
  )


plot_discontinuation_age <- patchwork::wrap_plots(
  cprd_dataset %>%
    filter(drugclass == "MFN" & !is.na(stopdrug_3m_6mFU)) %>%
    mutate(Discontinuation = factor(stopdrug_3m_6mFU, levels = c(0, 1), labels = c("No", "Yes"))) %>%
    ggplot() +
    geom_density(aes(x = dstartdate_age, fill = Discontinuation), alpha = 0.6) +
    lims(x = c(0, 110)) +
    labs(y = "MFN") +
    theme_minimal() +
    theme(axis.text = element_blank(),
          axis.ticks = element_blank(),
          axis.title.x = element_blank()),
  cprd_dataset %>%
    filter(drugclass == "MFN" & !is.na(stopdrug_6m_6mFU)) %>%
    mutate(Discontinuation = factor(stopdrug_6m_6mFU, levels = c(0, 1), labels = c("No", "Yes"))) %>%
    ggplot() +
    geom_density(aes(x = dstartdate_age, fill = Discontinuation), alpha = 0.6) +
    lims(x = c(0, 110)) +
    labs(y = "") +
    theme_minimal() +
    theme(axis.text = element_blank(),
          axis.ticks = element_blank(),
          axis.title.x = element_blank()),
  cprd_dataset %>%
    filter(drugclass == "DPP4" & !is.na(stopdrug_3m_6mFU)) %>%
    mutate(Discontinuation = factor(stopdrug_3m_6mFU, levels = c(0, 1), labels = c("No", "Yes"))) %>%
    ggplot() +
    geom_density(aes(x = dstartdate_age, fill = Discontinuation), alpha = 0.6) +
    lims(x = c(0, 110)) +
    labs(y = "DPP4") +
    theme_minimal() +
    theme(axis.text = element_blank(),
          axis.ticks = element_blank(),
          axis.title.x = element_blank()),
  cprd_dataset %>%
    filter(drugclass == "DPP4" & !is.na(stopdrug_6m_6mFU)) %>%
    mutate(Discontinuation = factor(stopdrug_6m_6mFU, levels = c(0, 1), labels = c("No", "Yes"))) %>%
    ggplot() +
    geom_density(aes(x = dstartdate_age, fill = Discontinuation), alpha = 0.6) +
    lims(x = c(0, 110)) +
    labs(y = "") +
    theme_minimal() +
    theme(axis.text = element_blank(),
          axis.ticks = element_blank(),
          axis.title.x = element_blank()),
  cprd_dataset %>%
    filter(drugclass == "SGLT2" & !is.na(stopdrug_3m_6mFU)) %>%
    mutate(Discontinuation = factor(stopdrug_3m_6mFU, levels = c(0, 1), labels = c("No", "Yes"))) %>%
    ggplot() +
    geom_density(aes(x = dstartdate_age, fill = Discontinuation), alpha = 0.6) +
    lims(x = c(0, 110)) +
    labs(y = "SGLT2") +
    theme_minimal() +
    theme(axis.text = element_blank(),
          axis.ticks = element_blank(),
          axis.title.x = element_blank()),
  cprd_dataset %>%
    filter(drugclass == "SGLT2" & !is.na(stopdrug_6m_6mFU)) %>%
    mutate(Discontinuation = factor(stopdrug_6m_6mFU, levels = c(0, 1), labels = c("No", "Yes"))) %>%
    ggplot() +
    geom_density(aes(x = dstartdate_age, fill = Discontinuation), alpha = 0.6) +
    lims(x = c(0, 110)) +
    labs(y = "") +
    theme_minimal() +
    theme(axis.text = element_blank(),
          axis.ticks = element_blank(),
          axis.title.x = element_blank()),
  cprd_dataset %>%
    filter(drugclass == "GLP1" & !is.na(stopdrug_3m_6mFU)) %>%
    mutate(Discontinuation = factor(stopdrug_3m_6mFU, levels = c(0, 1), labels = c("No", "Yes"))) %>%
    ggplot() +
    geom_density(aes(x = dstartdate_age, fill = Discontinuation), alpha = 0.6) +
    lims(x = c(0, 110)) +
    labs(y = "GLP1") +
    theme_minimal() +
    theme(axis.text = element_blank(),
          axis.ticks = element_blank(),
          axis.title.x = element_blank()),
  cprd_dataset %>%
    filter(drugclass == "GLP1" & !is.na(stopdrug_6m_6mFU)) %>%
    mutate(Discontinuation = factor(stopdrug_6m_6mFU, levels = c(0, 1), labels = c("No", "Yes"))) %>%
    ggplot() +
    geom_density(aes(x = dstartdate_age, fill = Discontinuation), alpha = 0.6) +
    lims(x = c(0, 110)) +
    labs(y = "") +
    theme_minimal() +
    theme(axis.text = element_blank(),
          axis.ticks = element_blank(),
          axis.title.x = element_blank()),
  cprd_dataset %>%
    filter(drugclass == "SU" & !is.na(stopdrug_3m_6mFU)) %>%
    mutate(Discontinuation = factor(stopdrug_3m_6mFU, levels = c(0, 1), labels = c("No", "Yes"))) %>%
    ggplot() +
    geom_density(aes(x = dstartdate_age, fill = Discontinuation), alpha = 0.6) +
    lims(x = c(0, 110)) +
    labs(y = "SU") +
    theme_minimal() +
    theme(axis.text = element_blank(),
          axis.ticks = element_blank(),
          axis.title.x = element_blank()),
  cprd_dataset %>%
    filter(drugclass == "SU" & !is.na(stopdrug_6m_6mFU)) %>%
    mutate(Discontinuation = factor(stopdrug_6m_6mFU, levels = c(0, 1), labels = c("No", "Yes"))) %>%
    ggplot() +
    geom_density(aes(x = dstartdate_age, fill = Discontinuation), alpha = 0.6) +
    lims(x = c(0, 110)) +
    labs(y = "") +
    theme_minimal() +
    theme(axis.text = element_blank(),
          axis.ticks = element_blank(),
          axis.title.x = element_blank()),
  cprd_dataset %>%
    filter(drugclass == "TZD" & !is.na(stopdrug_3m_6mFU)) %>%
    mutate(Discontinuation = factor(stopdrug_3m_6mFU, levels = c(0, 1), labels = c("No", "Yes"))) %>%
    ggplot() +
    geom_density(aes(x = dstartdate_age, fill = Discontinuation), alpha = 0.6) +
    lims(x = c(0, 110)) +
    labs(y = "TZD") +
    theme_minimal() +
    theme(axis.text = element_blank(),
          axis.ticks = element_blank(),
          axis.title.x = element_blank()),
  cprd_dataset %>%
    filter(drugclass == "TZD" & !is.na(stopdrug_6m_6mFU)) %>%
    mutate(Discontinuation = factor(stopdrug_6m_6mFU, levels = c(0, 1), labels = c("No", "Yes"))) %>%
    ggplot() +
    geom_density(aes(x = dstartdate_age, fill = Discontinuation), alpha = 0.6) +
    lims(x = c(0, 110)) +
    labs(y = "", x = "Baseline Age") +
    theme_minimal() +
    theme(axis.text.y = element_blank(),
          axis.ticks.y = element_blank()),
  ncol = 1
) +
  plot_layout(guides = 'collect') & 
  theme(
    legend.position = "bottom",
    axis.title.y = element_text(angle = 0, vjust = 0.5, hjust=1)
  )


plot_discontinuation_sex <- cprd_dataset %>%
  group_by(drugclass, gender) %>%
  summarise(
    `By 3-months` = sum(as.numeric(stopdrug_3m_6mFU)-1, na.rm = TRUE)/n(),
    `By 6-months` = sum(as.numeric(stopdrug_6m_6mFU)-1, na.rm = TRUE)/n()
  ) %>%
  ungroup() %>%
  gather("Discontinuation", "value", -drugclass, -gender) %>%
  mutate(
    category = paste(drugclass, Discontinuation, sep = "_"),
    category = factor(category, levels = rev(c("MFN_By 3-months", "MFN_By 6-months", 
                                               "DPP4_By 3-months", "DPP4_By 6-months", 
                                               "SGLT2_By 3-months", "SGLT2_By 6-months", 
                                               "GLP1_By 3-months", "GLP1_By 6-months", 
                                               "SU_By 3-months", "SU_By 6-months", 
                                               "TZD_By 3-months", "TZD_By 6-months"))
    ),
    gender = factor(gender, levels = c(1, 2), labels = c("Male", "Female"))
  ) %>%
  rename("Gender" = "gender") %>%
  ggplot() +
  geom_col(aes(y = category, x = value, fill = Gender), position = "dodge") +
  theme_bw() +
  scale_y_discrete("Treatments", labels = rev(c("Metformin", "", 
                                                "DPP4i", "", 
                                                "SGLT2i", "", 
                                                "GLP-1RA", "", 
                                                "SU", "", 
                                                "TZD", ""))) +
  scale_x_continuous("Discontinuation", labels = scales::percent) +
  theme(
    legend.position = "bottom",
    axis.ticks.y = element_blank()
  )



plot_discontinuation_bmi <- patchwork::wrap_plots(
  cprd_dataset %>%
    filter(drugclass == "MFN" & !is.na(stopdrug_3m_6mFU)) %>%
    mutate(Discontinuation = factor(stopdrug_3m_6mFU, levels = c(0, 1), labels = c("No", "Yes"))) %>%
    ggplot() +
    geom_density(aes(x = prebmi, fill = Discontinuation), alpha = 0.6) +
    lims(x = c(0, 110)) +
    labs(y = "MFN") +
    theme_minimal() +
    theme(axis.text = element_blank(),
          axis.ticks = element_blank(),
          axis.title.x = element_blank()),
  cprd_dataset %>%
    filter(drugclass == "MFN" & !is.na(stopdrug_6m_6mFU)) %>%
    mutate(Discontinuation = factor(stopdrug_6m_6mFU, levels = c(0, 1), labels = c("No", "Yes"))) %>%
    ggplot() +
    geom_density(aes(x = prebmi, fill = Discontinuation), alpha = 0.6) +
    lims(x = c(0, 110)) +
    labs(y = "") +
    theme_minimal() +
    theme(axis.text = element_blank(),
          axis.ticks = element_blank(),
          axis.title.x = element_blank()),
  cprd_dataset %>%
    filter(drugclass == "DPP4" & !is.na(stopdrug_3m_6mFU)) %>%
    mutate(Discontinuation = factor(stopdrug_3m_6mFU, levels = c(0, 1), labels = c("No", "Yes"))) %>%
    ggplot() +
    geom_density(aes(x = prebmi, fill = Discontinuation), alpha = 0.6) +
    lims(x = c(0, 110)) +
    labs(y = "DPP4") +
    theme_minimal() +
    theme(axis.text = element_blank(),
          axis.ticks = element_blank(),
          axis.title.x = element_blank()),
  cprd_dataset %>%
    filter(drugclass == "DPP4" & !is.na(stopdrug_6m_6mFU)) %>%
    mutate(Discontinuation = factor(stopdrug_6m_6mFU, levels = c(0, 1), labels = c("No", "Yes"))) %>%
    ggplot() +
    geom_density(aes(x = prebmi, fill = Discontinuation), alpha = 0.6) +
    lims(x = c(0, 110)) +
    labs(y = "") +
    theme_minimal() +
    theme(axis.text = element_blank(),
          axis.ticks = element_blank(),
          axis.title.x = element_blank()),
  cprd_dataset %>%
    filter(drugclass == "SGLT2" & !is.na(stopdrug_3m_6mFU)) %>%
    mutate(Discontinuation = factor(stopdrug_3m_6mFU, levels = c(0, 1), labels = c("No", "Yes"))) %>%
    ggplot() +
    geom_density(aes(x = prebmi, fill = Discontinuation), alpha = 0.6) +
    lims(x = c(0, 110)) +
    labs(y = "SGLT2") +
    theme_minimal() +
    theme(axis.text = element_blank(),
          axis.ticks = element_blank(),
          axis.title.x = element_blank()),
  cprd_dataset %>%
    filter(drugclass == "SGLT2" & !is.na(stopdrug_6m_6mFU)) %>%
    mutate(Discontinuation = factor(stopdrug_6m_6mFU, levels = c(0, 1), labels = c("No", "Yes"))) %>%
    ggplot() +
    geom_density(aes(x = prebmi, fill = Discontinuation), alpha = 0.6) +
    lims(x = c(0, 110)) +
    labs(y = "") +
    theme_minimal() +
    theme(axis.text = element_blank(),
          axis.ticks = element_blank(),
          axis.title.x = element_blank()),
  cprd_dataset %>%
    filter(drugclass == "GLP1" & !is.na(stopdrug_3m_6mFU)) %>%
    mutate(Discontinuation = factor(stopdrug_3m_6mFU, levels = c(0, 1), labels = c("No", "Yes"))) %>%
    ggplot() +
    geom_density(aes(x = prebmi, fill = Discontinuation), alpha = 0.6) +
    lims(x = c(0, 110)) +
    labs(y = "GLP1") +
    theme_minimal() +
    theme(axis.text = element_blank(),
          axis.ticks = element_blank(),
          axis.title.x = element_blank()),
  cprd_dataset %>%
    filter(drugclass == "GLP1" & !is.na(stopdrug_6m_6mFU)) %>%
    mutate(Discontinuation = factor(stopdrug_6m_6mFU, levels = c(0, 1), labels = c("No", "Yes"))) %>%
    ggplot() +
    geom_density(aes(x = prebmi, fill = Discontinuation), alpha = 0.6) +
    lims(x = c(0, 110)) +
    labs(y = "") +
    theme_minimal() +
    theme(axis.text = element_blank(),
          axis.ticks = element_blank(),
          axis.title.x = element_blank()),
  cprd_dataset %>%
    filter(drugclass == "SU" & !is.na(stopdrug_3m_6mFU)) %>%
    mutate(Discontinuation = factor(stopdrug_3m_6mFU, levels = c(0, 1), labels = c("No", "Yes"))) %>%
    ggplot() +
    geom_density(aes(x = prebmi, fill = Discontinuation), alpha = 0.6) +
    lims(x = c(0, 110)) +
    labs(y = "SU") +
    theme_minimal() +
    theme(axis.text = element_blank(),
          axis.ticks = element_blank(),
          axis.title.x = element_blank()),
  cprd_dataset %>%
    filter(drugclass == "SU" & !is.na(stopdrug_6m_6mFU)) %>%
    mutate(Discontinuation = factor(stopdrug_6m_6mFU, levels = c(0, 1), labels = c("No", "Yes"))) %>%
    ggplot() +
    geom_density(aes(x = prebmi, fill = Discontinuation), alpha = 0.6) +
    lims(x = c(0, 110)) +
    labs(y = "") +
    theme_minimal() +
    theme(axis.text = element_blank(),
          axis.ticks = element_blank(),
          axis.title.x = element_blank()),
  cprd_dataset %>%
    filter(drugclass == "TZD" & !is.na(stopdrug_3m_6mFU)) %>%
    mutate(Discontinuation = factor(stopdrug_3m_6mFU, levels = c(0, 1), labels = c("No", "Yes"))) %>%
    ggplot() +
    geom_density(aes(x = prebmi, fill = Discontinuation), alpha = 0.6) +
    lims(x = c(0, 110)) +
    labs(y = "TZD") +
    theme_minimal() +
    theme(axis.text = element_blank(),
          axis.ticks = element_blank(),
          axis.title.x = element_blank()),
  cprd_dataset %>%
    filter(drugclass == "TZD" & !is.na(stopdrug_6m_6mFU)) %>%
    mutate(Discontinuation = factor(stopdrug_6m_6mFU, levels = c(0, 1), labels = c("No", "Yes"))) %>%
    ggplot() +
    geom_density(aes(x = prebmi, fill = Discontinuation), alpha = 0.6) +
    lims(x = c(0, 110)) +
    labs(y = "", x = "BMI") +
    theme_minimal() +
    theme(axis.text.y = element_blank(),
          axis.ticks.y = element_blank()),
  ncol = 1
) +
  plot_layout(guides = 'collect') & 
  theme(
    legend.position = "bottom",
    axis.title.y = element_text(angle = 0, vjust = 0.5, hjust=1)
  )




plot_discontinuation_hba1c <- patchwork::wrap_plots(
  cprd_dataset %>%
    filter(drugclass == "MFN" & !is.na(stopdrug_3m_6mFU)) %>%
    mutate(Discontinuation = factor(stopdrug_3m_6mFU, levels = c(0, 1), labels = c("No", "Yes"))) %>%
    ggplot() +
    geom_density(aes(x = prehba1c, fill = Discontinuation), alpha = 0.6) +
    lims(x = c(20, 200)) +
    labs(y = "MFN") +
    theme_minimal() +
    theme(axis.text = element_blank(),
          axis.ticks = element_blank(),
          axis.title.x = element_blank()),
  cprd_dataset %>%
    filter(drugclass == "MFN" & !is.na(stopdrug_6m_6mFU)) %>%
    mutate(Discontinuation = factor(stopdrug_6m_6mFU, levels = c(0, 1), labels = c("No", "Yes"))) %>%
    ggplot() +
    geom_density(aes(x = prehba1c, fill = Discontinuation), alpha = 0.6) +
    lims(x = c(20, 200)) +
    labs(y = "") +
    theme_minimal() +
    theme(axis.text = element_blank(),
          axis.ticks = element_blank(),
          axis.title.x = element_blank()),
  cprd_dataset %>%
    filter(drugclass == "DPP4" & !is.na(stopdrug_3m_6mFU)) %>%
    mutate(Discontinuation = factor(stopdrug_3m_6mFU, levels = c(0, 1), labels = c("No", "Yes"))) %>%
    ggplot() +
    geom_density(aes(x = prehba1c, fill = Discontinuation), alpha = 0.6) +
    lims(x = c(20, 200)) +
    labs(y = "DPP4") +
    theme_minimal() +
    theme(axis.text = element_blank(),
          axis.ticks = element_blank(),
          axis.title.x = element_blank()),
  cprd_dataset %>%
    filter(drugclass == "DPP4" & !is.na(stopdrug_6m_6mFU)) %>%
    mutate(Discontinuation = factor(stopdrug_6m_6mFU, levels = c(0, 1), labels = c("No", "Yes"))) %>%
    ggplot() +
    geom_density(aes(x = prehba1c, fill = Discontinuation), alpha = 0.6) +
    lims(x = c(20, 200)) +
    labs(y = "") +
    theme_minimal() +
    theme(axis.text = element_blank(),
          axis.ticks = element_blank(),
          axis.title.x = element_blank()),
  cprd_dataset %>%
    filter(drugclass == "SGLT2" & !is.na(stopdrug_3m_6mFU)) %>%
    mutate(Discontinuation = factor(stopdrug_3m_6mFU, levels = c(0, 1), labels = c("No", "Yes"))) %>%
    ggplot() +
    geom_density(aes(x = prehba1c, fill = Discontinuation), alpha = 0.6) +
    lims(x = c(20, 200)) +
    labs(y = "SGLT2") +
    theme_minimal() +
    theme(axis.text = element_blank(),
          axis.ticks = element_blank(),
          axis.title.x = element_blank()),
  cprd_dataset %>%
    filter(drugclass == "SGLT2" & !is.na(stopdrug_6m_6mFU)) %>%
    mutate(Discontinuation = factor(stopdrug_6m_6mFU, levels = c(0, 1), labels = c("No", "Yes"))) %>%
    ggplot() +
    geom_density(aes(x = prehba1c, fill = Discontinuation), alpha = 0.6) +
    lims(x = c(20, 200)) +
    labs(y = "") +
    theme_minimal() +
    theme(axis.text = element_blank(),
          axis.ticks = element_blank(),
          axis.title.x = element_blank()),
  cprd_dataset %>%
    filter(drugclass == "GLP1" & !is.na(stopdrug_3m_6mFU)) %>%
    mutate(Discontinuation = factor(stopdrug_3m_6mFU, levels = c(0, 1), labels = c("No", "Yes"))) %>%
    ggplot() +
    geom_density(aes(x = prehba1c, fill = Discontinuation), alpha = 0.6) +
    lims(x = c(20, 200)) +
    labs(y = "GLP1") +
    theme_minimal() +
    theme(axis.text = element_blank(),
          axis.ticks = element_blank(),
          axis.title.x = element_blank()),
  cprd_dataset %>%
    filter(drugclass == "GLP1" & !is.na(stopdrug_6m_6mFU)) %>%
    mutate(Discontinuation = factor(stopdrug_6m_6mFU, levels = c(0, 1), labels = c("No", "Yes"))) %>%
    ggplot() +
    geom_density(aes(x = prehba1c, fill = Discontinuation), alpha = 0.6) +
    lims(x = c(20, 200)) +
    labs(y = "") +
    theme_minimal() +
    theme(axis.text = element_blank(),
          axis.ticks = element_blank(),
          axis.title.x = element_blank()),
  cprd_dataset %>%
    filter(drugclass == "SU" & !is.na(stopdrug_3m_6mFU)) %>%
    mutate(Discontinuation = factor(stopdrug_3m_6mFU, levels = c(0, 1), labels = c("No", "Yes"))) %>%
    ggplot() +
    geom_density(aes(x = prehba1c, fill = Discontinuation), alpha = 0.6) +
    lims(x = c(20, 200)) +
    labs(y = "SU") +
    theme_minimal() +
    theme(axis.text = element_blank(),
          axis.ticks = element_blank(),
          axis.title.x = element_blank()),
  cprd_dataset %>%
    filter(drugclass == "SU" & !is.na(stopdrug_6m_6mFU)) %>%
    mutate(Discontinuation = factor(stopdrug_6m_6mFU, levels = c(0, 1), labels = c("No", "Yes"))) %>%
    ggplot() +
    geom_density(aes(x = prehba1c, fill = Discontinuation), alpha = 0.6) +
    lims(x = c(20, 200)) +
    labs(y = "") +
    theme_minimal() +
    theme(axis.text = element_blank(),
          axis.ticks = element_blank(),
          axis.title.x = element_blank()),
  cprd_dataset %>%
    filter(drugclass == "TZD" & !is.na(stopdrug_3m_6mFU)) %>%
    mutate(Discontinuation = factor(stopdrug_3m_6mFU, levels = c(0, 1), labels = c("No", "Yes"))) %>%
    ggplot() +
    geom_density(aes(x = prehba1c, fill = Discontinuation), alpha = 0.6) +
    lims(x = c(20, 200)) +
    labs(y = "TZD") +
    theme_minimal() +
    theme(axis.text = element_blank(),
          axis.ticks = element_blank(),
          axis.title.x = element_blank()),
  cprd_dataset %>%
    filter(drugclass == "TZD" & !is.na(stopdrug_6m_6mFU)) %>%
    mutate(Discontinuation = factor(stopdrug_6m_6mFU, levels = c(0, 1), labels = c("No", "Yes"))) %>%
    ggplot() +
    geom_density(aes(x = prehba1c, fill = Discontinuation), alpha = 0.6) +
    lims(x = c(20, 200)) +
    labs(y = "", x = "HbA1c") +
    theme_minimal() +
    theme(axis.text.y = element_blank(),
          axis.ticks.y = element_blank()),
  ncol = 1
) +
  plot_layout(guides = 'collect') & 
  theme(
    legend.position = "bottom",
    axis.title.y = element_text(angle = 0, vjust = 0.5, hjust=1)
  )



plot_discontinuation_dm_dur <- patchwork::wrap_plots(
  cprd_dataset %>%
    filter(drugclass == "MFN" & !is.na(stopdrug_3m_6mFU)) %>%
    mutate(Discontinuation = factor(stopdrug_3m_6mFU, levels = c(0, 1), labels = c("No", "Yes"))) %>%
    ggplot() +
    geom_density(aes(x = dstartdate_dm_dur, fill = Discontinuation), alpha = 0.6) +
    lims(x = c(0,30)) +
    labs(y = "MFN") +
    theme_minimal() +
    theme(axis.text = element_blank(),
          axis.ticks = element_blank(),
          axis.title.x = element_blank()),
  cprd_dataset %>%
    filter(drugclass == "MFN" & !is.na(stopdrug_6m_6mFU)) %>%
    mutate(Discontinuation = factor(stopdrug_6m_6mFU, levels = c(0, 1), labels = c("No", "Yes"))) %>%
    ggplot() +
    geom_density(aes(x = dstartdate_dm_dur, fill = Discontinuation), alpha = 0.6) +
    lims(x = c(0,30)) +
    labs(y = "") +
    theme_minimal() +
    theme(axis.text = element_blank(),
          axis.ticks = element_blank(),
          axis.title.x = element_blank()),
  cprd_dataset %>%
    filter(drugclass == "DPP4" & !is.na(stopdrug_3m_6mFU)) %>%
    mutate(Discontinuation = factor(stopdrug_3m_6mFU, levels = c(0, 1), labels = c("No", "Yes"))) %>%
    ggplot() +
    geom_density(aes(x = dstartdate_dm_dur, fill = Discontinuation), alpha = 0.6) +
    lims(x = c(0,30)) +
    labs(y = "DPP4") +
    theme_minimal() +
    theme(axis.text = element_blank(),
          axis.ticks = element_blank(),
          axis.title.x = element_blank()),
  cprd_dataset %>%
    filter(drugclass == "DPP4" & !is.na(stopdrug_6m_6mFU)) %>%
    mutate(Discontinuation = factor(stopdrug_6m_6mFU, levels = c(0, 1), labels = c("No", "Yes"))) %>%
    ggplot() +
    geom_density(aes(x = dstartdate_dm_dur, fill = Discontinuation), alpha = 0.6) +
    lims(x = c(0,30)) +
    labs(y = "") +
    theme_minimal() +
    theme(axis.text = element_blank(),
          axis.ticks = element_blank(),
          axis.title.x = element_blank()),
  cprd_dataset %>%
    filter(drugclass == "SGLT2" & !is.na(stopdrug_3m_6mFU)) %>%
    mutate(Discontinuation = factor(stopdrug_3m_6mFU, levels = c(0, 1), labels = c("No", "Yes"))) %>%
    ggplot() +
    geom_density(aes(x = dstartdate_dm_dur, fill = Discontinuation), alpha = 0.6) +
    lims(x = c(0,30)) +
    labs(y = "SGLT2") +
    theme_minimal() +
    theme(axis.text = element_blank(),
          axis.ticks = element_blank(),
          axis.title.x = element_blank()),
  cprd_dataset %>%
    filter(drugclass == "SGLT2" & !is.na(stopdrug_6m_6mFU)) %>%
    mutate(Discontinuation = factor(stopdrug_6m_6mFU, levels = c(0, 1), labels = c("No", "Yes"))) %>%
    ggplot() +
    geom_density(aes(x = dstartdate_dm_dur, fill = Discontinuation), alpha = 0.6) +
    lims(x = c(0,30)) +
    labs(y = "") +
    theme_minimal() +
    theme(axis.text = element_blank(),
          axis.ticks = element_blank(),
          axis.title.x = element_blank()),
  cprd_dataset %>%
    filter(drugclass == "GLP1" & !is.na(stopdrug_3m_6mFU)) %>%
    mutate(Discontinuation = factor(stopdrug_3m_6mFU, levels = c(0, 1), labels = c("No", "Yes"))) %>%
    ggplot() +
    geom_density(aes(x = dstartdate_dm_dur, fill = Discontinuation), alpha = 0.6) +
    lims(x = c(0,30)) +
    labs(y = "GLP1") +
    theme_minimal() +
    theme(axis.text = element_blank(),
          axis.ticks = element_blank(),
          axis.title.x = element_blank()),
  cprd_dataset %>%
    filter(drugclass == "GLP1" & !is.na(stopdrug_6m_6mFU)) %>%
    mutate(Discontinuation = factor(stopdrug_6m_6mFU, levels = c(0, 1), labels = c("No", "Yes"))) %>%
    ggplot() +
    geom_density(aes(x = dstartdate_dm_dur, fill = Discontinuation), alpha = 0.6) +
    lims(x = c(0,30)) +
    labs(y = "") +
    theme_minimal() +
    theme(axis.text = element_blank(),
          axis.ticks = element_blank(),
          axis.title.x = element_blank()),
  cprd_dataset %>%
    filter(drugclass == "SU" & !is.na(stopdrug_3m_6mFU)) %>%
    mutate(Discontinuation = factor(stopdrug_3m_6mFU, levels = c(0, 1), labels = c("No", "Yes"))) %>%
    ggplot() +
    geom_density(aes(x = dstartdate_dm_dur, fill = Discontinuation), alpha = 0.6) +
    lims(x = c(0,30)) +
    labs(y = "SU") +
    theme_minimal() +
    theme(axis.text = element_blank(),
          axis.ticks = element_blank(),
          axis.title.x = element_blank()),
  cprd_dataset %>%
    filter(drugclass == "SU" & !is.na(stopdrug_6m_6mFU)) %>%
    mutate(Discontinuation = factor(stopdrug_6m_6mFU, levels = c(0, 1), labels = c("No", "Yes"))) %>%
    ggplot() +
    geom_density(aes(x = dstartdate_dm_dur, fill = Discontinuation), alpha = 0.6) +
    lims(x = c(0,30)) +
    labs(y = "") +
    theme_minimal() +
    theme(axis.text = element_blank(),
          axis.ticks = element_blank(),
          axis.title.x = element_blank()),
  cprd_dataset %>%
    filter(drugclass == "TZD" & !is.na(stopdrug_3m_6mFU)) %>%
    mutate(Discontinuation = factor(stopdrug_3m_6mFU, levels = c(0, 1), labels = c("No", "Yes"))) %>%
    ggplot() +
    geom_density(aes(x = dstartdate_dm_dur, fill = Discontinuation), alpha = 0.6) +
    lims(x = c(0,30)) +
    labs(y = "TZD") +
    theme_minimal() +
    theme(axis.text = element_blank(),
          axis.ticks = element_blank(),
          axis.title.x = element_blank()),
  cprd_dataset %>%
    filter(drugclass == "TZD" & !is.na(stopdrug_6m_6mFU)) %>%
    mutate(Discontinuation = factor(stopdrug_6m_6mFU, levels = c(0, 1), labels = c("No", "Yes"))) %>%
    ggplot() +
    geom_density(aes(x = dstartdate_dm_dur, fill = Discontinuation), alpha = 0.6) +
    lims(x = c(0,30)) +
    labs(y = "", x = "Diabetes duration") +
    theme_minimal() +
    theme(axis.text.y = element_blank(),
          axis.ticks.y = element_blank()),
  ncol = 1
) +
  plot_layout(guides = 'collect') & 
  theme(
    legend.position = "bottom",
    axis.title.y = element_text(angle = 0, vjust = 0.5, hjust=1)
  )

output_path <- "results/figures/"


pdf(paste0(output_path, "01.discontinuation_strata.pdf"), width = 7, height = 8)
plot_discontinuation_drugclass
plot_discontinuation_sex
plot_discontinuation_age
plot_discontinuation_bmi
plot_discontinuation_hba1c
plot_discontinuation_dm_dur
dev.off()

