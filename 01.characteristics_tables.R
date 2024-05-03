######################################################################
##
##  In this file, we create characteristics tables for all therapies
##
######################################################################

# load working directory
setwd("Samples/T2D_Discontinuation")

# load functions
source("code/00.set_up_data.R")

# load libraries
library(ggplot2)
library(tableone)
library(naniar)
library(patchwork)

# load dataset
cprd_dataset <- set_up_data(
  raw_data = "20240308_t2d_1stinstance",
  diagnosis = TRUE,
  therapies = c("DPP4", "GLP1", "MFN", "SGLT2", "SU", "TZD"),
  dataset = "full.dataset"
)


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

pdf("results/figures/plot_initiations_by_year.pdf", width = 11, height = 6)
plot_initiations_by_year
dev.off()


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

pdf("results/figures/plot_discontinuation_by_year.pdf", width = 11, height = 6)
plot_discontinuation_by_year
dev.off()



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

pdf("results/figures/plot_discontinuation_by_drugclass.pdf", width = 11, height = 6)
plot_discontinuation_by_drugclass
dev.off()





#########################################################################################################
#########################################################################################################
#
# Table of Characteristics
#
#########################################################################################################
#########################################################################################################

output_path <- "results/tables"


cprd_tables <- cprd_dataset %>%
  mutate(
    prealt_na = ifelse(!is.na(prealt), "No", NA),
    pretotalcholesterol_na = ifelse(!is.na(pretotalcholesterol), "No", NA),
    preegfr_na = ifelse(!is.na(preegfr), "No", NA),
    prebmi_na = ifelse(!is.na(prebmi), "No", NA),
    dstartdate_age_na = ifelse(!is.na(dstartdate_age), "No", NA),
    dstartdate_dm_dur_na = ifelse(!is.na(dstartdate_dm_dur), "No", NA)
  )


vars <- c(
  # Outcome
  "stopdrug_3m_6mFU",
  # Drug taken
  "drugclass",
  # Extra info
  "dstartdate_dm_dur", "dstartdate_dm_dur_na", "dstartdate_age", "dstartdate_age_na", "drugline", "numdrugs", "smoking_cat", "imd2015_10",
  "predrug_statins", "stopdrug_3m_3mFU_MFN_hist", "ethnicity_5cat", "gender", "predrug_bloodmed",
  # Biomarkers
  "prehba1c", "preegfr", "preegfr_na", "prebmi", "prebmi_na", "prealt", "prealt_na", 
  "pretotalcholesterol", "pretotalcholesterol_na",
  # preldl, prehdl, pretriglyceride,
  # Comorbidities
  ## Hist of cardiovascular
  "predrug_cardio_event",
  "predrug_angina", "predrug_myocardialinfarction", "predrug_ihd", "predrug_pad",
  "predrug_revasc", "predrug_stroke", 
  ## Heart problems
  "predrug_heart_event",
  "predrug_heartfailure", "predrug_hypertension",
  ### Microvascular
  "predrug_micro_event",
  "predrug_retinopathy", "predrug_diabeticnephropathy", "predrug_neuropathy",
  ## CKD
  "preckdstage", 
  ## CLD
  "predrug_cld"
)



vars_cat <- c(
  # Outcome
  "stopdrug_3m_6mFU",
  # Extra info
  "dstartdate_dm_dur_na", "dstartdate_age_na", "smoking_cat", "imd2015_10",
  "predrug_statins", "ethnicity_5cat", "gender", "predrug_bloodmed",
  # Biomarkers
  "preegfr_na", "prebmi_na", "prealt_na", 
  "pretotalcholesterol_na",
  # Comorbidities
  ## Hist of cardiovascular
  "predrug_cardio_event",
  "predrug_angina", "predrug_myocardialinfarction", "predrug_ihd", "predrug_pad",
  "predrug_revasc", "predrug_stroke",
  ## Heart problems
  "predrug_heart_event",
  "predrug_heartfailure", "predrug_hypertension",
  ### Microvascular
  "predrug_micro_event",
  "predrug_retinopathy", "predrug_diabeticnephropathy", "predrug_neuropathy",
  ## CKD
  "preckdstage", 
  ## CLD
  "predrug_cld"
)

### Table by drugclass

table_characteristics_drugs <- CreateTableOne(
  vars = vars,
  factorVars = vars_cat,
  includeNA = TRUE,
  strata = c("drugclass"),
  data = cprd_tables,
  test = TRUE
)

table_characteristics_drugs_print <- print(table_characteristics_drugs, exact = "stage", quote = FALSE, noSpaces = TRUE, printToggle = FALSE)

write.csv(table_characteristics_drugs_print, file = paste0(output_path, "/table_characteristics_drugs.csv"))


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


### DPP4

table_characteristics_dpp4 <- CreateTableOne(
  vars = vars,
  factorVars = vars_cat,
  includeNA = TRUE,
  strata = c("stopdrug_3m_6mFU", "drugclass"),
  data = cprd_tables %>% filter(drugclass == "DPP4"),
  test = TRUE
)

table_characteristics_dpp4_print <- print(table_characteristics_dpp4, exact = "stage", quote = FALSE, noSpaces = TRUE, printToggle = FALSE)

write.csv(table_characteristics_dpp4_print, file = paste0(output_path, "/table_characteristics_dpp4.csv"))

### GLP1

table_characteristics_glp1 <- CreateTableOne(
  vars = vars,
  factorVars = vars_cat,
  includeNA = TRUE,
  strata = c("stopdrug_3m_6mFU", "drugclass"),
  data = cprd_tables %>% filter(drugclass == "GLP1"),
  test = TRUE
)

table_characteristics_glp1_print <- print(table_characteristics_glp1, exact = "stage", quote = FALSE, noSpaces = TRUE, printToggle = FALSE)

write.csv(table_characteristics_glp1_print, file = paste0(output_path, "/table_characteristics_glp1.csv"))

### MFN

table_characteristics_mfn <- CreateTableOne(
  vars = vars,
  factorVars = vars_cat,
  includeNA = TRUE,
  strata = c("stopdrug_3m_6mFU", "drugclass"),
  data = cprd_tables %>% filter(drugclass == "MFN"),
  test = TRUE
)

table_characteristics_mfn_print <- print(table_characteristics_mfn, exact = "stage", quote = FALSE, noSpaces = TRUE, printToggle = FALSE)

write.csv(table_characteristics_mfn_print, file = paste0(output_path, "/table_characteristics_mfn.csv"))

### SGLT2

table_characteristics_sglt2 <- CreateTableOne(
  vars = vars,
  factorVars = vars_cat,
  includeNA = TRUE,
  strata = c("stopdrug_3m_6mFU", "drugclass"),
  data = cprd_tables %>% filter(drugclass == "SGLT2"),
  test = TRUE
)

table_characteristics_sglt2_print <- print(table_characteristics_sglt2, exact = "stage", quote = FALSE, noSpaces = TRUE, printToggle = FALSE)

write.csv(table_characteristics_sglt2_print, file = paste0(output_path, "/table_characteristics_sglt2.csv"))

### SU

table_characteristics_su <- CreateTableOne(
  vars = vars,
  factorVars = vars_cat,
  includeNA = TRUE,
  strata = c("stopdrug_3m_6mFU", "drugclass"),
  data = cprd_tables %>% filter(drugclass == "SU"),
  test = TRUE
)

table_characteristics_su_print <- print(table_characteristics_su, exact = "stage", quote = FALSE, noSpaces = TRUE, printToggle = FALSE)

write.csv(table_characteristics_su_print, file = paste0(output_path, "/table_characteristics_su.csv"))

### TZD

table_characteristics_tzd <- CreateTableOne(
  vars = vars,
  factorVars = vars_cat,
  includeNA = TRUE,
  strata = c("stopdrug_3m_6mFU", "drugclass"),
  data = cprd_tables %>% filter(drugclass == "TZD"),
  test = TRUE
)

table_characteristics_tzd_print <- print(table_characteristics_tzd, exact = "stage", quote = FALSE, noSpaces = TRUE, printToggle = FALSE)

write.csv(table_characteristics_tzd_print, file = paste0(output_path, "/table_characteristics_tzd.csv"))






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


pdf(paste0(output_path, "discontinuation_strata.pdf"), width = 7, height = 8)
plot_discontinuation_drugclass
plot_discontinuation_sex
plot_discontinuation_age
plot_discontinuation_bmi
plot_discontinuation_hba1c
plot_discontinuation_dm_dur
dev.off()



### Missingness investigation

pdf(paste0(output_path, "missingness_investigation.pdf"), width = 7, height = 12)

cprd_dataset %>%
  select(
    drugclass,
    # Biomarkers
    precreatinine_blood, 
    prealt, pretotalcholesterol, predbp, presbp, prehba1c, 
    preegfr, prebilirubin,
    # Commorbidities
    preckdstage, predrug_frailty_mild, predrug_frailty_moderate, 
    predrug_frailty_severe, predrug_primary_hhf, predrug_af, predrug_angina, 
    predrug_asthma, predrug_bronchiectasis, predrug_cld, predrug_copd, 
    predrug_cysticfibrosis, predrug_dementia, predrug_diabeticnephropathy, 
    predrug_fh_premature_cvd, 
    predrug_haem_cancer, predrug_heartfailure, 
    predrug_hypertension, predrug_ihd, predrug_myocardialinfarction, 
    predrug_neuropathy, predrug_otherneuroconditions, predrug_pad, 
    predrug_pulmonaryfibrosis, predrug_pulmonaryhypertension, 
    predrug_retinopathy, predrug_revasc, predrug_rheumatoidarthritis, 
    predrug_solid_cancer, predrug_solidorgantransplant, predrug_stroke, 
    predrug_tia, predrug_anxiety_disorders, predrug_medspecific_gi,
    predrug_benignprostatehyperplasia, predrug_micturition_control,
    predrug_volume_depletion, predrug_urinary_frequency, predrug_falls,
    predrug_lowerlimbfracture, predrug_incident_mi, predrug_incident_stroke,
    predrug_dka, predrug_osteoporosis, predrug_unstableangina, 
    predrug_amputation,
    hosp_admission_prev_year, hosp_admission_prev_year_count,
    # Extra info
    gender, prac_region, ethnicity_5cat, imd2015_10, dm_diag_age,
    ins_in_1_year, prebmi, smoking_cat, drugline, numdrugs,
    alcohol_cat, dstartdate_age, dstartdate_dm_dur, dstartmonth,
    CCI_index
  ) %>%
  gg_miss_upset()

cprd_dataset %>%
  filter(drugclass == "MFN") %>%
  gather(key, value) %>%
  group_by(key) %>%
  count(na = is.na(value)) %>%
  pivot_wider(names_from = na, values_from = n, values_fill = 0) %>%
  mutate(pct_missing = (`TRUE`/sum(`TRUE`, `FALSE`))*100) %>%
  ungroup() %>% 
  mutate(Present = 100 - pct_missing) %>% 
  gather(Key, value, 4:5) %>% 
  mutate(Key = recode(Key, pct_missing = "Missing")) %>% 
  ggplot(aes(x = key, y = value, fill = Key)) +
  geom_col(alpha = 0.85) +
  scale_fill_manual(name = "", 
                    values = c('tomato3', 'steelblue'), 
                    labels = c("Missing", "Present")) +
  coord_flip() +
  ggtitle("MFN") +
  labs(x = NULL, y = "Missing (%)")

cprd_dataset %>%
  filter(drugclass == "DPP4") %>%
  gather(key, value) %>%
  group_by(key) %>%
  count(na = is.na(value)) %>%
  pivot_wider(names_from = na, values_from = n, values_fill = 0) %>%
  mutate(pct_missing = (`TRUE`/sum(`TRUE`, `FALSE`))*100) %>%
  ungroup() %>% 
  mutate(Present = 100 - pct_missing) %>% 
  gather(Key, value, 4:5) %>% 
  mutate(Key = recode(Key, pct_missing = "Missing")) %>% 
  ggplot(aes(x = key, y = value, fill = Key)) +
  geom_col(alpha = 0.85) +
  scale_fill_manual(name = "", 
                    values = c('tomato3', 'steelblue'), 
                    labels = c("Missing", "Present")) +
  coord_flip() +
  ggtitle("DPP4") +
  labs(x = NULL, y = "Missing (%)")

cprd_dataset %>%
  filter(drugclass == "SGLT2") %>%
  gather(key, value) %>%
  group_by(key) %>%
  count(na = is.na(value)) %>%
  pivot_wider(names_from = na, values_from = n, values_fill = 0) %>%
  mutate(pct_missing = (`TRUE`/sum(`TRUE`, `FALSE`))*100) %>%
  ungroup() %>% 
  mutate(Present = 100 - pct_missing) %>% 
  gather(Key, value, 4:5) %>% 
  mutate(Key = recode(Key, pct_missing = "Missing")) %>% 
  ggplot(aes(x = key, y = value, fill = Key)) +
  geom_col(alpha = 0.85) +
  scale_fill_manual(name = "", 
                    values = c('tomato3', 'steelblue'), 
                    labels = c("Missing", "Present")) +
  coord_flip() +
  ggtitle("SGLT2") +
  labs(x = NULL, y = "Missing (%)")

cprd_dataset %>%
  filter(drugclass == "GLP1") %>%
  gather(key, value) %>%
  group_by(key) %>%
  count(na = is.na(value)) %>%
  pivot_wider(names_from = na, values_from = n, values_fill = 0) %>%
  mutate(pct_missing = (`TRUE`/sum(`TRUE`, `FALSE`))*100) %>%
  ungroup() %>% 
  mutate(Present = 100 - pct_missing) %>% 
  gather(Key, value, 4:5) %>% 
  mutate(Key = recode(Key, pct_missing = "Missing")) %>% 
  ggplot(aes(x = key, y = value, fill = Key)) +
  geom_col(alpha = 0.85) +
  scale_fill_manual(name = "", 
                    values = c('tomato3', 'steelblue'), 
                    labels = c("Missing", "Present")) +
  coord_flip() +
  ggtitle("GLP1") +
  labs(x = NULL, y = "Missing (%)")

cprd_dataset %>%
  filter(drugclass == "SU") %>%
  gather(key, value) %>%
  group_by(key) %>%
  count(na = is.na(value)) %>%
  pivot_wider(names_from = na, values_from = n, values_fill = 0) %>%
  mutate(pct_missing = (`TRUE`/sum(`TRUE`, `FALSE`))*100) %>%
  ungroup() %>% 
  mutate(Present = 100 - pct_missing) %>% 
  gather(Key, value, 4:5) %>% 
  mutate(Key = recode(Key, pct_missing = "Missing")) %>% 
  ggplot(aes(x = key, y = value, fill = Key)) +
  geom_col(alpha = 0.85) +
  scale_fill_manual(name = "", 
                    values = c('tomato3', 'steelblue'), 
                    labels = c("Missing", "Present")) +
  coord_flip() +
  ggtitle("SU") +
  labs(x = NULL, y = "Missing (%)")

cprd_dataset %>%
  filter(drugclass == "TZD") %>%
  gather(key, value) %>%
  group_by(key) %>%
  count(na = is.na(value)) %>%
  pivot_wider(names_from = na, values_from = n, values_fill = 0) %>%
  mutate(pct_missing = (`TRUE`/sum(`TRUE`, `FALSE`))*100) %>%
  ungroup() %>% 
  mutate(Present = 100 - pct_missing) %>% 
  gather(Key, value, 4:5) %>% 
  mutate(Key = recode(Key, pct_missing = "Missing")) %>% 
  ggplot(aes(x = key, y = value, fill = Key)) +
  geom_col(alpha = 0.85) +
  scale_fill_manual(name = "", 
                    values = c('tomato3', 'steelblue'), 
                    labels = c("Missing", "Present")) +
  coord_flip() +
  ggtitle("TZD") +
  labs(x = NULL, y = "Missing (%)")
dev.off()





