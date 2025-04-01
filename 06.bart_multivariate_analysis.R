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
library(rms)



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

# load dataset
cprd_dataset.6m <- set_up_data(
  raw_data = "20240814_t2d_1stinstance",
  diagnosis = FALSE,
  therapies = c("GLP1", "DPP4", "SGLT2", "TZD", "SU"),
  dataset = "6m.disc.dataset",
  follow_up = "6-months",
  full_prescribing_history = TRUE
) %>%
  drop_na(-prehdl, -stopdrug_12m_6mFU)

# load dataset
cprd_dataset.12m <- set_up_data(
  raw_data = "20240814_t2d_1stinstance",
  diagnosis = FALSE,
  therapies = c("GLP1", "DPP4", "SGLT2", "TZD", "SU"),
  dataset = "12m.disc.dataset",
  follow_up = "6-months",
  full_prescribing_history = TRUE
) %>%
  drop_na(-prehdl)

# load dataset
cprd_dataset.dev <- set_up_data(
  raw_data = "20240814_t2d_1stinstance",
  diagnosis = FALSE,
  therapies = c("MFN", "GLP1", "DPP4", "SGLT2", "TZD", "SU"),
  dataset = "3m.disc.dataset.dev",
  follow_up = "6-months",
  full_prescribing_history = TRUE
) %>%
  drop_na(-prehdl, -stopdrug_6m_6mFU, -stopdrug_12m_6mFU)

# load dataset
cprd_dataset.val <- set_up_data(
  raw_data = "20240814_t2d_1stinstance",
  diagnosis = FALSE,
  therapies = c("MFN", "GLP1", "DPP4", "SGLT2", "TZD", "SU"),
  dataset = "3m.disc.dataset.val",
  follow_up = "6-months",
  full_prescribing_history = TRUE
) %>%
  drop_na(-prehdl, -stopdrug_6m_6mFU, -stopdrug_12m_6mFU)


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


# Variable importance (sudo)
variable_importance_cohort_3m <- cprd_dataset.3m %>%
  mutate(
    pred = bartmachine_full_model_3m$p_hat
  )

m1 <- rms::ols(pred ~ rcs(dstartdate_age, 3) + gender + imd2015_10 + rcs(prebmi, 3) + rcs(dstartdate_dm_dur, 3) + rcs(prehba1c, 3) + predrug_frailty_proxy + ethnicity_5cat + predrug_bloodmed + smoking_cat + predrug_statins + rcs(preegfr, 3) + stopdrug_3m_3mFU_MFN_hist, data = variable_importance_cohort_3m, x = TRUE, y = TRUE)

values <- plot(anova(m1), what = 'proportion R2')

plot_variable_important_3m <- as.data.frame(values) %>%
  rownames_to_column() %>%
  left_join(
    data.frame(
      rowname = c("dstartdate_age", "gender", "imd2015_10", "prebmi", "dstartdate_dm_dur", "prehba1c", "drugline", "predrug_frailty_proxy", "ethnicity_5cat", "numdrugs", "predrug_bloodmed", "smoking_cat", "predrug_statins", "preegfr", "stopdrug_3m_3mFU_MFN_hist"),
      plotname = c("Age group", "Sex", "Index of multiple deprivation", "BMI", "Duration of diabetes", "HbA1c", "Number of therapy classes ever prescribed", "Frailty", "Ethnicity", "Number of other current glucose-lowering therapies", " Blood pressure medication", "Smoking", "Statins", "eGFR", "History of Metformin discontinuation within 3-months")
    ), by = c("rowname")
  ) %>%
  select(-rowname) %>%
  mutate(plotname = factor(plotname),
         values = values * 100) %>%
  ggplot(aes(y = forcats::fct_reorder(plotname, values), x = values)) +
  geom_segment(aes(x = 0, xend = values, yend = forcats::fct_reorder(plotname, values)), linetype = "dashed") +
  geom_point(size = 2, colour = "black") +
  ggtitle("Relative importance for discontinuation rate") +
  xlab("Relative Importance (%)") +
  theme_bw() +
  theme(axis.text.y = element_text(angle = 45, face = "bold"),
        axis.title.y = element_blank())


pdf("results/figures/06.var_importance_3m.pdf", width = 8, height = 6)
plot_variable_important_3m +
  theme(
    axis.text = element_text(size = 14),
    plot.title = element_text(size = 15, hjust = 0),
    plot.margin = unit(c(0.2, 0.2, 3, 0.2), "cm")
  )
dev.off()



plot_pooled_3m.overall <- calibration_plot(
  data = cprd_dataset.3m %>% mutate(pred = bartmachine_full_model_3m$p_hat) %>% mutate(stopdrug_3m_6mFU = as.numeric(stopdrug_3m_6mFU)-1) %>% mutate(grouping = "Pooled"), 
  obs = "stopdrug_3m_6mFU", 
  pred = "pred",
  group = "grouping",
  nTiles = 10)$calibration_plot +
  scale_x_continuous("Predicted discontinuation (%)", limits = c(0, 0.3), breaks = seq(0, 1, 0.1)) +
  scale_y_continuous("Observed discontinuation (%)", limits = c(0, 0.3), breaks = seq(0, 1, 0.1)) +
  theme_bw() +
  scale_colour_manual(values = c("Pooled" = "black", "SGLT2" = "#E69F00", "GLP1" = "#56B4E9", "SU" = "#CC79A7", "DPP4" = "#0072B2", "TZD" = "#D55E00", "MFN" = "grey"), breaks = rev(c("Pooled", "GLP1", "DPP4", "SGLT2", "TZD", "SU")), labels = rev(c("Pooled", "GLP-1RA", "DPP4i", "SGLT2i", "TZD", "SU")), name = "Therapy", guide = guide_legend(reverse = TRUE, nrow = 1)) +
  theme(
    legend.position = "bottom"
  )

# mean predictions for each decile
cprd_dataset.3m %>% 
  mutate(pred = bartmachine_full_model_3m$p_hat) %>%
  mutate(interim_test = ntile(pred, 10)) %>%
  group_by(interim_test) %>%
  mutate(mean_interim = median(pred, na.rm = TRUE)) %>%
  ungroup() %>%
  select(mean_interim) %>%
  distinct()

plot_pooled_3m.drugs <- calibration_plot(
  data = cprd_dataset.3m %>% mutate(pred = bartmachine_full_model_3m$p_hat) %>% mutate(stopdrug_3m_6mFU = as.numeric(stopdrug_3m_6mFU)-1), 
  obs = "stopdrug_3m_6mFU", 
  pred = "pred", 
  group = "drugclass",
  nTiles = 10)$calibration_plot +
  scale_x_continuous("Predicted discontinuation (%)", limits = c(0, 0.7), breaks = seq(0, 0.7, 0.1)) +
  scale_y_continuous("Observed discontinuation (%)", limits = c(0, 0.7), breaks = seq(0, 0.7, 0.1)) +
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


interim_dataset <- cprd_dataset.3m %>%
  mutate(pred = bartmachine_full_model_3m$p_hat) %>%
  select(drugclass, stopdrug_3m_6mFU, pred)


for (i in unique(interim_dataset$drugclass)) {
  
  interim_dataset_drugclass <- interim_dataset %>%
    filter(drugclass == i)
  
  interim_number <- pROC::roc(response = interim_dataset_drugclass %>% select(stopdrug_3m_6mFU) %>% unlist(),
            predictor = interim_dataset_drugclass %>% select(pred) %>% unlist(), ci = TRUE)$ci
  
  print(paste(i, signif(interim_number, digits = 3)))
  
}

#
# roc_coords_3m_by_drugs <- pROC::ci(pROC::roc(unlist(interim_dataset %>% filter(drugclass == "SGLT2") %>%select(stopdrug_3m_6mFU)), unlist(interim_dataset %>% filter(drugclass == "SGLT2") %>% select(pred))), of = "thresholds") %>%
#   as.data.frame() %>%
#   mutate(Therapy = "SGLT2i") %>%
#   rbind(
#     pROC::ci(pROC::roc(unlist(interim_dataset %>% filter(drugclass == "DPP4") %>% select(stopdrug_3m_6mFU)), unlist(interim_dataset %>% filter(drugclass == "DPP4") %>% select(pred))), of = "thresholds") %>%
#       as.data.frame() %>%
#       mutate(Therapy = "DPP4i"),
#     pROC::ci(pROC::roc(unlist(interim_dataset %>% filter(drugclass == "GLP1") %>%select(stopdrug_3m_6mFU)), unlist(interim_dataset %>% filter(drugclass == "GLP1") %>% select(pred))), of = "thresholds") %>%
#       as.data.frame() %>%
#       mutate(Therapy = "GLP1"),
#     pROC::ci(pROC::roc(unlist(interim_dataset %>% filter(drugclass == "TZD") %>%select(stopdrug_3m_6mFU)), unlist(interim_dataset %>% filter(drugclass == "TZD") %>% select(pred))), of = "thresholds") %>%
#       as.data.frame() %>%
#       mutate(Therapy = "TZD"),
#     pROC::ci(pROC::roc(unlist(interim_dataset %>% filter(drugclass == "SU") %>%select(stopdrug_3m_6mFU)), unlist(interim_dataset %>% filter(drugclass == "SU") %>% select(pred))), of = "thresholds") %>%
#       as.data.frame() %>%
#       mutate(Therapy = "SU")
#   )
# 
# saveRDS(roc_coords_3m_by_drugs, "results/Models/bartmachine/roc_coords_3m_by_drugs.rds")

roc_coords_3m_by_drugs <- readRDS("results/Models/bartmachine/roc_coords_3m_by_drugs.rds")



ppv_coords_3m <- pROC::coords(pROC::roc(
    cprd_dataset.3m %>%
      select(stopdrug_3m_6mFU) %>%
      mutate(stopdrug_3m_6mFU = factor(stopdrug_3m_6mFU)) %>%
      unlist(),
    bartmachine_full_model_3m$p_hat
  ), ret=c("threshold", "specificity", "ppv", "precision", "recall"))

plot_test <- ppv_coords_3m %>%
  as.data.frame() %>%
  ggplot() +
  geom_path(aes(x = threshold, y = ppv))


pdf("results/figures/06.precision_recall_3m.pdf", width = 5, height = 5.5)
ppv_coords_3m %>%
  as.data.frame() %>%
  ggplot() +
  geom_path(aes(x = recall, y = precision, colour = threshold)) +
  scale_colour_continuous(type = "viridis", limits = c(0, 1)) +
  theme_bw() +
  labs(x = "Recall", y = "Precision", colour = "Threshold") +
  theme(
    legend.position = "bottom",
    axis.text = element_text(size = 15),
    axis.title = element_text(size = 17),
    legend.title = element_text(size = 13)
  )
dev.off()


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
  annotate("label", x = 0.60, y = 0.125, size = 5, label = paste0("AUROC=", signif(roc_values_3m_overall[2], digits = 3), "\n95%CI 0.608-0.616")) + 
  scale_x_continuous("1 - Specificity", limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
  scale_y_continuous("Sensitivity", limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
  theme_bw() +
  theme(
    legend.position = "bottom",
    axis.title = element_text(size = 18),
    axis.text = element_text(size = 16)
  )


plot_roc_3m.drugs <- roc_coords_3m_by_drugs %>%
  as.data.frame() %>%
  ggplot() +
  geom_abline(aes(intercept = 0, slope = 1), colour = "black", linetype = "dashed") +
  # geom_path(aes(x = 1 - specificity.2.5., y = sensitivity.2.5., colour = Therapy), alpha = 0.4) +
  geom_path(aes(x = 1 - specificity.50., y = sensitivity.50., colour = Therapy)) +
  # geom_path(aes(x = 1 - specificity.97.5., y = sensitivity.97.5., colour = Therapy), alpha = 0.4) +
  scale_x_continuous("1 - Specificity", limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
  scale_y_continuous("Sensitivity", limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
  scale_colour_manual(values = c("Pooled" = "black", "SGLT2i" = "#E69F00", "GLP1" = "#56B4E9", "SU" = "#CC79A7", "DPP4i" = "#0072B2", "TZD" = "#D55E00", "MFN" = "grey"), breaks = rev(c("Pooled", "GLP1", "DPP4i", "SGLT2i", "TZD", "SU")), labels = rev(c("Pooled", "GLP-1RA", "DPP4i", "SGLT2i", "TZD", "SU")), name = "Therapy", guide = guide_legend(reverse = TRUE, nrow = 1)) +
  theme_bw() +
  theme(
    legend.position = "bottom",
    axis.title = element_text(size = 18),
    axis.text = element_text(size = 16)
  )
  
  
  
pdf("results/figures/06.roc_cal_overall_by_drug.pdf", width = 10, height = 10)
patchwork::wrap_plots(
  # calibration curve
  plot_pooled_3m.overall +
    scale_x_continuous("Predicted discontinuation (%)", labels = scales::percent) +
    scale_y_continuous("Observed discontinuation (%)", labels = scales::percent) +
    theme_bw() +
    theme(
      legend.position = "none",
      axis.title = element_text(size = 18),
      axis.text = element_text(size = 16)
    ),
  # discrimination curve
  plot_roc_3m +
    theme_bw() +
    theme(
      legend.position = "none",
      axis.title = element_text(size = 18),
      axis.text = element_text(size = 20)
    ),
  # calibration curve
  plot_pooled_3m.drugs +
    scale_x_continuous("Predicted discontinuation (%)", labels = scales::percent) +
    scale_y_continuous("Observed discontinuation (%)", labels = scales::percent) +
    theme_bw() +
    theme(
      legend.position = "none",
      axis.title = element_text(size = 18),
      axis.text = element_text(size = 16)
    ),
  # discrimination curve
  plot_roc_3m.drugs +
    theme_bw() +
    theme(
      legend.position = "none",
      axis.title = element_text(size = 18),
      axis.text = element_text(size = 20)
    ),
  ggpubr::as_ggplot(get_legend(plot_pooled_3m.drugs))
) +
  patchwork::plot_layout(design = "
                         AB
                         CD
                         EE
                         ", heights = c(1,1,0.1)) +
  patchwork::plot_annotation(tag_levels = list(c("A.1", "A.2", "B.1", "B.2", "", "")))
dev.off()


pdf("results/figures/06.roc_3m_overall.pdf", width = 5, height = 5)

plot_roc_3m

dev.off()


pdf("results/figures/06.roc_cal_3m_overall.pdf", width = 10, height = 5)

patchwork::wrap_plots(
  # calibration curve
  plot_pooled_3m.overall +
    scale_x_continuous("Predicted", labels = scales::percent) +
    scale_y_continuous("Observed", labels = scales::percent) +
    theme_bw(),
  # discrimination curve
  plot_roc_3m +
    theme_bw()
) +
  patchwork::plot_layout(guides = 'collect') +
  patchwork::plot_annotation(tag_levels = list(c("A", "B"))) &
  theme(
    legend.position = "none",
    plot.tag = element_text(size = 16),
    axis.title = element_text(size = 18),
    axis.text = element_text(size = 16)
  )

dev.off()

pdf("results/figures/06.easd_24_roc_3m_overall.pdf", width = 5, height = 5)

roc_coords_3m %>%
  as.data.frame() %>%
  ggplot() +
  geom_abline(aes(intercept = 0, slope = 1), colour = "black", linetype = "dashed") +
  geom_path(aes(x = 1 - specificity.2.5., y = sensitivity.2.5.), alpha = 0.4) +
  geom_path(aes(x = 1 - specificity.50., y = sensitivity.50.)) +
  geom_path(aes(x = 1 - specificity.97.5., y = sensitivity.97.5.), alpha = 0.4) +
  annotate("label", x = 0.60, y = 0.125, size = 7, label = paste0("AUROC=", signif(roc_values_3m_overall[2], digits = 3), ")")) + 
  scale_x_continuous("1 - Specificity", limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
  scale_y_continuous("Sensitivity", limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
  theme_bw() +
  theme(
    legend.position = "bottom",
    axis.title = element_text(size = 18),
    axis.text = element_text(size = 20)
  )

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

# Variable importance (sudo)
variable_importance_cohort_6m <- cprd_dataset.6m %>%
  mutate(
    pred = bartmachine_full_model_6m$p_hat
  )

m1 <- rms::ols(pred ~ rcs(dstartdate_age, 3) + gender + imd2015_10 + rcs(prebmi, 3) + rcs(dstartdate_dm_dur, 3) + rcs(prehba1c, 3) + predrug_frailty_proxy + ethnicity_5cat + predrug_bloodmed + smoking_cat + predrug_statins + rcs(preegfr, 3) + stopdrug_3m_3mFU_MFN_hist, data = variable_importance_cohort_6m, x = TRUE, y = TRUE)

values <- plot(anova(m1), what = 'proportion R2')

plot_variable_important_6m <- as.data.frame(values) %>%
  rownames_to_column() %>%
  left_join(
    data.frame(
      rowname = c("dstartdate_age", "gender", "imd2015_10", "prebmi", "dstartdate_dm_dur", "prehba1c", "drugline", "predrug_frailty_proxy", "ethnicity_5cat", "numdrugs", "predrug_bloodmed", "smoking_cat", "predrug_statins", "preegfr", "stopdrug_3m_3mFU_MFN_hist"),
      plotname = c("Age group", "Sex", "Index of multiple deprivation", "BMI", "Duration of diabetes", "HbA1c", "Number of therapy classes ever prescribed", "Frailty", "Ethnicity", "Number of other current glucose-lowering therapies", " Blood pressure medication", "Smoking", "Statins", "eGFR", "History of Metformin discontinuation within 3-months")
    ), by = c("rowname")
  ) %>%
  select(-rowname) %>%
  mutate(plotname = factor(plotname),
         values = values * 100) %>%
  ggplot(aes(y = forcats::fct_reorder(plotname, values), x = values)) +
  geom_segment(aes(x = 0, xend = values, yend = forcats::fct_reorder(plotname, values)), linetype = "dashed") +
  geom_point(size = 2, colour = "black") +
  ggtitle("Relative importance for discontinuation rate") +
  xlab("Relative Importance (%)") +
  theme_bw() +
  theme(axis.text.y = element_text(angle = 45, face = "bold"),
        axis.title.y = element_blank())




pdf("results/figures/06.var_importance_6m.pdf", width = 8, height = 6)
plot_variable_important_6m +
  theme(
    axis.text = element_text(size = 14),
    plot.title = element_text(size = 15, hjust = 0),
    plot.margin = unit(c(0.2, 0.2, 3.5, 0.2), "cm")
  )
dev.off()



plot_pooled_6m.overall <- calibration_plot(
  data = cprd_dataset.6m %>% mutate(pred = bartmachine_full_model_6m$p_hat) %>% mutate(stopdrug_6m_6mFU = as.numeric(stopdrug_6m_6mFU)-1) %>% mutate(grouping = "Pooled"), 
  obs = "stopdrug_6m_6mFU", 
  pred = "pred",
  group = "grouping",
  nTiles = 10)$calibration_plot +
  scale_x_continuous("Predicted discontinuation (%)", limits = c(0, 0.7), breaks = seq(0, 0.7, 0.1)) +
  scale_y_continuous("Observed discontinuation (%)", limits = c(0, 0.7), breaks = seq(0, 0.7, 0.1)) +
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
  scale_x_continuous("Predicted discontinuation (%)", limits = c(0, 0.7), breaks = seq(0, 0.7, 0.1)) +
  scale_y_continuous("Observed discontinuation (%)", limits = c(0, 0.7), breaks = seq(0, 0.7, 0.1)) +
  theme_bw() +
  scale_colour_manual(values = c("Pooled" = "black", "SGLT2" = "#E69F00", "GLP1" = "#56B4E9", "SU" = "#CC79A7", "DPP4" = "#0072B2", "TZD" = "#D55E00", "MFN" = "grey"), breaks = rev(c("Pooled", "GLP1", "DPP4", "SGLT2", "TZD", "SU")), labels = rev(c("Pooled", "GLP-1RA", "DPP4i", "SGLT2i", "TZD", "SU")), name = "Therapy", guide = guide_legend(reverse = TRUE, nrow = 1)) +
  theme(
    legend.position = "bottom"
  )


roc_values_6m_overall <- pROC::roc(response = cprd_dataset.6m %>%
                                     select(stopdrug_6m_6mFU) %>%
                                     mutate(stopdrug_6m_6mFU = factor(stopdrug_6m_6mFU)) %>%
                                     unlist(),
                                   predictor = bartmachine_full_model_6m$p_hat, ci = TRUE)$ci %>%
  as.data.frame() %>%
  unlist()


interim_dataset <- cprd_dataset.6m %>%
  mutate(pred = bartmachine_full_model_6m$p_hat) %>%
  select(drugclass, stopdrug_6m_6mFU, pred)


for (i in unique(interim_dataset$drugclass)) {
  
  interim_dataset_drugclass <- interim_dataset %>%
    filter(drugclass == i)
  
  interim_number <- pROC::roc(response = interim_dataset_drugclass %>% select(stopdrug_6m_6mFU) %>% unlist(),
                              predictor = interim_dataset_drugclass %>% select(pred) %>% unlist(), ci = TRUE)$ci
  
  print(paste(i, signif(interim_number, digits = 3)))
  
}



# roc_coords_6m_by_drugs <- pROC::ci(pROC::roc(unlist(interim_dataset %>% filter(drugclass == "SGLT2") %>%select(stopdrug_6m_6mFU)), unlist(interim_dataset %>% filter(drugclass == "SGLT2") %>% select(pred))), of = "thresholds") %>%
#   as.data.frame() %>%
#   mutate(Therapy = "SGLT2i") %>%
#   rbind(
#     pROC::ci(pROC::roc(unlist(interim_dataset %>% filter(drugclass == "DPP4") %>%select(stopdrug_6m_6mFU)), unlist(interim_dataset %>% filter(drugclass == "DPP4") %>% select(pred))), of = "thresholds") %>%
#       as.data.frame() %>%
#       mutate(Therapy = "DPP4i"),
#     pROC::ci(pROC::roc(unlist(interim_dataset %>% filter(drugclass == "GLP1") %>%select(stopdrug_6m_6mFU)), unlist(interim_dataset %>% filter(drugclass == "GLP1") %>% select(pred))), of = "thresholds") %>%
#       as.data.frame() %>%
#       mutate(Therapy = "GLP1"),
#     pROC::ci(pROC::roc(unlist(interim_dataset %>% filter(drugclass == "TZD") %>%select(stopdrug_6m_6mFU)), unlist(interim_dataset %>% filter(drugclass == "TZD") %>% select(pred))), of = "thresholds") %>%
#       as.data.frame() %>%
#       mutate(Therapy = "TZD"),
#     pROC::ci(pROC::roc(unlist(interim_dataset %>% filter(drugclass == "SU") %>%select(stopdrug_6m_6mFU)), unlist(interim_dataset %>% filter(drugclass == "SU") %>% select(pred))), of = "thresholds") %>%
#       as.data.frame() %>%
#       mutate(Therapy = "SU")
#   )
# 
# saveRDS(roc_coords_6m_by_drugs, "results/Models/bartmachine/roc_coords_6m_by_drugs.rds")

roc_coords_6m_by_drugs <- readRDS("results/Models/bartmachine/roc_coords_6m_by_drugs.rds")


ppv_coords_6m <- pROC::coords(pROC::roc(
  cprd_dataset.6m %>%
    select(stopdrug_6m_6mFU) %>%
    mutate(stopdrug_6m_6mFU = factor(stopdrug_6m_6mFU)) %>%
    unlist(),
  bartmachine_full_model_6m$p_hat
), ret=c("threshold", "specificity", "ppv", "precision", "recall"))

plot_test <- ppv_coords_6m %>%
  as.data.frame() %>%
  ggplot() +
  geom_path(aes(x = threshold, y = ppv))


pdf("results/figures/06.precision_recall_6m.pdf", width = 5, height = 5.5)
ppv_coords_6m %>%
  as.data.frame() %>%
  ggplot() +
  geom_path(aes(x = recall, y = precision, colour = threshold)) +
  scale_colour_continuous(type = "viridis", limits = c(0, 1)) +
  theme_bw() +
  labs(x = "Recall", y = "Precision", colour = "Threshold") +
  theme(
    legend.position = "bottom",
    axis.text = element_text(size = 15),
    axis.title = element_text(size = 17),
    legend.title = element_text(size = 13)
  )
dev.off()



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


plot_roc_6m <- roc_coords_6m %>%
  as.data.frame() %>%
  ggplot() +
  geom_abline(aes(intercept = 0, slope = 1), colour = "black", linetype = "dashed") +
  geom_path(aes(x = 1 - specificity.2.5., y = sensitivity.2.5.), alpha = 0.4) +
  geom_path(aes(x = 1 - specificity.50., y = sensitivity.50.)) +
  geom_path(aes(x = 1 - specificity.97.5., y = sensitivity.97.5.), alpha = 0.4) +
  annotate("label", x = 0.60, y = 0.125, size = 5, label = paste0("AUROC=", signif(roc_values_6m_overall[2], digits = 3), ")")) + 
  scale_x_continuous("1 - Specificity", limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
  scale_y_continuous("Sensitivity", limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
  theme_bw() +
  theme(
    legend.position = "bottom",
    axis.title = element_text(size = 18),
    axis.text = element_text(size = 16)
  )


plot_roc_6m.drugs <- roc_coords_6m_by_drugs %>%
  as.data.frame() %>%
  ggplot() +
  geom_abline(aes(intercept = 0, slope = 1), colour = "black", linetype = "dashed") +
  # geom_path(aes(x = 1 - specificity.2.5., y = sensitivity.2.5., colour = Therapy), alpha = 0.4) +
  geom_path(aes(x = 1 - specificity.50., y = sensitivity.50., colour = Therapy)) +
  # geom_path(aes(x = 1 - specificity.97.5., y = sensitivity.97.5., colour = Therapy), alpha = 0.4) +
  scale_x_continuous("1 - Specificity", limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
  scale_y_continuous("Sensitivity", limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
  scale_colour_manual(values = c("Pooled" = "black", "SGLT2i" = "#E69F00", "GLP1" = "#56B4E9", "SU" = "#CC79A7", "DPP4i" = "#0072B2", "TZD" = "#D55E00", "MFN" = "grey"), breaks = rev(c("Pooled", "GLP1", "DPP4i", "SGLT2i", "TZD", "SU")), labels = rev(c("Pooled", "GLP-1RA", "DPP4i", "SGLT2i", "TZD", "SU")), name = "Therapy", guide = guide_legend(reverse = TRUE, nrow = 1)) +
  theme_bw() +
  theme(
    legend.position = "bottom",
    axis.title = element_text(size = 18),
    axis.text = element_text(size = 16)
  )


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


# Variable importance (sudo)
variable_importance_cohort_12m <- cprd_dataset.12m %>%
  mutate(
    pred = bartmachine_full_model_12m$p_hat
  )

m1 <- rms::ols(pred ~ rcs(dstartdate_age, 3) + gender + imd2015_10 + rcs(prebmi, 3) + rcs(dstartdate_dm_dur, 3) + rcs(prehba1c, 3) + predrug_frailty_proxy + ethnicity_5cat + predrug_bloodmed + smoking_cat + predrug_statins + rcs(preegfr, 3) + stopdrug_3m_3mFU_MFN_hist, data = variable_importance_cohort_12m, x = TRUE, y = TRUE)

values <- plot(anova(m1), what = 'proportion R2')

plot_variable_important_12m <- as.data.frame(values) %>%
  rownames_to_column() %>%
  left_join(
    data.frame(
      rowname = c("dstartdate_age", "gender", "imd2015_10", "prebmi", "dstartdate_dm_dur", "prehba1c", "drugline", "predrug_frailty_proxy", "ethnicity_5cat", "numdrugs", "predrug_bloodmed", "smoking_cat", "predrug_statins", "preegfr", "stopdrug_3m_3mFU_MFN_hist"),
      plotname = c("Age group", "Sex", "Index of multiple deprivation", "BMI", "Duration of diabetes", "HbA1c", "Number of therapy classes ever prescribed", "Frailty", "Ethnicity", "Number of other current glucose-lowering therapies", " Blood pressure medication", "Smoking", "Statins", "eGFR", "History of Metformin discontinuation within 3-months")
    ), by = c("rowname")
  ) %>%
  select(-rowname) %>%
  mutate(plotname = factor(plotname),
         values = values * 100) %>%
  ggplot(aes(y = forcats::fct_reorder(plotname, values), x = values)) +
  geom_segment(aes(x = 0, xend = values, yend = forcats::fct_reorder(plotname, values)), linetype = "dashed") +
  geom_point(size = 2, colour = "black") +
  ggtitle("Relative importance for discontinuation rate") +
  xlab("Relative Importance (%)") +
  theme_bw() +
  theme(axis.text.y = element_text(angle = 45, face = "bold"),
        axis.title.y = element_blank())



pdf("results/figures/06.var_importance_12m.pdf", width = 8, height = 6)
plot_variable_important_12m +
  theme(
    axis.text = element_text(size = 14),
    plot.title = element_text(size = 15, hjust = 0),
    plot.margin = unit(c(0.2, 0.2, 3.5, 0.2), "cm")
  )
dev.off()




roc_values_12m_overall <- pROC::roc(response = cprd_dataset.12m %>%
                                     select(stopdrug_12m_6mFU) %>%
                                     mutate(stopdrug_12m_6mFU = factor(stopdrug_12m_6mFU)) %>%
                                     unlist(),
                                   predictor = bartmachine_full_model_12m$p_hat, ci = TRUE)$ci %>%
  as.data.frame() %>%
  unlist()


plot_pooled_12m.overall <- calibration_plot(
  data = cprd_dataset.12m %>% mutate(pred = bartmachine_full_model_12m$p_hat) %>% mutate(stopdrug_12m_6mFU = as.numeric(stopdrug_12m_6mFU)-1) %>% mutate(grouping = "Pooled"), 
  obs = "stopdrug_12m_6mFU", 
  pred = "pred",
  group = "grouping",
  nTiles = 10)$calibration_plot +
  scale_x_continuous("Predicted discontinuation (%)", limits = c(0, 0.7), breaks = seq(0, 0.7, 0.1)) +
  scale_y_continuous("Observed discontinuation (%)", limits = c(0, 0.7), breaks = seq(0, 0.7, 0.1)) +
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
  scale_x_continuous("Predicted discontinuation (%)", limits = c(0, 0.7), breaks = seq(0, 0.7, 0.1)) +
  scale_y_continuous("Observed discontinuation (%)", limits = c(0, 0.7), breaks = seq(0, 0.7, 0.1)) +
  theme_bw() +
  scale_colour_manual(values = c("Pooled" = "black", "SGLT2" = "#E69F00", "GLP1" = "#56B4E9", "SU" = "#CC79A7", "DPP4" = "#0072B2", "TZD" = "#D55E00", "MFN" = "grey"), breaks = rev(c("Pooled", "GLP1", "DPP4", "SGLT2", "TZD", "SU")), labels = rev(c("Pooled", "GLP-1RA", "DPP4i", "SGLT2i", "TZD", "SU")), name = "Therapy", guide = guide_legend(reverse = TRUE, nrow = 1)) +
  theme(
    legend.position = "bottom"
  )


pdf("results/figures/06.calibration_overall_and_per_drug.pdf", width = 10, height = 13)

patchwork::wrap_plots(
  
  plot_pooled_3m.overall +
    scale_x_continuous("Predicted", limits = c(0, 0.7), breaks = seq(0, 1, 0.1), labels = scales::percent) +
    scale_y_continuous("Observed", limits = c(0, 0.7), breaks = seq(0, 1, 0.1), labels = scales::percent) +
    theme(legend.position = "none"),
  
  plot_pooled_3m.drugs +
    scale_x_continuous("Predicted", limits = c(0, 0.7), breaks = seq(0, 1, 0.1), labels = scales::percent) +
    scale_y_continuous("Observed", limits = c(0, 0.7), breaks = seq(0, 1, 0.1), labels = scales::percent) +
    theme(legend.position = "none"),
  
  plot_pooled_6m.overall +
    scale_x_continuous("Predicted", limits = c(0, 0.7), breaks = seq(0, 1, 0.1), labels = scales::percent) +
    scale_y_continuous("Observed", limits = c(0, 0.7), breaks = seq(0, 1, 0.1), labels = scales::percent) +
    theme(legend.position = "none"),
  
  plot_pooled_6m.drugs +
    scale_x_continuous("Predicted", limits = c(0, 0.7), breaks = seq(0, 1, 0.1), labels = scales::percent) +
    scale_y_continuous("Observed", limits = c(0, 0.7), breaks = seq(0, 1, 0.1), labels = scales::percent) +
    theme(legend.position = "none"),
  
  plot_pooled_12m.overall +
    scale_x_continuous("Predicted", limits = c(0, 0.7), breaks = seq(0, 1, 0.1), labels = scales::percent) +
    scale_y_continuous("Observed", limits = c(0, 0.7), breaks = seq(0, 1, 0.1), labels = scales::percent) +
    theme(legend.position = "none"),
  
  plot_pooled_12m.drugs +
    scale_x_continuous("Predicted", limits = c(0, 0.7), breaks = seq(0, 1, 0.1), labels = scales::percent) +
    scale_y_continuous("Observed", limits = c(0, 0.7), breaks = seq(0, 1, 0.1), labels = scales::percent) +
    theme(legend.position = "none"),
  
  ggpubr::as_ggplot(ggpubr::get_legend(plot_pooled_12m.drugs +
                                         theme(
                                           legend.text = element_text(size = 14),
                                           legend.title = element_text(size = 16)
                                         )))
  
) +
  patchwork::plot_layout(
    design = "
    AB
    CD
    EF
    GG
    ",
    heights = c(5,5,5,1)
  ) +
  patchwork::plot_annotation(tag_levels = list(c("A.1", "A.2", "B.1", "B.2", "C.1", "C.2", ""))) &
  theme(
    # legend.position = "bottom",
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 14)
  )

dev.off()



# pROC::roc(response = cprd_dataset.12m %>%
#             select(stopdrug_12m_6mFU) %>%
#             mutate(stopdrug_12m_6mFU = factor(stopdrug_12m_6mFU)) %>%
#             unlist(),
#           predictor = bartmachine_full_model_12m$p_hat, ci = TRUE)





ppv_coords_12m <- pROC::coords(pROC::roc(
  cprd_dataset.12m %>%
    select(stopdrug_12m_6mFU) %>%
    mutate(stopdrug_12m_6mFU = factor(stopdrug_12m_6mFU)) %>%
    unlist(),
  bartmachine_full_model_12m$p_hat
), ret=c("threshold", "specificity", "ppv", "precision", "recall"))

plot_test <- ppv_coords_12m %>%
  as.data.frame() %>%
  ggplot() +
  geom_path(aes(x = threshold, y = ppv))


pdf("results/figures/06.precision_recall_12m.pdf", width = 5, height = 5.5)
ppv_coords_12m %>%
  as.data.frame() %>%
  ggplot() +
  geom_path(aes(x = recall, y = precision, colour = threshold)) +
  scale_colour_continuous(type = "viridis", limits = c(0, 1)) +
  theme_bw() +
  labs(x = "Recall", y = "Precision", colour = "Threshold") +
  theme(
    legend.position = "bottom",
    axis.text = element_text(size = 15),
    axis.title = element_text(size = 17),
    legend.title = element_text(size = 13)
  )
dev.off()



interim_dataset <- cprd_dataset.12m %>%
  mutate(pred = bartmachine_full_model_12m$p_hat) %>%
  select(drugclass, stopdrug_12m_6mFU, pred)


for (i in unique(interim_dataset$drugclass)) {
  
  interim_dataset_drugclass <- interim_dataset %>%
    filter(drugclass == i)
  
  interim_number <- pROC::roc(response = interim_dataset_drugclass %>% select(stopdrug_12m_6mFU) %>% unlist(),
                              predictor = interim_dataset_drugclass %>% select(pred) %>% unlist(), ci = TRUE)$ci
  
  print(paste(i, signif(interim_number, digits = 3)))
  
}



# roc_coords_12m_by_drugs <- pROC::ci(pROC::roc(unlist(interim_dataset %>% filter(drugclass == "SGLT2") %>%select(stopdrug_12m_6mFU)), unlist(interim_dataset %>% filter(drugclass == "SGLT2") %>% select(pred))), of = "thresholds") %>%
#   as.data.frame() %>%
#   mutate(Therapy = "SGLT2i") %>%
#   rbind(
#     pROC::ci(pROC::roc(unlist(interim_dataset %>% filter(drugclass == "DPP4") %>%select(stopdrug_12m_6mFU)), unlist(interim_dataset %>% filter(drugclass == "DPP4") %>% select(pred))), of = "thresholds") %>%
#       as.data.frame() %>%
#       mutate(Therapy = "DPP4i"),
#     pROC::ci(pROC::roc(unlist(interim_dataset %>% filter(drugclass == "GLP1") %>%select(stopdrug_12m_6mFU)), unlist(interim_dataset %>% filter(drugclass == "GLP1") %>% select(pred))), of = "thresholds") %>%
#       as.data.frame() %>%
#       mutate(Therapy = "GLP1"),
#     pROC::ci(pROC::roc(unlist(interim_dataset %>% filter(drugclass == "TZD") %>%select(stopdrug_12m_6mFU)), unlist(interim_dataset %>% filter(drugclass == "TZD") %>% select(pred))), of = "thresholds") %>%
#       as.data.frame() %>%
#       mutate(Therapy = "TZD"),
#     pROC::ci(pROC::roc(unlist(interim_dataset %>% filter(drugclass == "SU") %>%select(stopdrug_12m_6mFU)), unlist(interim_dataset %>% filter(drugclass == "SU") %>% select(pred))), of = "thresholds") %>%
#       as.data.frame() %>%
#       mutate(Therapy = "SU")
#   )
# 
# saveRDS(roc_coords_12m_by_drugs, "results/Models/bartmachine/roc_coords_12m_by_drugs.rds")

roc_coords_12m_by_drugs <- readRDS("results/Models/bartmachine/roc_coords_12m_by_drugs.rds")



# roc_coords_12m <- pROC::ci(pROC::roc(
#   cprd_dataset.12m %>%
#     select(stopdrug_12m_6mFU) %>%
#     mutate(stopdrug_12m_6mFU = factor(stopdrug_12m_6mFU)) %>%
#     unlist(),
#   bartmachine_full_model_12m$p_hat
# ), of = "thresholds")
# 
# saveRDS(roc_coords_12m, "results/Models/bartmachine/roc_coords_12m.rds")

roc_coords_12m <- readRDS("results/Models/bartmachine/roc_coords_12m.rds")



plot_roc_12m <- roc_coords_12m %>%
  as.data.frame() %>%
  ggplot() +
  geom_abline(aes(intercept = 0, slope = 1), colour = "black", linetype = "dashed") +
  geom_path(aes(x = 1 - specificity.2.5., y = sensitivity.2.5.), alpha = 0.4) +
  geom_path(aes(x = 1 - specificity.50., y = sensitivity.50.)) +
  geom_path(aes(x = 1 - specificity.97.5., y = sensitivity.97.5.), alpha = 0.4) +
  annotate("label", x = 0.60, y = 0.125, size = 5, label = paste0("AUROC=", signif(roc_values_12m_overall[2], digits = 3), ")")) + 
  scale_x_continuous("1 - Specificity", limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
  scale_y_continuous("Sensitivity", limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
  theme_bw() +
  theme(
    legend.position = "bottom",
    axis.title = element_text(size = 18),
    axis.text = element_text(size = 16)
  )


plot_roc_12m.drugs <- roc_coords_12m_by_drugs %>%
  as.data.frame() %>%
  ggplot() +
  geom_abline(aes(intercept = 0, slope = 1), colour = "black", linetype = "dashed") +
  # geom_path(aes(x = 1 - specificity.2.5., y = sensitivity.2.5., colour = Therapy), alpha = 0.4) +
  geom_path(aes(x = 1 - specificity.50., y = sensitivity.50., colour = Therapy)) +
  # geom_path(aes(x = 1 - specificity.97.5., y = sensitivity.97.5., colour = Therapy), alpha = 0.4) +
  scale_x_continuous("1 - Specificity", limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
  scale_y_continuous("Sensitivity", limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
  scale_colour_manual(values = c("Pooled" = "black", "SGLT2i" = "#E69F00", "GLP1" = "#56B4E9", "SU" = "#CC79A7", "DPP4i" = "#0072B2", "TZD" = "#D55E00", "MFN" = "grey"), breaks = rev(c("Pooled", "GLP1", "DPP4i", "SGLT2i", "TZD", "SU")), labels = rev(c("Pooled", "GLP-1RA", "DPP4i", "SGLT2i", "TZD", "SU")), name = "Therapy", guide = guide_legend(reverse = TRUE, nrow = 1)) +
  theme_bw() +
  theme(
    legend.position = "bottom",
    axis.title = element_text(size = 18),
    axis.text = element_text(size = 16)
  )


pdf("results/figures/06.roc_curve_overall_and_per_drug.pdf", width = 10, height = 13)

patchwork::wrap_plots(
  
  plot_roc_3m +
    theme(legend.position = "none"),
  
  plot_roc_3m.drugs +
    theme(legend.position = "none"),
  
  plot_roc_6m +
    theme(legend.position = "none"),
  
  plot_roc_6m.drugs +
    theme(legend.position = "none"),
  
  plot_roc_12m +
    theme(legend.position = "none"),
  
  plot_roc_12m.drugs +
    theme(legend.position = "none"),
  
  ggpubr::as_ggplot(ggpubr::get_legend(plot_roc_12m.drugs +
                                         theme(
                                           legend.text = element_text(size = 14),
                                           legend.title = element_text(size = 16)
                                         )))
  
) +
  patchwork::plot_layout(
    design = "
    AB
    CD
    EF
    GG
    ",
    heights = c(5,5,5,1)
  ) +
  patchwork::plot_annotation(tag_levels = list(c("A.1", "A.2", "B.1", "B.2", "C.1", "C.2", ""))) &
  theme(
    # legend.position = "bottom",
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 14)
  )

dev.off()













###############################################################################
###############################################################################
######################## Fit the model: dev and val ###########################
###############################################################################
###############################################################################


# # Run model
# if (class(try(
#   
#   bartmachine_full_model <- readRDS("results/Models/bartmachine/bartmachine_full_model.rds")
#   
#   , silent = TRUE)) == "try-error") {
#   
#   bartmachine_full_model <- bartMachine::build_bart_machine(X = cprd_dataset.dev %>%
#                                                 select(drugclass, 
#                                                        dstartdate_age, 
#                                                        gender, 
#                                                        imd2015_10, 
#                                                        prebmi, 
#                                                        dstartdate_dm_dur, 
#                                                        prehba1c, 
#                                                        drugline, 
#                                                        predrug_frailty_proxy,
#                                                        ethnicity_5cat,
#                                                        numdrugs,
#                                                        predrug_bloodmed,
#                                                        smoking_cat,
#                                                        predrug_statins,
#                                                        preegfr,
#                                                        prehdl,
#                                                        stopdrug_3m_3mFU_MFN_hist),
#                                               y = cprd_dataset.dev %>% 
#                                                 select(stopdrug_3m_6mFU) %>%
#                                                 mutate(stopdrug_3m_6mFU = factor(stopdrug_3m_6mFU)) %>%
#                                                 unlist(),
#                                               num_trees = 200,
#                                               use_missing_data = FALSE,
#                                               num_burn_in = 15000,
#                                               num_iterations_after_burn_in = 10000,
#                                               serialize = TRUE)
#   
#   saveRDS(bartmachine_full_model, "results/Models/bartmachine/bartmachine_full_model.rds")
#   
# }
# 
# 
# # Check convergence
# # plot_convergence_diagnostics(bartmachine_full_model, plots =  c("sigsqs", "mh_acceptance", "num_nodes", "tree_depths"))
# 
# 
# 
# # Predict Development dataset
# if (class(try(
#   
#   bartmachine_pred_dev <- readRDS("results/Models/bartmachine/bartmachine_pred_dev.rds")
#   
#   , silent = TRUE)) == "try-error") {
#   
#   bartmachine_pred_dev <- predict(bartmachine_full_model,
#                                   cprd_dataset.dev %>%
#                                                               select(drugclass, 
#                                                                      dstartdate_age, 
#                                                                      gender, 
#                                                                      imd2015_10, 
#                                                                      prebmi, 
#                                                                      dstartdate_dm_dur, 
#                                                                      prehba1c, 
#                                                                      drugline, 
#                                                                      predrug_frailty_proxy,
#                                                                      ethnicity_5cat,
#                                                                      numdrugs,
#                                                                      predrug_bloodmed,
#                                                                      smoking_cat,
#                                                                      predrug_statins,
#                                                                      preegfr,
#                                                                      prehdl,
#                                                                      stopdrug_3m_3mFU_MFN_hist))
#   
#   saveRDS(bartmachine_pred_dev, "results/Models/bartmachine/bartmachine_pred_dev.rds")
#   
# }
# 
# 
# 
# # Predict Validation dataset
# if (class(try(
#   
#   bartmachine_pred_val <- readRDS("results/Models/bartmachine/bartmachine_pred_val.rds")
#   
#   , silent = TRUE)) == "try-error") {
#   
#   bartmachine_pred_val <- predict(bartmachine_full_model,
#                                   cprd_dataset.val %>%
#                                     select(drugclass, 
#                                            dstartdate_age, 
#                                            gender, 
#                                            imd2015_10, 
#                                            prebmi, 
#                                            dstartdate_dm_dur, 
#                                            prehba1c, 
#                                            drugline, 
#                                            predrug_frailty_proxy,
#                                            ethnicity_5cat,
#                                            numdrugs,
#                                            predrug_bloodmed,
#                                            smoking_cat,
#                                            predrug_statins,
#                                            preegfr,
#                                            prehdl,
#                                            stopdrug_3m_3mFU_MFN_hist))
#   
#   saveRDS(bartmachine_pred_val, "results/Models/bartmachine/bartmachine_pred_val.rds")
#   
# }
# 
# 
#  
# # pROC::roc(response = cprd_dataset.dev %>%
# #             select(stopdrug_3m_6mFU) %>%
# #             mutate(stopdrug_3m_6mFU = factor(stopdrug_3m_6mFU)) %>%
# #             unlist(),
# #           predictor = bartmachine_pred_dev, ci = TRUE)
# 
# 
# 
# 
# ################################################################################
# ################################################################################
# 
# 
# ################################################################################
# ################################################################################
# 
# 
# ## Try BART
# # Run model
# if (class(try(
#   
#   bart_full_model <- readRDS("results/Models/BART/bart_full_model.rds")
#   
#   , silent = TRUE)) == "try-error") {
#   
#   bart_full_model <- wbart(x.train = model.matrix(~., cprd_dataset.dev %>%
#                                                     select(drugclass, 
#                                                            dstartdate_age, 
#                                                            gender, 
#                                                            imd2015_10, 
#                                                            prebmi, 
#                                                            dstartdate_dm_dur, 
#                                                            prehba1c, 
#                                                            drugline, 
#                                                            predrug_frailty_proxy,
#                                                            ethnicity_5cat,
#                                                            numdrugs,
#                                                            predrug_bloodmed,
#                                                            smoking_cat,
#                                                            predrug_statins,
#                                                            preegfr,
#                                                            prehdl,
#                                                            stopdrug_3m_3mFU_MFN_hist)) %>% 
#                              as.data.frame() %>% 
#                              select(-`(Intercept)`),
#                            y.train = cprd_dataset.dev %>%
#                              select(stopdrug_3m_6mFU) %>%
#                              mutate(stopdrug_3m_6mFU = as.numeric(stopdrug_3m_6mFU) - 1) %>% 
#                              unlist(),
#                            nskip = 15000,
#                            ndpost = 10000)
#   
#   saveRDS(bart_full_model, "results/Models/BART/bart_full_model.rds")
#   
# }
# 
# 
# if (class(try(
#   
#   bart_pred_dev <- readRDS("results/Models/BART/bart_pred_dev.rds")
#   
#   , silent = TRUE)) == "try-error") {
#   
#   bart_pred_dev = predict(bart_full_model, model.matrix(~., cprd_dataset.dev %>%
#                                                      select(drugclass, 
#                                                             dstartdate_age, 
#                                                             gender, 
#                                                             imd2015_10, 
#                                                             prebmi, 
#                                                             dstartdate_dm_dur, 
#                                                             prehba1c, 
#                                                             drugline, 
#                                                             predrug_frailty_proxy,
#                                                             ethnicity_5cat,
#                                                             numdrugs,
#                                                             predrug_bloodmed,
#                                                             smoking_cat,
#                                                             predrug_statins,
#                                                             preegfr,
#                                                             prehdl,
#                                                             stopdrug_3m_3mFU_MFN_hist)) %>% 
#                        as.data.frame() %>% 
#                        select(-`(Intercept)`))
#   
#   saveRDS(bart_pred_dev, "results/Models/BART/bart_pred_dev.rds")
#   
# }
# 
# 
# 
# if (class(try(
#   
#   bart_pred_val <- readRDS("results/Models/BART/bart_pred_val.rds")
#   
#   , silent = TRUE)) == "try-error") {
#   
#   bart_pred_val = predict(bart_full_model, model.matrix(~., cprd_dataset.val %>%
#                                                      select(drugclass, 
#                                                             dstartdate_age, 
#                                                             gender, 
#                                                             imd2015_10, 
#                                                             prebmi, 
#                                                             dstartdate_dm_dur, 
#                                                             prehba1c, 
#                                                             drugline, 
#                                                             predrug_frailty_proxy,
#                                                             ethnicity_5cat,
#                                                             numdrugs,
#                                                             predrug_bloodmed,
#                                                             smoking_cat,
#                                                             predrug_statins,
#                                                             preegfr,
#                                                             prehdl,
#                                                             stopdrug_3m_3mFU_MFN_hist)) %>% 
#                        as.data.frame() %>% 
#                        select(-`(Intercept)`))
#   
#   saveRDS(bart_pred_val, "results/Models/BART/bart_pred_val.rds")
#   
# }
# 
# 
# 
# 
# # pROC::roc(response = cprd_dataset.val %>%
# #             select(stopdrug_3m_6mFU) %>%
# #             mutate(stopdrug_3m_6mFU = factor(stopdrug_3m_6mFU)) %>%
# #             unlist(),
# #           predictor = colMeans(bart_pred_val), ci = TRUE)
# 
# 




