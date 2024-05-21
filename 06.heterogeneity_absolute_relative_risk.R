####################
## Description: 
##  - In this file we:
##    - Investigate heterogeneity of therapy discontinuation
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
cprd_dataset.dev <- set_up_data(
  raw_data = "20240308_t2d_1stinstance",
  diagnosis = FALSE,
  therapies = c("DPP4", "GLP1", "MFN", "SGLT2", "SU", "TZD"),
  dataset = "3m.disc.dataset.dev"
) %>%
  drop_na(-stopdrug_6m_6mFU, -stopdrug_12m_6mFU)


# load dataset
cprd_dataset.val <- set_up_data(
  raw_data = "20240308_t2d_1stinstance",
  diagnosis = FALSE,
  therapies = c("DPP4", "GLP1", "MFN", "SGLT2", "SU", "TZD"),
  dataset = "3m.disc.dataset.val"
) %>%
  drop_na(-stopdrug_6m_6mFU, -stopdrug_12m_6mFU)


probabilities.only_dataset <- readRDS("results/Models/Predictions/model_predictions_3m_all_drugs.rds")


# join predictions into dataset
cprd_dataset.dev <- cprd_dataset.dev %>%
  left_join(
    probabilities.only_dataset,
    by = c("patid", "dstartdate")
  )
cprd_dataset.val <- cprd_dataset.val %>%
  left_join(
    probabilities.only_dataset,
    by = c("patid", "dstartdate")
  )



###############################################################################
###############################################################################
####################### Per combination comparison ############################
###############################################################################
###############################################################################


breakpoints <- c(-0.1, -0.05, -0.025, 0, 0.025, 0.05, 0.1)
# Odds ratios 
odds_ratios.no_weight.all <- calc_odds_ratios(cprd_dataset.dev %>% rbind(cprd_dataset.val), drugs = c("MFN", "DPP4", "GLP1", "SGLT2", "SU", "TZD"), pred.variable = "no_weight", break_points = breakpoints)



pdf("results/figures/06.plot_1.pdf", width = 10, height = 7)

odds_ratios.no_weight.all %>%
  drop_na() %>%
  mutate(
    obs = exp(as.numeric(obs)),
    lci = exp(as.numeric(lci)),
    uci = exp(as.numeric(uci))
  ) %>%
  filter(drug_1 == "MFN" | drug_2 == "MFN") %>%
  ggplot() + 
  geom_vline(aes(xintercept = 1), colour = "grey") +
  geom_point(aes(y = group, x = obs)) + 
  geom_errorbar(aes(y = group, xmin = lci, xmax = uci, x = obs)) + 
  xlab("Odds Ratio") +
  scale_x_log10(breaks = c(0.25, 0.5, 1, 2, 4)) +
  coord_cartesian(xlim = c(0.2, 4.1)) +
  ggtitle("Development + Validation") +
  facet_wrap(~Type, scales = "free_x") +
  theme_bw() +
  theme(
    axis.title.y = element_blank()
  )

dev.off()



###############################################################################
###############################################################################
####################### Per combination comparison ############################
###############################################################################
###############################################################################



# load dataset
cprd_dataset.dev <- set_up_data(
  raw_data = "20240308_t2d_1stinstance",
  diagnosis = FALSE,
  therapies = c("DPP4", "GLP1", "MFN", "SGLT2", "SU", "TZD"),
  dataset = "3m.disc.dataset.dev"
) %>%
  drop_na(-stopdrug_6m_6mFU, -stopdrug_12m_6mFU)


# load dataset
cprd_dataset.val <- set_up_data(
  raw_data = "20240308_t2d_1stinstance",
  diagnosis = FALSE,
  therapies = c("DPP4", "GLP1", "MFN", "SGLT2", "SU", "TZD"),
  dataset = "3m.disc.dataset.val"
) %>%
  drop_na(-stopdrug_6m_6mFU, -stopdrug_12m_6mFU)


probabilities.only_dataset <- readRDS("results/Models/Predictions/model_predictions_3m_2nd_line_drugs.rds")


# join predictions into dataset
cprd_dataset.dev <- cprd_dataset.dev %>%
  left_join(
    probabilities.only_dataset,
    by = c("patid", "dstartdate")
  )
cprd_dataset.val <- cprd_dataset.val %>%
  left_join(
    probabilities.only_dataset,
    by = c("patid", "dstartdate")
  )



breakpoints <- c(-0.1, -0.05, -0.025, 0, 0.025, 0.05, 0.1)
# Odds ratios 
odds_ratios.no_weight.2nd_line <- calc_odds_ratios(cprd_dataset.dev %>% rbind(cprd_dataset.val), drugs = c("DPP4", "GLP1", "SGLT2", "SU", "TZD"), pred.variable = "no_weight", break_points = breakpoints)



pdf("results/figures/06.plot_2.pdf", width = 12, height = 9)

odds_ratios.no_weight.2nd_line %>%
  drop_na() %>%
  mutate(
    obs = exp(as.numeric(obs)),
    lci = exp(as.numeric(lci)),
    uci = exp(as.numeric(uci))
  ) %>%
  filter(drug_1 != "MFN" & drug_2 != "MFN") %>%
  ggplot() + 
  geom_vline(aes(xintercept = 1), colour = "grey") +
  geom_point(aes(y = group, x = obs)) + 
  geom_errorbar(aes(y = group, xmin = lci, xmax = uci, x = obs)) + 
  xlab("Odds Ratio") +
  scale_x_log10(breaks = c(0.25, 0.5, 1, 2, 4)) +
  coord_cartesian(xlim = c(0.2, 4.1)) +
  ggtitle("Development + Validation") +
  facet_wrap(~Type, scales = "free_x") +
  theme_bw() +
  theme(
    axis.title.y = element_blank()
  )

dev.off()


#:----------------------------
# Combining heterogeneity and odds ratios

breakpoints <- c(-0.1, -0.05, -0.025, 0, 0.025, 0.05, 0.1)

# ATE not adjusted
ATE.var_adj.no_adj_2nd_line <- calc_ATE(cprd_dataset.dev, break_points = breakpoints, drugs = c("GLP1", "SGLT2"), pred.variable = "no_weight", n_bootstrap = 10)


# Odds ratios 
odds_ratios.no_weight.2nd_line <- calc_odds_ratios(cprd_dataset.dev %>% rbind(cprd_dataset.val), drugs = c("GLP1", "SGLT2"), pred.variable = "no_weight", break_points = breakpoints)

title <- ATE.var_adj.no_adj_2nd_line[1,]
title[,] <- NA

interim_1 <- title %>% as.data.frame() %>%mutate(group = "Predicted benefit on GLP1:") %>%
  rbind(
    ATE.var_adj.no_adj_2nd_line %>%
      slice(which(grepl("(-", ATE.var_adj.no_adj_2nd_line$group, fixed =TRUE) == TRUE)),
    title %>% as.data.frame() %>% mutate(group = "Predicted benefit on SGLT2:"),
    ATE.var_adj.no_adj_2nd_line %>%
      slice(which(!grepl("(-", ATE.var_adj.no_adj_2nd_line$group, fixed =TRUE) == TRUE))
    
  )

title <- odds_ratios.no_weight.2nd_line[1,]
title[,] <- NA


interim_2 <- title %>% as.data.frame() %>%mutate(group = "Predicted benefit on GLP1:") %>%
  rbind(
    odds_ratios.no_weight.2nd_line %>%
      slice(which(grepl("(-", odds_ratios.no_weight.2nd_line$group, fixed =TRUE) == TRUE)),
    title %>% as.data.frame() %>% mutate(group = "Predicted benefit on SGLT2:"),
    odds_ratios.no_weight.2nd_line %>%
      slice(which(!grepl("(-", odds_ratios.no_weight.2nd_line$group, fixed =TRUE) == TRUE))
    
  )

patchwork::wrap_plots(
  
  interim_1 %>%
    ggplot(aes(y = group, x = obs, xmin = lci, xmax = uci)) +
    geom_vline(aes(xintercept = 0), colour = "grey") +
    geom_point() +
    geom_errorbar() + 
    scale_y_discrete(
      limits= rev(interim_1$group),
      labels = rev(c(
        "Predicted benefit on GLP1:",
        "Benefit >10% (n=112, events = 18)",
        "Benefit 5-10% (n=3,216, events = 644)",
        "Benefit 2.5-5% (n=9,122, events = 1,522)",
        "Benefit 0-2.5% (n=19,848, events = 2,749)",
        "Predicted benefit on SGLT2:",
        "Benefit 0-2.5% (n=19,518, events = 2,382)",
        "Benefit 2.5-5% (n=4,264, events = 558)",
        "Benefit 5-10% (n=407, events = 68)",
        "Benefit >10% (n=4, events = 0)"
      ))
    ) +
    scale_x_continuous("Predicted Conditional Average Treatment Discontinuation Difference", labels = scales::percent, breaks = seq(-1, 1, 0.02)) +
    coord_cartesian(xlim = c(-0.12, 0.06)) +
    ggtitle("Absolute benefit (below 0 = benefit on GLP1)") +
    theme_bw(),
  
  
  
  interim_2 %>%
    mutate(
      obs = exp(as.numeric(obs)),
      lci = exp(as.numeric(lci)),
      uci = exp(as.numeric(uci))
    ) %>%
    ggplot() + 
    geom_vline(aes(xintercept = 1), colour = "grey") +
    geom_point(aes(y = group, x = obs)) + 
    geom_errorbar(aes(y = group, xmin = lci, xmax = uci, x = obs)) + 
    scale_y_discrete(
      limits= rev(interim_2$group),
      labels = rev(c(
        "Predicted benefit on GLP1:",
        "Benefit >10% (n=112, events = 18)",
        "Benefit 5-10% (n=3,216, events = 644)",
        "Benefit 2.5-5% (n=9,122, events = 1,522)",
        "Benefit 0-2.5% (n=19,848, events = 2,749)",
        "Predicted benefit on SGLT2:",
        "Benefit 0-2.5% (n=19,518, events = 2,382)",
        "Benefit 2.5-5% (n=4,264, events = 558)",
        "Benefit 5-10% (n=407, events = 68)",
        "Benefit >10% (n=4, events = 0)"
      ))
    ) +
    xlab("Odds Ratio") +
    scale_x_log10(breaks = c(0.25, 0.5, 1, 2, 4)) +
    coord_cartesian(xlim = c(0.2, 4.1)) +
    ggtitle("Relative benefit (below 1 = benefit on GLP1)") +
    theme_bw()
  
) + plot_layout(axes = "collect") +
  plot_annotation(
    title = 'GLP1 vs SGLT2'
  ) &
  theme(
    axis.title.y = element_blank(),
    axis.title.x = element_text(size = 13),
    axis.text = element_text(size = 15)
  )





#:-----------------------------------------------------------------
# Combining pdfs

qpdf::pdf_combine(input = c("results/figures/06.plot_1.pdf", "results/figures/06.plot_2.pdf"),
                  output = "results/figures/heterogeneity_odds_ratios.pdf")


file.remove(c("results/figures/06.plot_1.pdf", "results/figures/06.plot_2.pdf"))

