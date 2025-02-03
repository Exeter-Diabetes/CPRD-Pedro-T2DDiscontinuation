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
  dataset = "3m.disc.dataset.dev",
  follow_up = "6-months",
  full_prescribing_history = TRUE
) %>%
  drop_na(-stopdrug_6m_6mFU, -stopdrug_12m_6mFU)


# load dataset
cprd_dataset.val <- set_up_data(
  raw_data = "20240308_t2d_1stinstance",
  diagnosis = FALSE,
  therapies = c("DPP4", "GLP1", "MFN", "SGLT2", "SU", "TZD"),
  dataset = "3m.disc.dataset.val",
  follow_up = "6-months",
  full_prescribing_history = TRUE
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
odds_ratios.no_weight.all <- calc_odds_ratios(cprd_dataset.val, drugs = c("MFN", "DPP4", "GLP1", "SGLT2", "SU", "TZD"), pred.variable = ".no_weight", break_points = breakpoints)



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
  dataset = "3m.disc.dataset.dev",
  follow_up = "6-months",
  full_prescribing_history = TRUE
) %>%
  drop_na(-stopdrug_6m_6mFU, -stopdrug_12m_6mFU)


# load dataset
cprd_dataset.val <- set_up_data(
  raw_data = "20240308_t2d_1stinstance",
  diagnosis = FALSE,
  therapies = c("DPP4", "GLP1", "MFN", "SGLT2", "SU", "TZD"),
  dataset = "3m.disc.dataset.val",
  follow_up = "6-months",
  full_prescribing_history = TRUE
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
odds_ratios.no_weight.2nd_line <- calc_odds_ratios(cprd_dataset.val, drugs = c("DPP4", "GLP1", "SGLT2", "SU", "TZD"), pred.variable = ".no_weight", break_points = breakpoints)



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



#:-----------------------------------------------------------------
# Combining pdfs

qpdf::pdf_combine(input = c("results/figures/06.plot_1.pdf", "results/figures/06.plot_2.pdf"),
                  output = "results/figures/heterogeneity_odds_ratios.pdf")


file.remove(c("results/figures/06.plot_1.pdf", "results/figures/06.plot_2.pdf"))


#:----------------------------
#:----------------------------
#:----------------------------


# Combining heterogeneity and odds ratios

breakpoints <- c(-0.1, -0.05, -0.025, 0, 0.025, 0.05, 0.1)

# ATE not adjusted
ATE.var_adj.no_adj_2nd_line <- calc_ATE(cprd_dataset.dev, break_points = breakpoints, drugs = c("DPP4", "GLP1", "SGLT2", "SU", "TZD"), pred.variable = ".no_weight", n_bootstrap = 10) %>%
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
absolute.var_adj.no_adj_2n_line <- calc_predicted_risk(cprd_dataset.dev, break_points = breakpoints, drugs = c("DPP4", "GLP1", "SGLT2", "SU", "TZD"), pred.variable = ".no_weight") %>%
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
odds_ratios.no_weight.2nd_line <- calc_odds_ratios(cprd_dataset.dev, drugs = c("DPP4", "GLP1", "SGLT2", "SU", "TZD"), pred.variable = ".no_weight", break_points = breakpoints) %>%
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



pdf("results/figures/absolute_relative_development.pdf", width = 20, height = 7)

##: Per drug combination: DPP4vsGLP1
drug1 = "DPP4"; drug2 = "GLP1"

title <- ATE.var_adj.no_adj_2nd_line[1,]
title[,] <- NA

interim_1 <- title %>% as.data.frame() %>%mutate(group = paste0("Predicted benefit on ", drug1, ":"), label = paste0("Predicted benefit on ", drug1, ":")) %>%
  rbind(
    ATE.var_adj.no_adj_2nd_line %>%
      filter(drug_1 == drug1 & drug_2 == drug2) %>%
      slice(which(grepl("(-", ATE.var_adj.no_adj_2nd_line$group, fixed =TRUE) == TRUE)),
    title %>% as.data.frame() %>% mutate(group = paste0("Predicted benefit on ", drug2, ":"), label = paste0("Predicted benefit on ", drug2, ":")),
    ATE.var_adj.no_adj_2nd_line %>%
      filter(drug_1 == drug1 & drug_2 == drug2) %>%
      slice(which(!grepl("(-", ATE.var_adj.no_adj_2nd_line$group, fixed =TRUE) == TRUE))
    
  )

title <- odds_ratios.no_weight.2nd_line[1,]
title[,] <- NA


interim_2 <- title %>% as.data.frame() %>%mutate(group = paste0("Predicted benefit on ", drug1, ":"), label = paste0("Predicted benefit on ", drug1, ":")) %>%
  rbind(
    odds_ratios.no_weight.2nd_line %>%
      filter(drug_1 == drug1 & drug_2 == drug2) %>%
      slice(which(grepl("(-", odds_ratios.no_weight.2nd_line$group, fixed =TRUE) == TRUE)),
    title %>% as.data.frame() %>% mutate(group = paste0("Predicted benefit on ", drug2, ":"), label = paste0("Predicted benefit on ", drug2, ":")),
    odds_ratios.no_weight.2nd_line %>%
      filter(drug_1 == drug1 & drug_2 == drug2) %>%
      slice(which(!grepl("(-", odds_ratios.no_weight.2nd_line$group, fixed =TRUE) == TRUE))
    
  )

title <- absolute.var_adj.no_adj_2n_line[1,]
title[,] <- NA

interim_3 <- title %>% as.data.frame() %>% mutate(group = paste0("Predicted benefit on ", drug1, ":"), label = paste0("Predicted benefit on ", drug1, ":")) %>%
  mutate(drug_type = drug1) %>%
  rbind(
    absolute.var_adj.no_adj_2n_line %>%
      filter(drug_1 == drug1 & drug_2 == drug2) %>%
      slice(which(grepl("(-", absolute.var_adj.no_adj_2n_line$group, fixed =TRUE) == TRUE)),
    title %>% as.data.frame() %>% mutate(group = paste0("Predicted benefit on ", drug2, ":"), label = paste0("Predicted benefit on ", drug2, ":")) %>%
      mutate(drug_type = drug2),
    absolute.var_adj.no_adj_2n_line %>%
      filter(drug_1 == drug1 & drug_2 == drug2) %>%
      slice(which(!grepl("(-", absolute.var_adj.no_adj_2n_line$group, fixed =TRUE) == TRUE))
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
  
) + plot_layout(axes = "collect") +
  plot_annotation(
    title = paste0(drug1, " vs ", drug2)
  ) &
  theme(
    axis.title.y = element_blank(),
    axis.title.x = element_text(size = 13),
    axis.text = element_text(size = 13)
  )


##: Per drug combination: DPP4vsSGLT2
drug1 = "DPP4"; drug2 = "SGLT2"

title <- ATE.var_adj.no_adj_2nd_line[1,]
title[,] <- NA

interim_1 <- title %>% as.data.frame() %>%mutate(group = paste0("Predicted benefit on ", drug1, ":"), label = paste0("Predicted benefit on ", drug1, ":")) %>%
  rbind(
    ATE.var_adj.no_adj_2nd_line %>%
      filter(drug_1 == drug1 & drug_2 == drug2) %>%
      slice(which(grepl("(-", ATE.var_adj.no_adj_2nd_line$group, fixed =TRUE) == TRUE)),
    title %>% as.data.frame() %>% mutate(group = paste0("Predicted benefit on ", drug2, ":"), label = paste0("Predicted benefit on ", drug2, ":")),
    ATE.var_adj.no_adj_2nd_line %>%
      filter(drug_1 == drug1 & drug_2 == drug2) %>%
      slice(which(!grepl("(-", ATE.var_adj.no_adj_2nd_line$group, fixed =TRUE) == TRUE))
    
  )

title <- odds_ratios.no_weight.2nd_line[1,]
title[,] <- NA


interim_2 <- title %>% as.data.frame() %>%mutate(group = paste0("Predicted benefit on ", drug1, ":"), label = paste0("Predicted benefit on ", drug1, ":")) %>%
  rbind(
    odds_ratios.no_weight.2nd_line %>%
      filter(drug_1 == drug1 & drug_2 == drug2) %>%
      slice(which(grepl("(-", odds_ratios.no_weight.2nd_line$group, fixed =TRUE) == TRUE)),
    title %>% as.data.frame() %>% mutate(group = paste0("Predicted benefit on ", drug2, ":"), label = paste0("Predicted benefit on ", drug2, ":")),
    odds_ratios.no_weight.2nd_line %>%
      filter(drug_1 == drug1 & drug_2 == drug2) %>%
      slice(which(!grepl("(-", odds_ratios.no_weight.2nd_line$group, fixed =TRUE) == TRUE))
    
  )

title <- absolute.var_adj.no_adj_2n_line[1,]
title[,] <- NA

interim_3 <- title %>% as.data.frame() %>% mutate(group = paste0("Predicted benefit on ", drug1, ":"), label = paste0("Predicted benefit on ", drug1, ":")) %>%
  mutate(drug_type = drug1) %>%
  rbind(
    absolute.var_adj.no_adj_2n_line %>%
      filter(drug_1 == drug1 & drug_2 == drug2) %>%
      slice(which(grepl("(-", absolute.var_adj.no_adj_2n_line$group, fixed =TRUE) == TRUE)),
    title %>% as.data.frame() %>% mutate(group = paste0("Predicted benefit on ", drug2, ":"), label = paste0("Predicted benefit on ", drug2, ":")) %>%
      mutate(drug_type = drug2),
    absolute.var_adj.no_adj_2n_line %>%
      filter(drug_1 == drug1 & drug_2 == drug2) %>%
      slice(which(!grepl("(-", absolute.var_adj.no_adj_2n_line$group, fixed =TRUE) == TRUE))
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
    coord_cartesian(xlim = c(-0.3, 0.1)) +
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
    coord_cartesian(xlim = c(0.05, 0.4)) +
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
  
) + plot_layout(axes = "collect") +
  plot_annotation(
    title = paste0(drug1, " vs ", drug2)
  ) &
  theme(
    axis.title.y = element_blank(),
    axis.title.x = element_text(size = 13),
    axis.text = element_text(size = 13)
  )


##: Per drug combination: DPP4vsSU
drug1 = "DPP4"; drug2 = "SU"

title <- ATE.var_adj.no_adj_2nd_line[1,]
title[,] <- NA

interim_1 <- title %>% as.data.frame() %>%mutate(group = paste0("Predicted benefit on ", drug1, ":"), label = paste0("Predicted benefit on ", drug1, ":")) %>%
  rbind(
    ATE.var_adj.no_adj_2nd_line %>%
      filter(drug_1 == drug1 & drug_2 == drug2) %>%
      slice(which(grepl("(-", ATE.var_adj.no_adj_2nd_line$group, fixed =TRUE) == TRUE)),
    title %>% as.data.frame() %>% mutate(group = paste0("Predicted benefit on ", drug2, ":"), label = paste0("Predicted benefit on ", drug2, ":")),
    ATE.var_adj.no_adj_2nd_line %>%
      filter(drug_1 == drug1 & drug_2 == drug2) %>%
      slice(which(!grepl("(-", ATE.var_adj.no_adj_2nd_line$group, fixed =TRUE) == TRUE))
    
  )

title <- odds_ratios.no_weight.2nd_line[1,]
title[,] <- NA


interim_2 <- title %>% as.data.frame() %>%mutate(group = paste0("Predicted benefit on ", drug1, ":"), label = paste0("Predicted benefit on ", drug1, ":")) %>%
  rbind(
    odds_ratios.no_weight.2nd_line %>%
      filter(drug_1 == drug1 & drug_2 == drug2) %>%
      slice(which(grepl("(-", odds_ratios.no_weight.2nd_line$group, fixed =TRUE) == TRUE)),
    title %>% as.data.frame() %>% mutate(group = paste0("Predicted benefit on ", drug2, ":"), label = paste0("Predicted benefit on ", drug2, ":")),
    odds_ratios.no_weight.2nd_line %>%
      filter(drug_1 == drug1 & drug_2 == drug2) %>%
      slice(which(!grepl("(-", odds_ratios.no_weight.2nd_line$group, fixed =TRUE) == TRUE))
    
  )

title <- absolute.var_adj.no_adj_2n_line[1,]
title[,] <- NA

interim_3 <- title %>% as.data.frame() %>% mutate(group = paste0("Predicted benefit on ", drug1, ":"), label = paste0("Predicted benefit on ", drug1, ":")) %>%
  mutate(drug_type = drug1) %>%
  rbind(
    absolute.var_adj.no_adj_2n_line %>%
      filter(drug_1 == drug1 & drug_2 == drug2) %>%
      slice(which(grepl("(-", absolute.var_adj.no_adj_2n_line$group, fixed =TRUE) == TRUE)),
    title %>% as.data.frame() %>% mutate(group = paste0("Predicted benefit on ", drug2, ":"), label = paste0("Predicted benefit on ", drug2, ":")) %>%
      mutate(drug_type = drug2),
    absolute.var_adj.no_adj_2n_line %>%
      filter(drug_1 == drug1 & drug_2 == drug2) %>%
      slice(which(!grepl("(-", absolute.var_adj.no_adj_2n_line$group, fixed =TRUE) == TRUE))
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
  
) + plot_layout(axes = "collect") +
  plot_annotation(
    title = paste0(drug1, " vs ", drug2)
  ) &
  theme(
    axis.title.y = element_blank(),
    axis.title.x = element_text(size = 13),
    axis.text = element_text(size = 13)
  )


##: Per drug combination: DPP4vsTZD
drug1 = "DPP4"; drug2 = "TZD"

title <- ATE.var_adj.no_adj_2nd_line[1,]
title[,] <- NA

interim_1 <- title %>% as.data.frame() %>%mutate(group = paste0("Predicted benefit on ", drug1, ":"), label = paste0("Predicted benefit on ", drug1, ":")) %>%
  rbind(
    ATE.var_adj.no_adj_2nd_line %>%
      filter(drug_1 == drug1 & drug_2 == drug2) %>%
      slice(which(grepl("(-", ATE.var_adj.no_adj_2nd_line$group, fixed =TRUE) == TRUE)),
    title %>% as.data.frame() %>% mutate(group = paste0("Predicted benefit on ", drug2, ":"), label = paste0("Predicted benefit on ", drug2, ":")),
    ATE.var_adj.no_adj_2nd_line %>%
      filter(drug_1 == drug1 & drug_2 == drug2) %>%
      slice(which(!grepl("(-", ATE.var_adj.no_adj_2nd_line$group, fixed =TRUE) == TRUE))
    
  )

title <- odds_ratios.no_weight.2nd_line[1,]
title[,] <- NA


interim_2 <- title %>% as.data.frame() %>%mutate(group = paste0("Predicted benefit on ", drug1, ":"), label = paste0("Predicted benefit on ", drug1, ":")) %>%
  rbind(
    odds_ratios.no_weight.2nd_line %>%
      filter(drug_1 == drug1 & drug_2 == drug2) %>%
      slice(which(grepl("(-", odds_ratios.no_weight.2nd_line$group, fixed =TRUE) == TRUE)),
    title %>% as.data.frame() %>% mutate(group = paste0("Predicted benefit on ", drug2, ":"), label = paste0("Predicted benefit on ", drug2, ":")),
    odds_ratios.no_weight.2nd_line %>%
      filter(drug_1 == drug1 & drug_2 == drug2) %>%
      slice(which(!grepl("(-", odds_ratios.no_weight.2nd_line$group, fixed =TRUE) == TRUE))
    
  )

title <- absolute.var_adj.no_adj_2n_line[1,]
title[,] <- NA

interim_3 <- title %>% as.data.frame() %>% mutate(group = paste0("Predicted benefit on ", drug1, ":"), label = paste0("Predicted benefit on ", drug1, ":")) %>%
  mutate(drug_type = drug1) %>%
  rbind(
    absolute.var_adj.no_adj_2n_line %>%
      filter(drug_1 == drug1 & drug_2 == drug2) %>%
      slice(which(grepl("(-", absolute.var_adj.no_adj_2n_line$group, fixed =TRUE) == TRUE)),
    title %>% as.data.frame() %>% mutate(group = paste0("Predicted benefit on ", drug2, ":"), label = paste0("Predicted benefit on ", drug2, ":")) %>%
      mutate(drug_type = drug2),
    absolute.var_adj.no_adj_2n_line %>%
      filter(drug_1 == drug1 & drug_2 == drug2) %>%
      slice(which(!grepl("(-", absolute.var_adj.no_adj_2n_line$group, fixed =TRUE) == TRUE))
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
  
) + plot_layout(axes = "collect") +
  plot_annotation(
    title = paste0(drug1, " vs ", drug2)
  ) &
  theme(
    axis.title.y = element_blank(),
    axis.title.x = element_text(size = 13),
    axis.text = element_text(size = 13)
  )


##: Per drug combination: GLP1vsSGLT2
drug1 = "GLP1"; drug2 = "SGLT2"

title <- ATE.var_adj.no_adj_2nd_line[1,]
title[,] <- NA

interim_1 <- title %>% as.data.frame() %>%mutate(group = paste0("Predicted benefit on ", drug1, ":"), label = paste0("Predicted benefit on ", drug1, ":")) %>%
  rbind(
    ATE.var_adj.no_adj_2nd_line %>%
      filter(drug_1 == drug1 & drug_2 == drug2) %>%
      slice(which(grepl("(-", ATE.var_adj.no_adj_2nd_line$group, fixed =TRUE) == TRUE)),
    title %>% as.data.frame() %>% mutate(group = paste0("Predicted benefit on ", drug2, ":"), label = paste0("Predicted benefit on ", drug2, ":")),
    ATE.var_adj.no_adj_2nd_line %>%
      filter(drug_1 == drug1 & drug_2 == drug2) %>%
      slice(which(!grepl("(-", ATE.var_adj.no_adj_2nd_line$group, fixed =TRUE) == TRUE))
    
  )

title <- odds_ratios.no_weight.2nd_line[1,]
title[,] <- NA


interim_2 <- title %>% as.data.frame() %>%mutate(group = paste0("Predicted benefit on ", drug1, ":"), label = paste0("Predicted benefit on ", drug1, ":")) %>%
  rbind(
    odds_ratios.no_weight.2nd_line %>%
      filter(drug_1 == drug1 & drug_2 == drug2) %>%
      slice(which(grepl("(-", odds_ratios.no_weight.2nd_line$group, fixed =TRUE) == TRUE)),
    title %>% as.data.frame() %>% mutate(group = paste0("Predicted benefit on ", drug2, ":"), label = paste0("Predicted benefit on ", drug2, ":")),
    odds_ratios.no_weight.2nd_line %>%
      filter(drug_1 == drug1 & drug_2 == drug2) %>%
      slice(which(!grepl("(-", odds_ratios.no_weight.2nd_line$group, fixed =TRUE) == TRUE))
    
  )

title <- absolute.var_adj.no_adj_2n_line[1,]
title[,] <- NA

interim_3 <- title %>% as.data.frame() %>% mutate(group = paste0("Predicted benefit on ", drug1, ":"), label = paste0("Predicted benefit on ", drug1, ":")) %>%
  mutate(drug_type = drug1) %>%
  rbind(
    absolute.var_adj.no_adj_2n_line %>%
      filter(drug_1 == drug1 & drug_2 == drug2) %>%
      slice(which(grepl("(-", absolute.var_adj.no_adj_2n_line$group, fixed =TRUE) == TRUE)),
    title %>% as.data.frame() %>% mutate(group = paste0("Predicted benefit on ", drug2, ":"), label = paste0("Predicted benefit on ", drug2, ":")) %>%
    mutate(drug_type = drug2),
    absolute.var_adj.no_adj_2n_line %>%
      filter(drug_1 == drug1 & drug_2 == drug2) %>%
      slice(which(!grepl("(-", absolute.var_adj.no_adj_2n_line$group, fixed =TRUE) == TRUE))
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
  
) + plot_layout(axes = "collect") +
  plot_annotation(
    title = paste0(drug1, " vs ", drug2)
  ) &
  theme(
    axis.title.y = element_blank(),
    axis.title.x = element_text(size = 13),
    axis.text = element_text(size = 13)
  )



##: Per drug combination: GLP1vsSU
drug1 = "GLP1"; drug2 = "SU"

title <- ATE.var_adj.no_adj_2nd_line[1,]
title[,] <- NA

interim_1 <- title %>% as.data.frame() %>%mutate(group = paste0("Predicted benefit on ", drug1, ":"), label = paste0("Predicted benefit on ", drug1, ":")) %>%
  rbind(
    ATE.var_adj.no_adj_2nd_line %>%
      filter(drug_1 == drug1 & drug_2 == drug2) %>%
      slice(which(grepl("(-", ATE.var_adj.no_adj_2nd_line$group, fixed =TRUE) == TRUE)),
    title %>% as.data.frame() %>% mutate(group = paste0("Predicted benefit on ", drug2, ":"), label = paste0("Predicted benefit on ", drug2, ":")),
    ATE.var_adj.no_adj_2nd_line %>%
      filter(drug_1 == drug1 & drug_2 == drug2) %>%
      slice(which(!grepl("(-", ATE.var_adj.no_adj_2nd_line$group, fixed =TRUE) == TRUE))
    
  )

title <- odds_ratios.no_weight.2nd_line[1,]
title[,] <- NA


interim_2 <- title %>% as.data.frame() %>%mutate(group = paste0("Predicted benefit on ", drug1, ":"), label = paste0("Predicted benefit on ", drug1, ":")) %>%
  rbind(
    odds_ratios.no_weight.2nd_line %>%
      filter(drug_1 == drug1 & drug_2 == drug2) %>%
      slice(which(grepl("(-", odds_ratios.no_weight.2nd_line$group, fixed =TRUE) == TRUE)),
    title %>% as.data.frame() %>% mutate(group = paste0("Predicted benefit on ", drug2, ":"), label = paste0("Predicted benefit on ", drug2, ":")),
    odds_ratios.no_weight.2nd_line %>%
      filter(drug_1 == drug1 & drug_2 == drug2) %>%
      slice(which(!grepl("(-", odds_ratios.no_weight.2nd_line$group, fixed =TRUE) == TRUE))
    
  )

title <- absolute.var_adj.no_adj_2n_line[1,]
title[,] <- NA

interim_3 <- title %>% as.data.frame() %>% mutate(group = paste0("Predicted benefit on ", drug1, ":"), label = paste0("Predicted benefit on ", drug1, ":")) %>%
  mutate(drug_type = drug1) %>%
  rbind(
    absolute.var_adj.no_adj_2n_line %>%
      filter(drug_1 == drug1 & drug_2 == drug2) %>%
      slice(which(grepl("(-", absolute.var_adj.no_adj_2n_line$group, fixed =TRUE) == TRUE)),
    title %>% as.data.frame() %>% mutate(group = paste0("Predicted benefit on ", drug2, ":"), label = paste0("Predicted benefit on ", drug2, ":")) %>%
      mutate(drug_type = drug2),
    absolute.var_adj.no_adj_2n_line %>%
      filter(drug_1 == drug1 & drug_2 == drug2) %>%
      slice(which(!grepl("(-", absolute.var_adj.no_adj_2n_line$group, fixed =TRUE) == TRUE))
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
  
) + plot_layout(axes = "collect") +
  plot_annotation(
    title = paste0(drug1, " vs ", drug2)
  ) &
  theme(
    axis.title.y = element_blank(),
    axis.title.x = element_text(size = 13),
    axis.text = element_text(size = 13)
  )



##: Per drug combination: GLP1vsTZD
drug1 = "GLP1"; drug2 = "TZD"

title <- ATE.var_adj.no_adj_2nd_line[1,]
title[,] <- NA

interim_1 <- title %>% as.data.frame() %>%mutate(group = paste0("Predicted benefit on ", drug1, ":"), label = paste0("Predicted benefit on ", drug1, ":")) %>%
  rbind(
    ATE.var_adj.no_adj_2nd_line %>%
      filter(drug_1 == drug1 & drug_2 == drug2) %>%
      slice(which(grepl("(-", ATE.var_adj.no_adj_2nd_line$group, fixed =TRUE) == TRUE)),
    title %>% as.data.frame() %>% mutate(group = paste0("Predicted benefit on ", drug2, ":"), label = paste0("Predicted benefit on ", drug2, ":")),
    ATE.var_adj.no_adj_2nd_line %>%
      filter(drug_1 == drug1 & drug_2 == drug2) %>%
      slice(which(!grepl("(-", ATE.var_adj.no_adj_2nd_line$group, fixed =TRUE) == TRUE))
    
  )

title <- odds_ratios.no_weight.2nd_line[1,]
title[,] <- NA


interim_2 <- title %>% as.data.frame() %>%mutate(group = paste0("Predicted benefit on ", drug1, ":"), label = paste0("Predicted benefit on ", drug1, ":")) %>%
  rbind(
    odds_ratios.no_weight.2nd_line %>%
      filter(drug_1 == drug1 & drug_2 == drug2) %>%
      slice(which(grepl("(-", odds_ratios.no_weight.2nd_line$group, fixed =TRUE) == TRUE)),
    title %>% as.data.frame() %>% mutate(group = paste0("Predicted benefit on ", drug2, ":"), label = paste0("Predicted benefit on ", drug2, ":")),
    odds_ratios.no_weight.2nd_line %>%
      filter(drug_1 == drug1 & drug_2 == drug2) %>%
      slice(which(!grepl("(-", odds_ratios.no_weight.2nd_line$group, fixed =TRUE) == TRUE))
    
  )

title <- absolute.var_adj.no_adj_2n_line[1,]
title[,] <- NA

interim_3 <- title %>% as.data.frame() %>% mutate(group = paste0("Predicted benefit on ", drug1, ":"), label = paste0("Predicted benefit on ", drug1, ":")) %>%
  mutate(drug_type = drug1) %>%
  rbind(
    absolute.var_adj.no_adj_2n_line %>%
      filter(drug_1 == drug1 & drug_2 == drug2) %>%
      slice(which(grepl("(-", absolute.var_adj.no_adj_2n_line$group, fixed =TRUE) == TRUE)),
    title %>% as.data.frame() %>% mutate(group = paste0("Predicted benefit on ", drug2, ":"), label = paste0("Predicted benefit on ", drug2, ":")) %>%
      mutate(drug_type = drug2),
    absolute.var_adj.no_adj_2n_line %>%
      filter(drug_1 == drug1 & drug_2 == drug2) %>%
      slice(which(!grepl("(-", absolute.var_adj.no_adj_2n_line$group, fixed =TRUE) == TRUE))
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
  
) + plot_layout(axes = "collect") +
  plot_annotation(
    title = paste0(drug1, " vs ", drug2)
  ) &
  theme(
    axis.title.y = element_blank(),
    axis.title.x = element_text(size = 13),
    axis.text = element_text(size = 13)
  )


##: Per drug combination: SGLT2vsSU
drug1 = "SGLT2"; drug2 = "SU"

title <- ATE.var_adj.no_adj_2nd_line[1,]
title[,] <- NA

interim_1 <- title %>% as.data.frame() %>%mutate(group = paste0("Predicted benefit on ", drug1, ":"), label = paste0("Predicted benefit on ", drug1, ":")) %>%
  rbind(
    ATE.var_adj.no_adj_2nd_line %>%
      filter(drug_1 == drug1 & drug_2 == drug2) %>%
      slice(which(grepl("(-", ATE.var_adj.no_adj_2nd_line$group, fixed =TRUE) == TRUE)),
    title %>% as.data.frame() %>% mutate(group = paste0("Predicted benefit on ", drug2, ":"), label = paste0("Predicted benefit on ", drug2, ":")),
    ATE.var_adj.no_adj_2nd_line %>%
      filter(drug_1 == drug1 & drug_2 == drug2) %>%
      slice(which(!grepl("(-", ATE.var_adj.no_adj_2nd_line$group, fixed =TRUE) == TRUE))
    
  )

title <- odds_ratios.no_weight.2nd_line[1,]
title[,] <- NA


interim_2 <- title %>% as.data.frame() %>%mutate(group = paste0("Predicted benefit on ", drug1, ":"), label = paste0("Predicted benefit on ", drug1, ":")) %>%
  rbind(
    odds_ratios.no_weight.2nd_line %>%
      filter(drug_1 == drug1 & drug_2 == drug2) %>%
      slice(which(grepl("(-", odds_ratios.no_weight.2nd_line$group, fixed =TRUE) == TRUE)),
    title %>% as.data.frame() %>% mutate(group = paste0("Predicted benefit on ", drug2, ":"), label = paste0("Predicted benefit on ", drug2, ":")),
    odds_ratios.no_weight.2nd_line %>%
      filter(drug_1 == drug1 & drug_2 == drug2) %>%
      slice(which(!grepl("(-", odds_ratios.no_weight.2nd_line$group, fixed =TRUE) == TRUE))
    
  )

title <- absolute.var_adj.no_adj_2n_line[1,]
title[,] <- NA

interim_3 <- title %>% as.data.frame() %>% mutate(group = paste0("Predicted benefit on ", drug1, ":"), label = paste0("Predicted benefit on ", drug1, ":")) %>%
  mutate(drug_type = drug1) %>%
  rbind(
    absolute.var_adj.no_adj_2n_line %>%
      filter(drug_1 == drug1 & drug_2 == drug2) %>%
      slice(which(grepl("(-", absolute.var_adj.no_adj_2n_line$group, fixed =TRUE) == TRUE)),
    title %>% as.data.frame() %>% mutate(group = paste0("Predicted benefit on ", drug2, ":"), label = paste0("Predicted benefit on ", drug2, ":")) %>%
      mutate(drug_type = drug2),
    absolute.var_adj.no_adj_2n_line %>%
      filter(drug_1 == drug1 & drug_2 == drug2) %>%
      slice(which(!grepl("(-", absolute.var_adj.no_adj_2n_line$group, fixed =TRUE) == TRUE))
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
  
) + plot_layout(axes = "collect") +
  plot_annotation(
    title = paste0(drug1, " vs ", drug2)
  ) &
  theme(
    axis.title.y = element_blank(),
    axis.title.x = element_text(size = 13),
    axis.text = element_text(size = 13)
  )



##: Per drug combination: SGLT2vsTZD
drug1 = "SGLT2"; drug2 = "TZD"

title <- ATE.var_adj.no_adj_2nd_line[1,]
title[,] <- NA

interim_1 <- title %>% as.data.frame() %>%mutate(group = paste0("Predicted benefit on ", drug1, ":"), label = paste0("Predicted benefit on ", drug1, ":")) %>%
  rbind(
    ATE.var_adj.no_adj_2nd_line %>%
      filter(drug_1 == drug1 & drug_2 == drug2) %>%
      slice(which(grepl("(-", ATE.var_adj.no_adj_2nd_line$group, fixed =TRUE) == TRUE)),
    title %>% as.data.frame() %>% mutate(group = paste0("Predicted benefit on ", drug2, ":"), label = paste0("Predicted benefit on ", drug2, ":")),
    ATE.var_adj.no_adj_2nd_line %>%
      filter(drug_1 == drug1 & drug_2 == drug2) %>%
      slice(which(!grepl("(-", ATE.var_adj.no_adj_2nd_line$group, fixed =TRUE) == TRUE))
    
  )

title <- odds_ratios.no_weight.2nd_line[1,]
title[,] <- NA


interim_2 <- title %>% as.data.frame() %>%mutate(group = paste0("Predicted benefit on ", drug1, ":"), label = paste0("Predicted benefit on ", drug1, ":")) %>%
  rbind(
    odds_ratios.no_weight.2nd_line %>%
      filter(drug_1 == drug1 & drug_2 == drug2) %>%
      slice(which(grepl("(-", odds_ratios.no_weight.2nd_line$group, fixed =TRUE) == TRUE)),
    title %>% as.data.frame() %>% mutate(group = paste0("Predicted benefit on ", drug2, ":"), label = paste0("Predicted benefit on ", drug2, ":")),
    odds_ratios.no_weight.2nd_line %>%
      filter(drug_1 == drug1 & drug_2 == drug2) %>%
      slice(which(!grepl("(-", odds_ratios.no_weight.2nd_line$group, fixed =TRUE) == TRUE))
    
  )

title <- absolute.var_adj.no_adj_2n_line[1,]
title[,] <- NA

interim_3 <- title %>% as.data.frame() %>% mutate(group = paste0("Predicted benefit on ", drug1, ":"), label = paste0("Predicted benefit on ", drug1, ":")) %>%
  mutate(drug_type = drug1) %>%
  rbind(
    absolute.var_adj.no_adj_2n_line %>%
      filter(drug_1 == drug1 & drug_2 == drug2) %>%
      slice(which(grepl("(-", absolute.var_adj.no_adj_2n_line$group, fixed =TRUE) == TRUE)),
    title %>% as.data.frame() %>% mutate(group = paste0("Predicted benefit on ", drug2, ":"), label = paste0("Predicted benefit on ", drug2, ":")) %>%
      mutate(drug_type = drug2),
    absolute.var_adj.no_adj_2n_line %>%
      filter(drug_1 == drug1 & drug_2 == drug2) %>%
      slice(which(!grepl("(-", absolute.var_adj.no_adj_2n_line$group, fixed =TRUE) == TRUE))
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
  
) + plot_layout(axes = "collect") +
  plot_annotation(
    title = paste0(drug1, " vs ", drug2)
  ) &
  theme(
    axis.title.y = element_blank(),
    axis.title.x = element_text(size = 13),
    axis.text = element_text(size = 13)
  )




##: Per drug combination: SUvsTZD
drug1 = "SU"; drug2 = "TZD"

title <- ATE.var_adj.no_adj_2nd_line[1,]
title[,] <- NA

interim_1 <- title %>% as.data.frame() %>%mutate(group = paste0("Predicted benefit on ", drug1, ":"), label = paste0("Predicted benefit on ", drug1, ":")) %>%
  rbind(
    ATE.var_adj.no_adj_2nd_line %>%
      filter(drug_1 == drug1 & drug_2 == drug2) %>%
      slice(which(grepl("(-", ATE.var_adj.no_adj_2nd_line$group, fixed =TRUE) == TRUE)),
    title %>% as.data.frame() %>% mutate(group = paste0("Predicted benefit on ", drug2, ":"), label = paste0("Predicted benefit on ", drug2, ":")),
    ATE.var_adj.no_adj_2nd_line %>%
      filter(drug_1 == drug1 & drug_2 == drug2) %>%
      slice(which(!grepl("(-", ATE.var_adj.no_adj_2nd_line$group, fixed =TRUE) == TRUE))
    
  )

title <- odds_ratios.no_weight.2nd_line[1,]
title[,] <- NA


interim_2 <- title %>% as.data.frame() %>%mutate(group = paste0("Predicted benefit on ", drug1, ":"), label = paste0("Predicted benefit on ", drug1, ":")) %>%
  rbind(
    odds_ratios.no_weight.2nd_line %>%
      filter(drug_1 == drug1 & drug_2 == drug2) %>%
      slice(which(grepl("(-", odds_ratios.no_weight.2nd_line$group, fixed =TRUE) == TRUE)),
    title %>% as.data.frame() %>% mutate(group = paste0("Predicted benefit on ", drug2, ":"), label = paste0("Predicted benefit on ", drug2, ":")),
    odds_ratios.no_weight.2nd_line %>%
      filter(drug_1 == drug1 & drug_2 == drug2) %>%
      slice(which(!grepl("(-", odds_ratios.no_weight.2nd_line$group, fixed =TRUE) == TRUE))
    
  )

title <- absolute.var_adj.no_adj_2n_line[1,]
title[,] <- NA

interim_3 <- title %>% as.data.frame() %>% mutate(group = paste0("Predicted benefit on ", drug1, ":"), label = paste0("Predicted benefit on ", drug1, ":")) %>%
  mutate(drug_type = drug1) %>%
  rbind(
    absolute.var_adj.no_adj_2n_line %>%
      filter(drug_1 == drug1 & drug_2 == drug2) %>%
      slice(which(grepl("(-", absolute.var_adj.no_adj_2n_line$group, fixed =TRUE) == TRUE)),
    title %>% as.data.frame() %>% mutate(group = paste0("Predicted benefit on ", drug2, ":"), label = paste0("Predicted benefit on ", drug2, ":")) %>%
      mutate(drug_type = drug2),
    absolute.var_adj.no_adj_2n_line %>%
      filter(drug_1 == drug1 & drug_2 == drug2) %>%
      slice(which(!grepl("(-", absolute.var_adj.no_adj_2n_line$group, fixed =TRUE) == TRUE))
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
  
) + plot_layout(axes = "collect") +
  plot_annotation(
    title = paste0(drug1, " vs ", drug2)
  ) &
  theme(
    axis.title.y = element_blank(),
    axis.title.x = element_text(size = 13),
    axis.text = element_text(size = 13)
  )



dev.off()





#:----------------------------
#:----------------------------
#:----------------------------


# Combining heterogeneity and odds ratios

breakpoints <- c(-0.1, -0.05, -0.025, 0, 0.025, 0.05, 0.1)

# ATE not adjusted
ATE.var_adj.no_adj_2nd_line <- calc_ATE(cprd_dataset.val, break_points = breakpoints, drugs = c("DPP4", "GLP1", "SGLT2", "SU", "TZD"), pred.variable = ".no_weight", n_bootstrap = 10) %>%
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
absolute.var_adj.no_adj_2n_line <- calc_predicted_risk(cprd_dataset.val, break_points = breakpoints, drugs = c("DPP4", "GLP1", "SGLT2", "SU", "TZD"), pred.variable = ".no_weight") %>%
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
odds_ratios.no_weight.2nd_line <- calc_odds_ratios(cprd_dataset.val, drugs = c("DPP4", "GLP1", "SGLT2", "SU", "TZD"), pred.variable = ".no_weight", break_points = breakpoints) %>%
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



pdf("results/figures/absolute_relative_validation.pdf", width = 20, height = 7)

##: Per drug combination: DPP4vsGLP1
drug1 = "DPP4"; drug2 = "GLP1"

title <- ATE.var_adj.no_adj_2nd_line[1,]
title[,] <- NA

interim_1 <- title %>% as.data.frame() %>%mutate(group = paste0("Predicted benefit on ", drug1, ":"), label = paste0("Predicted benefit on ", drug1, ":")) %>%
  rbind(
    ATE.var_adj.no_adj_2nd_line %>%
      filter(drug_1 == drug1 & drug_2 == drug2) %>%
      slice(which(grepl("(-", ATE.var_adj.no_adj_2nd_line$group, fixed =TRUE) == TRUE)),
    title %>% as.data.frame() %>% mutate(group = paste0("Predicted benefit on ", drug2, ":"), label = paste0("Predicted benefit on ", drug2, ":")),
    ATE.var_adj.no_adj_2nd_line %>%
      filter(drug_1 == drug1 & drug_2 == drug2) %>%
      slice(which(!grepl("(-", ATE.var_adj.no_adj_2nd_line$group, fixed =TRUE) == TRUE))
    
  )

title <- odds_ratios.no_weight.2nd_line[1,]
title[,] <- NA


interim_2 <- title %>% as.data.frame() %>%mutate(group = paste0("Predicted benefit on ", drug1, ":"), label = paste0("Predicted benefit on ", drug1, ":")) %>%
  rbind(
    odds_ratios.no_weight.2nd_line %>%
      filter(drug_1 == drug1 & drug_2 == drug2) %>%
      slice(which(grepl("(-", odds_ratios.no_weight.2nd_line$group, fixed =TRUE) == TRUE)),
    title %>% as.data.frame() %>% mutate(group = paste0("Predicted benefit on ", drug2, ":"), label = paste0("Predicted benefit on ", drug2, ":")),
    odds_ratios.no_weight.2nd_line %>%
      filter(drug_1 == drug1 & drug_2 == drug2) %>%
      slice(which(!grepl("(-", odds_ratios.no_weight.2nd_line$group, fixed =TRUE) == TRUE))
    
  )

title <- absolute.var_adj.no_adj_2n_line[1,]
title[,] <- NA

interim_3 <- title %>% as.data.frame() %>% mutate(group = paste0("Predicted benefit on ", drug1, ":"), label = paste0("Predicted benefit on ", drug1, ":")) %>%
  mutate(drug_type = drug1) %>%
  rbind(
    absolute.var_adj.no_adj_2n_line %>%
      filter(drug_1 == drug1 & drug_2 == drug2) %>%
      slice(which(grepl("(-", absolute.var_adj.no_adj_2n_line$group, fixed =TRUE) == TRUE)),
    title %>% as.data.frame() %>% mutate(group = paste0("Predicted benefit on ", drug2, ":"), label = paste0("Predicted benefit on ", drug2, ":")) %>%
      mutate(drug_type = drug2),
    absolute.var_adj.no_adj_2n_line %>%
      filter(drug_1 == drug1 & drug_2 == drug2) %>%
      slice(which(!grepl("(-", absolute.var_adj.no_adj_2n_line$group, fixed =TRUE) == TRUE))
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
  
) + plot_layout(axes = "collect") +
  plot_annotation(
    title = paste0(drug1, " vs ", drug2)
  ) &
  theme(
    axis.title.y = element_blank(),
    axis.title.x = element_text(size = 13),
    axis.text = element_text(size = 13)
  )


##: Per drug combination: DPP4vsSGLT2
drug1 = "DPP4"; drug2 = "SGLT2"

title <- ATE.var_adj.no_adj_2nd_line[1,]
title[,] <- NA

interim_1 <- title %>% as.data.frame() %>%mutate(group = paste0("Predicted benefit on ", drug1, ":"), label = paste0("Predicted benefit on ", drug1, ":")) %>%
  rbind(
    ATE.var_adj.no_adj_2nd_line %>%
      filter(drug_1 == drug1 & drug_2 == drug2) %>%
      slice(which(grepl("(-", ATE.var_adj.no_adj_2nd_line$group, fixed =TRUE) == TRUE)),
    title %>% as.data.frame() %>% mutate(group = paste0("Predicted benefit on ", drug2, ":"), label = paste0("Predicted benefit on ", drug2, ":")),
    ATE.var_adj.no_adj_2nd_line %>%
      filter(drug_1 == drug1 & drug_2 == drug2) %>%
      slice(which(!grepl("(-", ATE.var_adj.no_adj_2nd_line$group, fixed =TRUE) == TRUE))
    
  )

title <- odds_ratios.no_weight.2nd_line[1,]
title[,] <- NA


interim_2 <- title %>% as.data.frame() %>%mutate(group = paste0("Predicted benefit on ", drug1, ":"), label = paste0("Predicted benefit on ", drug1, ":")) %>%
  rbind(
    odds_ratios.no_weight.2nd_line %>%
      filter(drug_1 == drug1 & drug_2 == drug2) %>%
      slice(which(grepl("(-", odds_ratios.no_weight.2nd_line$group, fixed =TRUE) == TRUE)),
    title %>% as.data.frame() %>% mutate(group = paste0("Predicted benefit on ", drug2, ":"), label = paste0("Predicted benefit on ", drug2, ":")),
    odds_ratios.no_weight.2nd_line %>%
      filter(drug_1 == drug1 & drug_2 == drug2) %>%
      slice(which(!grepl("(-", odds_ratios.no_weight.2nd_line$group, fixed =TRUE) == TRUE))
    
  )

title <- absolute.var_adj.no_adj_2n_line[1,]
title[,] <- NA

interim_3 <- title %>% as.data.frame() %>% mutate(group = paste0("Predicted benefit on ", drug1, ":"), label = paste0("Predicted benefit on ", drug1, ":")) %>%
  mutate(drug_type = drug1) %>%
  rbind(
    absolute.var_adj.no_adj_2n_line %>%
      filter(drug_1 == drug1 & drug_2 == drug2) %>%
      slice(which(grepl("(-", absolute.var_adj.no_adj_2n_line$group, fixed =TRUE) == TRUE)),
    title %>% as.data.frame() %>% mutate(group = paste0("Predicted benefit on ", drug2, ":"), label = paste0("Predicted benefit on ", drug2, ":")) %>%
      mutate(drug_type = drug2),
    absolute.var_adj.no_adj_2n_line %>%
      filter(drug_1 == drug1 & drug_2 == drug2) %>%
      slice(which(!grepl("(-", absolute.var_adj.no_adj_2n_line$group, fixed =TRUE) == TRUE))
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
    coord_cartesian(xlim = c(-0.3, 0.1)) +
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
    coord_cartesian(xlim = c(0.05, 0.4)) +
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
  
) + plot_layout(axes = "collect") +
  plot_annotation(
    title = paste0(drug1, " vs ", drug2)
  ) &
  theme(
    axis.title.y = element_blank(),
    axis.title.x = element_text(size = 13),
    axis.text = element_text(size = 13)
  )


##: Per drug combination: DPP4vsSU
drug1 = "DPP4"; drug2 = "SU"

title <- ATE.var_adj.no_adj_2nd_line[1,]
title[,] <- NA

interim_1 <- title %>% as.data.frame() %>%mutate(group = paste0("Predicted benefit on ", drug1, ":"), label = paste0("Predicted benefit on ", drug1, ":")) %>%
  rbind(
    ATE.var_adj.no_adj_2nd_line %>%
      filter(drug_1 == drug1 & drug_2 == drug2) %>%
      slice(which(grepl("(-", ATE.var_adj.no_adj_2nd_line$group, fixed =TRUE) == TRUE)),
    title %>% as.data.frame() %>% mutate(group = paste0("Predicted benefit on ", drug2, ":"), label = paste0("Predicted benefit on ", drug2, ":")),
    ATE.var_adj.no_adj_2nd_line %>%
      filter(drug_1 == drug1 & drug_2 == drug2) %>%
      slice(which(!grepl("(-", ATE.var_adj.no_adj_2nd_line$group, fixed =TRUE) == TRUE))
    
  )

title <- odds_ratios.no_weight.2nd_line[1,]
title[,] <- NA


interim_2 <- title %>% as.data.frame() %>%mutate(group = paste0("Predicted benefit on ", drug1, ":"), label = paste0("Predicted benefit on ", drug1, ":")) %>%
  rbind(
    odds_ratios.no_weight.2nd_line %>%
      filter(drug_1 == drug1 & drug_2 == drug2) %>%
      slice(which(grepl("(-", odds_ratios.no_weight.2nd_line$group, fixed =TRUE) == TRUE)),
    title %>% as.data.frame() %>% mutate(group = paste0("Predicted benefit on ", drug2, ":"), label = paste0("Predicted benefit on ", drug2, ":")),
    odds_ratios.no_weight.2nd_line %>%
      filter(drug_1 == drug1 & drug_2 == drug2) %>%
      slice(which(!grepl("(-", odds_ratios.no_weight.2nd_line$group, fixed =TRUE) == TRUE))
    
  )

title <- absolute.var_adj.no_adj_2n_line[1,]
title[,] <- NA

interim_3 <- title %>% as.data.frame() %>% mutate(group = paste0("Predicted benefit on ", drug1, ":"), label = paste0("Predicted benefit on ", drug1, ":")) %>%
  mutate(drug_type = drug1) %>%
  rbind(
    absolute.var_adj.no_adj_2n_line %>%
      filter(drug_1 == drug1 & drug_2 == drug2) %>%
      slice(which(grepl("(-", absolute.var_adj.no_adj_2n_line$group, fixed =TRUE) == TRUE)),
    title %>% as.data.frame() %>% mutate(group = paste0("Predicted benefit on ", drug2, ":"), label = paste0("Predicted benefit on ", drug2, ":")) %>%
      mutate(drug_type = drug2),
    absolute.var_adj.no_adj_2n_line %>%
      filter(drug_1 == drug1 & drug_2 == drug2) %>%
      slice(which(!grepl("(-", absolute.var_adj.no_adj_2n_line$group, fixed =TRUE) == TRUE))
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
  
) + plot_layout(axes = "collect") +
  plot_annotation(
    title = paste0(drug1, " vs ", drug2)
  ) &
  theme(
    axis.title.y = element_blank(),
    axis.title.x = element_text(size = 13),
    axis.text = element_text(size = 13)
  )


##: Per drug combination: DPP4vsTZD
drug1 = "DPP4"; drug2 = "TZD"

title <- ATE.var_adj.no_adj_2nd_line[1,]
title[,] <- NA

interim_1 <- title %>% as.data.frame() %>%mutate(group = paste0("Predicted benefit on ", drug1, ":"), label = paste0("Predicted benefit on ", drug1, ":")) %>%
  rbind(
    ATE.var_adj.no_adj_2nd_line %>%
      filter(drug_1 == drug1 & drug_2 == drug2) %>%
      slice(which(grepl("(-", ATE.var_adj.no_adj_2nd_line$group, fixed =TRUE) == TRUE)),
    title %>% as.data.frame() %>% mutate(group = paste0("Predicted benefit on ", drug2, ":"), label = paste0("Predicted benefit on ", drug2, ":")),
    ATE.var_adj.no_adj_2nd_line %>%
      filter(drug_1 == drug1 & drug_2 == drug2) %>%
      slice(which(!grepl("(-", ATE.var_adj.no_adj_2nd_line$group, fixed =TRUE) == TRUE))
    
  )

title <- odds_ratios.no_weight.2nd_line[1,]
title[,] <- NA


interim_2 <- title %>% as.data.frame() %>%mutate(group = paste0("Predicted benefit on ", drug1, ":"), label = paste0("Predicted benefit on ", drug1, ":")) %>%
  rbind(
    odds_ratios.no_weight.2nd_line %>%
      filter(drug_1 == drug1 & drug_2 == drug2) %>%
      slice(which(grepl("(-", odds_ratios.no_weight.2nd_line$group, fixed =TRUE) == TRUE)),
    title %>% as.data.frame() %>% mutate(group = paste0("Predicted benefit on ", drug2, ":"), label = paste0("Predicted benefit on ", drug2, ":")),
    odds_ratios.no_weight.2nd_line %>%
      filter(drug_1 == drug1 & drug_2 == drug2) %>%
      slice(which(!grepl("(-", odds_ratios.no_weight.2nd_line$group, fixed =TRUE) == TRUE))
    
  )

title <- absolute.var_adj.no_adj_2n_line[1,]
title[,] <- NA

interim_3 <- title %>% as.data.frame() %>% mutate(group = paste0("Predicted benefit on ", drug1, ":"), label = paste0("Predicted benefit on ", drug1, ":")) %>%
  mutate(drug_type = drug1) %>%
  rbind(
    absolute.var_adj.no_adj_2n_line %>%
      filter(drug_1 == drug1 & drug_2 == drug2) %>%
      slice(which(grepl("(-", absolute.var_adj.no_adj_2n_line$group, fixed =TRUE) == TRUE)),
    title %>% as.data.frame() %>% mutate(group = paste0("Predicted benefit on ", drug2, ":"), label = paste0("Predicted benefit on ", drug2, ":")) %>%
      mutate(drug_type = drug2),
    absolute.var_adj.no_adj_2n_line %>%
      filter(drug_1 == drug1 & drug_2 == drug2) %>%
      slice(which(!grepl("(-", absolute.var_adj.no_adj_2n_line$group, fixed =TRUE) == TRUE))
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
  
) + plot_layout(axes = "collect") +
  plot_annotation(
    title = paste0(drug1, " vs ", drug2)
  ) &
  theme(
    axis.title.y = element_blank(),
    axis.title.x = element_text(size = 13),
    axis.text = element_text(size = 13)
  )


##: Per drug combination: GLP1vsSGLT2
drug1 = "GLP1"; drug2 = "SGLT2"

title <- ATE.var_adj.no_adj_2nd_line[1,]
title[,] <- NA

interim_1 <- title %>% as.data.frame() %>%mutate(group = paste0("Predicted benefit on ", drug1, ":"), label = paste0("Predicted benefit on ", drug1, ":")) %>%
  rbind(
    ATE.var_adj.no_adj_2nd_line %>%
      filter(drug_1 == drug1 & drug_2 == drug2) %>%
      slice(which(grepl("(-", ATE.var_adj.no_adj_2nd_line$group, fixed =TRUE) == TRUE)),
    title %>% as.data.frame() %>% mutate(group = paste0("Predicted benefit on ", drug2, ":"), label = paste0("Predicted benefit on ", drug2, ":")),
    ATE.var_adj.no_adj_2nd_line %>%
      filter(drug_1 == drug1 & drug_2 == drug2) %>%
      slice(which(!grepl("(-", ATE.var_adj.no_adj_2nd_line$group, fixed =TRUE) == TRUE))
    
  )

title <- odds_ratios.no_weight.2nd_line[1,]
title[,] <- NA


interim_2 <- title %>% as.data.frame() %>%mutate(group = paste0("Predicted benefit on ", drug1, ":"), label = paste0("Predicted benefit on ", drug1, ":")) %>%
  rbind(
    odds_ratios.no_weight.2nd_line %>%
      filter(drug_1 == drug1 & drug_2 == drug2) %>%
      slice(which(grepl("(-", odds_ratios.no_weight.2nd_line$group, fixed =TRUE) == TRUE)),
    title %>% as.data.frame() %>% mutate(group = paste0("Predicted benefit on ", drug2, ":"), label = paste0("Predicted benefit on ", drug2, ":")),
    odds_ratios.no_weight.2nd_line %>%
      filter(drug_1 == drug1 & drug_2 == drug2) %>%
      slice(which(!grepl("(-", odds_ratios.no_weight.2nd_line$group, fixed =TRUE) == TRUE))
    
  )

title <- absolute.var_adj.no_adj_2n_line[1,]
title[,] <- NA

interim_3 <- title %>% as.data.frame() %>% mutate(group = paste0("Predicted benefit on ", drug1, ":"), label = paste0("Predicted benefit on ", drug1, ":")) %>%
  mutate(drug_type = drug1) %>%
  rbind(
    absolute.var_adj.no_adj_2n_line %>%
      filter(drug_1 == drug1 & drug_2 == drug2) %>%
      slice(which(grepl("(-", absolute.var_adj.no_adj_2n_line$group, fixed =TRUE) == TRUE)),
    title %>% as.data.frame() %>% mutate(group = paste0("Predicted benefit on ", drug2, ":"), label = paste0("Predicted benefit on ", drug2, ":")) %>%
      mutate(drug_type = drug2),
    absolute.var_adj.no_adj_2n_line %>%
      filter(drug_1 == drug1 & drug_2 == drug2) %>%
      slice(which(!grepl("(-", absolute.var_adj.no_adj_2n_line$group, fixed =TRUE) == TRUE))
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
  
) + plot_layout(axes = "collect") +
  plot_annotation(
    title = paste0(drug1, " vs ", drug2)
  ) &
  theme(
    axis.title.y = element_blank(),
    axis.title.x = element_text(size = 13),
    axis.text = element_text(size = 13)
  )



##: Per drug combination: GLP1vsSU
drug1 = "GLP1"; drug2 = "SU"

title <- ATE.var_adj.no_adj_2nd_line[1,]
title[,] <- NA

interim_1 <- title %>% as.data.frame() %>%mutate(group = paste0("Predicted benefit on ", drug1, ":"), label = paste0("Predicted benefit on ", drug1, ":")) %>%
  rbind(
    ATE.var_adj.no_adj_2nd_line %>%
      filter(drug_1 == drug1 & drug_2 == drug2) %>%
      slice(which(grepl("(-", ATE.var_adj.no_adj_2nd_line$group, fixed =TRUE) == TRUE)),
    title %>% as.data.frame() %>% mutate(group = paste0("Predicted benefit on ", drug2, ":"), label = paste0("Predicted benefit on ", drug2, ":")),
    ATE.var_adj.no_adj_2nd_line %>%
      filter(drug_1 == drug1 & drug_2 == drug2) %>%
      slice(which(!grepl("(-", ATE.var_adj.no_adj_2nd_line$group, fixed =TRUE) == TRUE))
    
  )

title <- odds_ratios.no_weight.2nd_line[1,]
title[,] <- NA


interim_2 <- title %>% as.data.frame() %>%mutate(group = paste0("Predicted benefit on ", drug1, ":"), label = paste0("Predicted benefit on ", drug1, ":")) %>%
  rbind(
    odds_ratios.no_weight.2nd_line %>%
      filter(drug_1 == drug1 & drug_2 == drug2) %>%
      slice(which(grepl("(-", odds_ratios.no_weight.2nd_line$group, fixed =TRUE) == TRUE)),
    title %>% as.data.frame() %>% mutate(group = paste0("Predicted benefit on ", drug2, ":"), label = paste0("Predicted benefit on ", drug2, ":")),
    odds_ratios.no_weight.2nd_line %>%
      filter(drug_1 == drug1 & drug_2 == drug2) %>%
      slice(which(!grepl("(-", odds_ratios.no_weight.2nd_line$group, fixed =TRUE) == TRUE))
    
  )

title <- absolute.var_adj.no_adj_2n_line[1,]
title[,] <- NA

interim_3 <- title %>% as.data.frame() %>% mutate(group = paste0("Predicted benefit on ", drug1, ":"), label = paste0("Predicted benefit on ", drug1, ":")) %>%
  mutate(drug_type = drug1) %>%
  rbind(
    absolute.var_adj.no_adj_2n_line %>%
      filter(drug_1 == drug1 & drug_2 == drug2) %>%
      slice(which(grepl("(-", absolute.var_adj.no_adj_2n_line$group, fixed =TRUE) == TRUE)),
    title %>% as.data.frame() %>% mutate(group = paste0("Predicted benefit on ", drug2, ":"), label = paste0("Predicted benefit on ", drug2, ":")) %>%
      mutate(drug_type = drug2),
    absolute.var_adj.no_adj_2n_line %>%
      filter(drug_1 == drug1 & drug_2 == drug2) %>%
      slice(which(!grepl("(-", absolute.var_adj.no_adj_2n_line$group, fixed =TRUE) == TRUE))
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
  
) + plot_layout(axes = "collect") +
  plot_annotation(
    title = paste0(drug1, " vs ", drug2)
  ) &
  theme(
    axis.title.y = element_blank(),
    axis.title.x = element_text(size = 13),
    axis.text = element_text(size = 13)
  )



##: Per drug combination: GLP1vsTZD
drug1 = "GLP1"; drug2 = "TZD"

title <- ATE.var_adj.no_adj_2nd_line[1,]
title[,] <- NA

interim_1 <- title %>% as.data.frame() %>%mutate(group = paste0("Predicted benefit on ", drug1, ":"), label = paste0("Predicted benefit on ", drug1, ":")) %>%
  rbind(
    ATE.var_adj.no_adj_2nd_line %>%
      filter(drug_1 == drug1 & drug_2 == drug2) %>%
      slice(which(grepl("(-", ATE.var_adj.no_adj_2nd_line$group, fixed =TRUE) == TRUE)),
    title %>% as.data.frame() %>% mutate(group = paste0("Predicted benefit on ", drug2, ":"), label = paste0("Predicted benefit on ", drug2, ":")),
    ATE.var_adj.no_adj_2nd_line %>%
      filter(drug_1 == drug1 & drug_2 == drug2) %>%
      slice(which(!grepl("(-", ATE.var_adj.no_adj_2nd_line$group, fixed =TRUE) == TRUE))
    
  )

title <- odds_ratios.no_weight.2nd_line[1,]
title[,] <- NA


interim_2 <- title %>% as.data.frame() %>%mutate(group = paste0("Predicted benefit on ", drug1, ":"), label = paste0("Predicted benefit on ", drug1, ":")) %>%
  rbind(
    odds_ratios.no_weight.2nd_line %>%
      filter(drug_1 == drug1 & drug_2 == drug2) %>%
      slice(which(grepl("(-", odds_ratios.no_weight.2nd_line$group, fixed =TRUE) == TRUE)),
    title %>% as.data.frame() %>% mutate(group = paste0("Predicted benefit on ", drug2, ":"), label = paste0("Predicted benefit on ", drug2, ":")),
    odds_ratios.no_weight.2nd_line %>%
      filter(drug_1 == drug1 & drug_2 == drug2) %>%
      slice(which(!grepl("(-", odds_ratios.no_weight.2nd_line$group, fixed =TRUE) == TRUE))
    
  )

title <- absolute.var_adj.no_adj_2n_line[1,]
title[,] <- NA

interim_3 <- title %>% as.data.frame() %>% mutate(group = paste0("Predicted benefit on ", drug1, ":"), label = paste0("Predicted benefit on ", drug1, ":")) %>%
  mutate(drug_type = drug1) %>%
  rbind(
    absolute.var_adj.no_adj_2n_line %>%
      filter(drug_1 == drug1 & drug_2 == drug2) %>%
      slice(which(grepl("(-", absolute.var_adj.no_adj_2n_line$group, fixed =TRUE) == TRUE)),
    title %>% as.data.frame() %>% mutate(group = paste0("Predicted benefit on ", drug2, ":"), label = paste0("Predicted benefit on ", drug2, ":")) %>%
      mutate(drug_type = drug2),
    absolute.var_adj.no_adj_2n_line %>%
      filter(drug_1 == drug1 & drug_2 == drug2) %>%
      slice(which(!grepl("(-", absolute.var_adj.no_adj_2n_line$group, fixed =TRUE) == TRUE))
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
  
) + plot_layout(axes = "collect") +
  plot_annotation(
    title = paste0(drug1, " vs ", drug2)
  ) &
  theme(
    axis.title.y = element_blank(),
    axis.title.x = element_text(size = 13),
    axis.text = element_text(size = 13)
  )


##: Per drug combination: SGLT2vsSU
drug1 = "SGLT2"; drug2 = "SU"

title <- ATE.var_adj.no_adj_2nd_line[1,]
title[,] <- NA

interim_1 <- title %>% as.data.frame() %>%mutate(group = paste0("Predicted benefit on ", drug1, ":"), label = paste0("Predicted benefit on ", drug1, ":")) %>%
  rbind(
    ATE.var_adj.no_adj_2nd_line %>%
      filter(drug_1 == drug1 & drug_2 == drug2) %>%
      slice(which(grepl("(-", ATE.var_adj.no_adj_2nd_line$group, fixed =TRUE) == TRUE)),
    title %>% as.data.frame() %>% mutate(group = paste0("Predicted benefit on ", drug2, ":"), label = paste0("Predicted benefit on ", drug2, ":")),
    ATE.var_adj.no_adj_2nd_line %>%
      filter(drug_1 == drug1 & drug_2 == drug2) %>%
      slice(which(!grepl("(-", ATE.var_adj.no_adj_2nd_line$group, fixed =TRUE) == TRUE))
    
  )

title <- odds_ratios.no_weight.2nd_line[1,]
title[,] <- NA


interim_2 <- title %>% as.data.frame() %>%mutate(group = paste0("Predicted benefit on ", drug1, ":"), label = paste0("Predicted benefit on ", drug1, ":")) %>%
  rbind(
    odds_ratios.no_weight.2nd_line %>%
      filter(drug_1 == drug1 & drug_2 == drug2) %>%
      slice(which(grepl("(-", odds_ratios.no_weight.2nd_line$group, fixed =TRUE) == TRUE)),
    title %>% as.data.frame() %>% mutate(group = paste0("Predicted benefit on ", drug2, ":"), label = paste0("Predicted benefit on ", drug2, ":")),
    odds_ratios.no_weight.2nd_line %>%
      filter(drug_1 == drug1 & drug_2 == drug2) %>%
      slice(which(!grepl("(-", odds_ratios.no_weight.2nd_line$group, fixed =TRUE) == TRUE))
    
  )

title <- absolute.var_adj.no_adj_2n_line[1,]
title[,] <- NA

interim_3 <- title %>% as.data.frame() %>% mutate(group = paste0("Predicted benefit on ", drug1, ":"), label = paste0("Predicted benefit on ", drug1, ":")) %>%
  mutate(drug_type = drug1) %>%
  rbind(
    absolute.var_adj.no_adj_2n_line %>%
      filter(drug_1 == drug1 & drug_2 == drug2) %>%
      slice(which(grepl("(-", absolute.var_adj.no_adj_2n_line$group, fixed =TRUE) == TRUE)),
    title %>% as.data.frame() %>% mutate(group = paste0("Predicted benefit on ", drug2, ":"), label = paste0("Predicted benefit on ", drug2, ":")) %>%
      mutate(drug_type = drug2),
    absolute.var_adj.no_adj_2n_line %>%
      filter(drug_1 == drug1 & drug_2 == drug2) %>%
      slice(which(!grepl("(-", absolute.var_adj.no_adj_2n_line$group, fixed =TRUE) == TRUE))
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
  
) + plot_layout(axes = "collect") +
  plot_annotation(
    title = paste0(drug1, " vs ", drug2)
  ) &
  theme(
    axis.title.y = element_blank(),
    axis.title.x = element_text(size = 13),
    axis.text = element_text(size = 13)
  )



##: Per drug combination: SGLT2vsTZD
drug1 = "SGLT2"; drug2 = "TZD"

title <- ATE.var_adj.no_adj_2nd_line[1,]
title[,] <- NA

interim_1 <- title %>% as.data.frame() %>%mutate(group = paste0("Predicted benefit on ", drug1, ":"), label = paste0("Predicted benefit on ", drug1, ":")) %>%
  rbind(
    ATE.var_adj.no_adj_2nd_line %>%
      filter(drug_1 == drug1 & drug_2 == drug2) %>%
      slice(which(grepl("(-", ATE.var_adj.no_adj_2nd_line$group, fixed =TRUE) == TRUE)),
    title %>% as.data.frame() %>% mutate(group = paste0("Predicted benefit on ", drug2, ":"), label = paste0("Predicted benefit on ", drug2, ":")),
    ATE.var_adj.no_adj_2nd_line %>%
      filter(drug_1 == drug1 & drug_2 == drug2) %>%
      slice(which(!grepl("(-", ATE.var_adj.no_adj_2nd_line$group, fixed =TRUE) == TRUE))
    
  )

title <- odds_ratios.no_weight.2nd_line[1,]
title[,] <- NA


interim_2 <- title %>% as.data.frame() %>%mutate(group = paste0("Predicted benefit on ", drug1, ":"), label = paste0("Predicted benefit on ", drug1, ":")) %>%
  rbind(
    odds_ratios.no_weight.2nd_line %>%
      filter(drug_1 == drug1 & drug_2 == drug2) %>%
      slice(which(grepl("(-", odds_ratios.no_weight.2nd_line$group, fixed =TRUE) == TRUE)),
    title %>% as.data.frame() %>% mutate(group = paste0("Predicted benefit on ", drug2, ":"), label = paste0("Predicted benefit on ", drug2, ":")),
    odds_ratios.no_weight.2nd_line %>%
      filter(drug_1 == drug1 & drug_2 == drug2) %>%
      slice(which(!grepl("(-", odds_ratios.no_weight.2nd_line$group, fixed =TRUE) == TRUE))
    
  )

title <- absolute.var_adj.no_adj_2n_line[1,]
title[,] <- NA

interim_3 <- title %>% as.data.frame() %>% mutate(group = paste0("Predicted benefit on ", drug1, ":"), label = paste0("Predicted benefit on ", drug1, ":")) %>%
  mutate(drug_type = drug1) %>%
  rbind(
    absolute.var_adj.no_adj_2n_line %>%
      filter(drug_1 == drug1 & drug_2 == drug2) %>%
      slice(which(grepl("(-", absolute.var_adj.no_adj_2n_line$group, fixed =TRUE) == TRUE)),
    title %>% as.data.frame() %>% mutate(group = paste0("Predicted benefit on ", drug2, ":"), label = paste0("Predicted benefit on ", drug2, ":")) %>%
      mutate(drug_type = drug2),
    absolute.var_adj.no_adj_2n_line %>%
      filter(drug_1 == drug1 & drug_2 == drug2) %>%
      slice(which(!grepl("(-", absolute.var_adj.no_adj_2n_line$group, fixed =TRUE) == TRUE))
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
  
) + plot_layout(axes = "collect") +
  plot_annotation(
    title = paste0(drug1, " vs ", drug2)
  ) &
  theme(
    axis.title.y = element_blank(),
    axis.title.x = element_text(size = 13),
    axis.text = element_text(size = 13)
  )




##: Per drug combination: SUvsTZD
drug1 = "SU"; drug2 = "TZD"

title <- ATE.var_adj.no_adj_2nd_line[1,]
title[,] <- NA

interim_1 <- title %>% as.data.frame() %>%mutate(group = paste0("Predicted benefit on ", drug1, ":"), label = paste0("Predicted benefit on ", drug1, ":")) %>%
  rbind(
    ATE.var_adj.no_adj_2nd_line %>%
      filter(drug_1 == drug1 & drug_2 == drug2) %>%
      slice(which(grepl("(-", ATE.var_adj.no_adj_2nd_line$group, fixed =TRUE) == TRUE)),
    title %>% as.data.frame() %>% mutate(group = paste0("Predicted benefit on ", drug2, ":"), label = paste0("Predicted benefit on ", drug2, ":")),
    ATE.var_adj.no_adj_2nd_line %>%
      filter(drug_1 == drug1 & drug_2 == drug2) %>%
      slice(which(!grepl("(-", ATE.var_adj.no_adj_2nd_line$group, fixed =TRUE) == TRUE))
    
  )

title <- odds_ratios.no_weight.2nd_line[1,]
title[,] <- NA


interim_2 <- title %>% as.data.frame() %>%mutate(group = paste0("Predicted benefit on ", drug1, ":"), label = paste0("Predicted benefit on ", drug1, ":")) %>%
  rbind(
    odds_ratios.no_weight.2nd_line %>%
      filter(drug_1 == drug1 & drug_2 == drug2) %>%
      slice(which(grepl("(-", odds_ratios.no_weight.2nd_line$group, fixed =TRUE) == TRUE)),
    title %>% as.data.frame() %>% mutate(group = paste0("Predicted benefit on ", drug2, ":"), label = paste0("Predicted benefit on ", drug2, ":")),
    odds_ratios.no_weight.2nd_line %>%
      filter(drug_1 == drug1 & drug_2 == drug2) %>%
      slice(which(!grepl("(-", odds_ratios.no_weight.2nd_line$group, fixed =TRUE) == TRUE))
    
  )

title <- absolute.var_adj.no_adj_2n_line[1,]
title[,] <- NA

interim_3 <- title %>% as.data.frame() %>% mutate(group = paste0("Predicted benefit on ", drug1, ":"), label = paste0("Predicted benefit on ", drug1, ":")) %>%
  mutate(drug_type = drug1) %>%
  rbind(
    absolute.var_adj.no_adj_2n_line %>%
      filter(drug_1 == drug1 & drug_2 == drug2) %>%
      slice(which(grepl("(-", absolute.var_adj.no_adj_2n_line$group, fixed =TRUE) == TRUE)),
    title %>% as.data.frame() %>% mutate(group = paste0("Predicted benefit on ", drug2, ":"), label = paste0("Predicted benefit on ", drug2, ":")) %>%
      mutate(drug_type = drug2),
    absolute.var_adj.no_adj_2n_line %>%
      filter(drug_1 == drug1 & drug_2 == drug2) %>%
      slice(which(!grepl("(-", absolute.var_adj.no_adj_2n_line$group, fixed =TRUE) == TRUE))
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
  
) + plot_layout(axes = "collect") +
  plot_annotation(
    title = paste0(drug1, " vs ", drug2)
  ) &
  theme(
    axis.title.y = element_blank(),
    axis.title.x = element_text(size = 13),
    axis.text = element_text(size = 13)
  )



dev.off()
