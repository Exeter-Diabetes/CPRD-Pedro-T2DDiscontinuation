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






pdf("results/figures/08.absolute_relative_development.pdf", width = 20, height = 7)


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
  
) + plot_layout(axes = "collect") +
  plot_annotation(
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
  
) + plot_layout(axes = "collect") +
  plot_annotation(
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
  
) + plot_layout(axes = "collect") +
  plot_annotation(
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
  
) + plot_layout(axes = "collect") +
  plot_annotation(
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
  
) + plot_layout(axes = "collect") +
  plot_annotation(
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




