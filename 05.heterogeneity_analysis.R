####################
## Description: 
##  - In this file we:
##    - Investigate heterogeneity of therapy discontinuation
#################### 


# load working directory
setwd("Samples/T2D_Discontinuation")

# load functions
source("code/00.set_up_data.R")

# load libraries
library(tidyverse)
library(MatchIt)
library(qpdf)

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



# load propensity scores
ps.only_dataset <- readRDS("results/PS_model/ps.dataset_lm_all.rds")


# join propensity scores into dataset
cprd_dataset.dev <- cprd_dataset.dev %>%
  left_join(
    ps.only_dataset %>%
      select(patid, dstartdate, prop.score.overlap, prop.score.IPW),
    by = c("patid", "dstartdate")
  )
cprd_dataset.val <- cprd_dataset.val %>%
  left_join(
    ps.only_dataset %>%
      select(patid, dstartdate, prop.score.overlap, prop.score.IPW),
    by = c("patid", "dstartdate")
  )


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


calc_ATE <- function(data, drugs, pred.variable = "no_weight", weight.variable = NULL, breakdown = NULL, matching = FALSE, ntiles = 10, caliper = 0.05, replace = FALSE, order = "random", n_bootstrap = 100) {
  
  # output object
  final_dataset <- NULL
  
  # calculation of all drug combinations (from provided)
  combinations_matrix <- combn(drugs, 2)
  
  # number of combinations
  n_comb <- ncol(combinations_matrix)
  
  # iterate through each combinations
  for (comb_i in 1:n_comb) {
    
    # current drugs being analysed
    current_drugs <- combinations_matrix[,comb_i]
    
    # collect only the drugs we are interested
    data.initial <- data %>%
      filter(drugclass %in% current_drugs) %>%
      rename(
        "pred.drug_1" = paste0("pred.", pred.variable, ".", current_drugs[1]),
        "pred.drug_2" = paste0("pred.", pred.variable, ".", current_drugs[2]),
      ) %>%
      mutate(diff = pred.drug_1 - pred.drug_2) %>%
      mutate(group = ntile(diff, ntiles))
    
    # if weight variable is provided, rename propensity score variable to default
    if (!is.null(weight.variable)) {
      data.initial <- data.initial %>%
        rename(
          "propensity.score" = paste0("prop.score.", weight.variable)
        )
    }
    
    # calculate the summary statistics for each ntile
    predicted_treatment_effect <- data.initial %>%
      plyr::ddply("group", dplyr::summarise,
                  N = length(diff),
                  diff.pred = mean(diff)) %>%
      left_join(
        data.initial %>%
          group_by(group) %>%
          mutate(
            events = sum(as.numeric(stopdrug_3m_6mFU) - 1, na.rm = TRUE)
          ) %>%
          ungroup() %>%
          select(group, events) %>%
          unique() %>%
          mutate(
            drug_1 = current_drugs[1],
            drug_2 = current_drugs[2]
          ),
        by = c("group")
      )
    
    
    # maximum number of deciles being tested
    quantiles <- length(unique(data.initial[,"group"]))
    
    # create lists with results
    mnumber = c(1:quantiles)
    models  <- as.list(1:quantiles)
    obs <- vector(); lci <- vector(); uci <- vector()
    
    # iterate through deciles
    for (i in mnumber) {
      
      # do this differently if "group" is categorical
      if (is.factor(data.initial[,"group"])) {
        # dataset being used in this quantile
        data.new <- data.initial[data.initial[,"group"] == levels(data.initial[,"group"])[i],]
      } else {
        # dataset being used in this quantile
        data.new <- data.initial[data.initial[,"group"] == i,]
      }
      
      # If breakdown is not given (breakdown = variable to adjust by)
      if (is.null(breakdown)) {
        
        # formula
        formula <- "stopdrug_3m_6mFU ~ factor(drugclass)"
        
      } else {
        
        # variables used in adjustment
        breakdown_adjust <- breakdown
        # variables with only one variable represented
        checker <- which(sapply(data.new[,breakdown_adjust], function(col) length(unique(col))) > 1)
        
        # formula
        formula <- paste0("stopdrug_3m_6mFU ~ factor(drugclass) +", paste(breakdown_adjust[checker], collapse = "+"))
        
      }
      
      # if (!is.null(weight.variable)) {
      #   
      #   if (matching == FALSE) {
      #     
      #     # fit linear regression for decile in the matched dataset
      #     models[[i]] <- glm(as.formula("stopdrug_3m_6mFU ~ factor(drugclass)"), data=data.new, weights = propensity.score, family = quasibinomial())
      #     
      #   } else {
      #     
      #     # model if propensity scores are provided
      #     matching_package_result <- MatchIt::matchit(
      #       formula = formula("drugclass ~ prehba1c"), # shouldn't be used since we are specifying 'distance' (propensity scores)
      #       data = data.new[which(data.new[,"group"] == i),], # select people in the quantile
      #       method = "nearest",
      #       distance = data.new[which(data.new[,"group"] == i),"propensity.score"],
      #       replace = replace,
      #       m.order = order,
      #       caliper = caliper,
      #       mahvars = NULL, estimand = "ATT", exact = NULL, antiexact = NULL, discard = "none", reestimate = FALSE, s.weights = NULL, std.caliper = TRUE, ratio = 1, verbose = FALSE, include.obj = FALSE,
      #     )
      #     
      #     # fit linear regression for decile in the matched dataset
      #     models[[i]] <- glm(as.formula("stopdrug_3m_6mFU ~ factor(drugclass)"), data=data.new, weights = matching_package_result$weights, family = binomial())
      #     
      #   }
      #   
      # } else {
      #   
      #   # fit linear regression for decile in the matched dataset
      #   models[[i]] <- glm(as.formula(formula),data=data.new, family = binomial())
      #   
      # }
      
      # iterate through bootstrapped datasets to get CI
      
      ## WRONG: bootstrap before model being run
      
      bootstrap_differences <- NULL
      
      for (bootstrap_i in 1:n_bootstrap) {
        
        # bootstrap dataset
        data.bootstrap <- data.new[sample(nrow(data.new), size = nrow(data.new), replace = TRUE),] 
        
        
        if (!is.null(weight.variable)) {
          
          if (matching == FALSE) {
            
            # fit linear regression for decile in the matched dataset
            models[[i]] <- glm(as.formula("stopdrug_3m_6mFU ~ factor(drugclass)"), data=data.bootstrap, weights = propensity.score, family = quasibinomial())
            
          } else {
            
            # model if propensity scores are provided
            matching_package_result <- MatchIt::matchit(
              formula = formula("drugclass ~ prehba1c"), # shouldn't be used since we are specifying 'distance' (propensity scores)
              data = data.bootstrap[which(data.bootstrap[,"group"] == i),], # select people in the quantile
              method = "nearest",
              distance = data.bootstrap[which(data.bootstrap[,"group"] == i),"propensity.score"],
              replace = replace,
              m.order = order,
              caliper = caliper,
              mahvars = NULL, estimand = "ATT", exact = NULL, antiexact = NULL, discard = "none", reestimate = FALSE, s.weights = NULL, std.caliper = TRUE, ratio = 1, verbose = FALSE, include.obj = FALSE,
            )
            
            # fit linear regression for decile in the matched dataset
            models[[i]] <- glm(as.formula("stopdrug_3m_6mFU ~ factor(drugclass)"), data=data.bootstrap, weights = matching_package_result$weights, family = binomial())
            
          }
          
        } else {
          
          # fit linear regression for decile in the matched dataset
          models[[i]] <- glm(as.formula(formula),data=data.bootstrap, family = binomial())
          
        }
        
        # make predictions
        preds.drug1 <- predict(models[[i]], newdata = data.bootstrap %>% filter(drugclass == current_drugs[1]), type = "response")
        
        preds.drug2 <- predict(models[[i]], newdata = data.bootstrap %>% filter(drugclass == current_drugs[2]), type = "response")
        
        # save differences
        bootstrap_differences <- append(bootstrap_differences, mean(preds.drug1) - mean(preds.drug2))
        
      }
      
      # collect mean
      obs <- append(obs, quantile(bootstrap_differences, probs = c(0.5)))
      
      # collect lower bound CI
      lci <- append(lci, quantile(bootstrap_differences, probs = c(0.025)))
      
      # collect upper bound CI
      uci <- append(uci, quantile(bootstrap_differences, probs = c(0.975)))
      
    }
    
    # join treatment effects for deciles in a data.frame
    effects <- data.frame(predicted_treatment_effect,cbind(obs, lci, uci))
    
    # returned list with fitted propensity model + decile treatment effects
    final_dataset <- rbind(
      final_dataset,
      effects %>% mutate(Type = paste(current_drugs[1], "vs", current_drugs[2]))
    )
    
  }
  
  return(final_dataset)
  
}



#
#:------------ All drugs model
#

# variables to adjust
breakdown <- c(
  # Extra info
  "dstartdate_dm_dur", "dstartdate_age", "drugline", "numdrugs",
  "smoking_cat", "imd2015_10", "gender",
  # Biomarkers
  "prehba1c", "preegfr", "prebmi", "prealt",
  "pretotalcholesterol"
)


#:--------------
## UBER model variable adjusted

# ATE not adjusted

ATE.var_adj.no_adj_all_drugs <- calc_ATE(cprd_dataset.dev, drugs = c("MFN", "DPP4", "GLP1", "SGLT2", "SU", "TZD"), pred.variable = "no_weight", n_bootstrap = 50)

saveRDS(ATE.var_adj.no_adj_all_drugs, "results/Heterogeneity/ATE.var_adj.no_adj_all_drugs.rds")

# ATE overlap matching

ATE.var_adj.overlap_match_all_drugs <- calc_ATE(cprd_dataset.dev, drugs = c("MFN", "DPP4", "GLP1", "SGLT2", "SU", "TZD"), pred.variable = "no_weight", weight.variable = "overlap", matching = TRUE, n_bootstrap = 50)

saveRDS(ATE.var_adj.overlap_match_all_drugs, "results/Heterogeneity/ATE.var_adj.overlap_match_all_drugs.rds")

# ATE IPW matching

ATE.var_adj.IPW_match_all_drugs <- calc_ATE(cprd_dataset.dev, drugs = c("MFN", "DPP4", "GLP1", "SGLT2", "SU", "TZD"), pred.variable = "no_weight", weight.variable = "IPW", matching = TRUE, n_bootstrap = 50)

saveRDS(ATE.var_adj.IPW_match_all_drugs, "results/Heterogeneity/ATE.var_adj.IPW_match_all_drugs.rds")







pdf("results/figures/plot_1.pdf", width = 10, height = 7)

ATE.var_adj.no_adj_all_drugs %>%
  filter(drug_1 == "MFN" | drug_2 == "MFN") %>%
  ggplot() +
  geom_abline(aes(intercept = 0, slope = 1), colour = "red") +
  geom_hline(aes(yintercept = 0), colour = "grey") +
  geom_vline(aes(xintercept = 0), colour = "grey") +
  geom_point(aes(y = obs, x = diff.pred)) +
  geom_errorbar(aes(y = obs, ymin = lci, ymax = uci, x = diff.pred)) +
  scale_x_continuous("Predicted Conditional Average Treatment Discontinuation Difference", labels = scales::percent, breaks = seq(-1, 1, 0.02)) +
  scale_y_continuous("Decile Average Treatment Discontinuation Difference", labels = scales::percent, breaks = seq(-1, 1, 0.02)) +
  ggtitle("Comparison of MFN against other therapies") +
  facet_wrap(~Type, scales = "free", nrow = 2) +
  theme_bw()

dev.off()



pdf("results/figures/plot_2.pdf", width = 12, height = 9)

ATE.var_adj.no_adj_all_drugs %>%
  filter(drug_1 != "MFN" & drug_2 != "MFN") %>%
  ggplot() +
  geom_abline(aes(intercept = 0, slope = 1), colour = "red") +
  geom_hline(aes(yintercept = 0), colour = "grey") +
  geom_vline(aes(xintercept = 0), colour = "grey") +
  geom_point(aes(y = obs, x = diff.pred)) +
  geom_errorbar(aes(y = obs, ymin = lci, ymax = uci, x = diff.pred)) +
  scale_x_continuous("Predicted Conditional Average Treatment Discontinuation Difference", labels = scales::percent, breaks = seq(-1, 1, 0.02)) +
  scale_y_continuous("Decile Average Treatment Discontinuation Difference (no adjustment)", labels = scales::percent, breaks = seq(-1, 1, 0.02)) +
  ggtitle("Therapy discontinuation heterogeneity (no adjustment)") +
  facet_wrap(~Type, scales = "free") +
  theme_bw()

ATE.var_adj.overlap_match_all_drugs %>%
  filter(drug_1 != "MFN" & drug_2 != "MFN") %>%
  ggplot() +
  geom_abline(aes(intercept = 0, slope = 1), colour = "red") +
  geom_hline(aes(yintercept = 0), colour = "grey") +
  geom_vline(aes(xintercept = 0), colour = "grey") +
  geom_point(aes(y = obs, x = diff.pred)) +
  geom_errorbar(aes(y = obs, ymin = lci, ymax = uci, x = diff.pred)) +
  scale_x_continuous("Predicted Conditional Average Treatment Discontinuation Difference", labels = scales::percent, breaks = seq(-1, 1, 0.02)) +
  scale_y_continuous("Decile Average Treatment Discontinuation Difference", labels = scales::percent, breaks = seq(-1, 1, 0.02)) +
  ggtitle("Therapy discontinuation heterogeneity (Overlap weight matching)") +
  facet_wrap(~Type, scales = "free") +
  theme_bw()

ATE.var_adj.IPW_match_all_drugs %>%
  filter(drug_1 != "MFN" & drug_2 != "MFN") %>%
  ggplot() +
  geom_abline(aes(intercept = 0, slope = 1), colour = "red") +
  geom_hline(aes(yintercept = 0), colour = "grey") +
  geom_vline(aes(xintercept = 0), colour = "grey") +
  geom_point(aes(y = obs, x = diff.pred)) +
  geom_errorbar(aes(y = obs, ymin = lci, ymax = uci, x = diff.pred)) +
  scale_x_continuous("Predicted Conditional Average Treatment Discontinuation Difference", labels = scales::percent, breaks = seq(-1, 1, 0.02)) +
  scale_y_continuous("Decile Average Treatment Discontinuation Difference", labels = scales::percent, breaks = seq(-1, 1, 0.02)) +
  ggtitle("Therapy discontinuation heterogeneity (Overlap weight matching)") +
  facet_wrap(~Type, scales = "free") +
  theme_bw()

dev.off()












#
#:------------ 2nd line drugs model
#

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



# load propensity scores
ps.only_dataset <- readRDS("results/PS_model/ps.dataset_lm_all.rds")


# join propensity scores into dataset
cprd_dataset.dev <- cprd_dataset.dev %>%
  left_join(
    ps.only_dataset %>%
      select(patid, dstartdate, prop.score.overlap, prop.score.IPW),
    by = c("patid", "dstartdate")
  )
cprd_dataset.val <- cprd_dataset.val %>%
  left_join(
    ps.only_dataset %>%
      select(patid, dstartdate, prop.score.overlap, prop.score.IPW),
    by = c("patid", "dstartdate")
  )


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



# variables to adjust
breakdown <- c(
  # Extra info
  "dstartdate_dm_dur", "dstartdate_age", "drugline", "numdrugs",
  "smoking_cat", "imd2015_10", "gender",
  # Biomarkers
  "prehba1c", "preegfr", "prebmi", "prealt",
  "pretotalcholesterol"
)


#:--------------
## UBER model variable adjusted

# ATE not adjusted

ATE.var_adj.no_adj_2nd_line <- calc_ATE(cprd_dataset.dev, drugs = c("DPP4", "GLP1", "SGLT2", "SU", "TZD"), pred.variable = "no_weight", n_bootstrap = 50)

saveRDS(ATE.var_adj.no_adj_2nd_line, "results/Heterogeneity/ATE.var_adj.no_adj_2nd_line.rds")

# ATE overlap matching

ATE.var_adj.overlap_match_2nd_line <- calc_ATE(cprd_dataset.dev, drugs = c("DPP4", "GLP1", "SGLT2", "SU", "TZD"), pred.variable = "no_weight", weight.variable = "overlap", matching = TRUE, n_bootstrap = 50)

saveRDS(ATE.var_adj.overlap_match_2nd_line, "results/Heterogeneity/ATE.var_adj.overlap_match_2nd_line.rds")

# ATE IPW matching

ATE.var_adj.IPW_match_2nd_line <- calc_ATE(cprd_dataset.dev, drugs = c("DPP4", "GLP1", "SGLT2", "SU", "TZD"), pred.variable = "no_weight", weight.variable = "IPW", matching = TRUE, n_bootstrap = 50)

saveRDS(ATE.var_adj.IPW_match_2nd_line, "results/Heterogeneity/ATE.var_adj.IPW_match_2nd_line.rds")





pdf("results/figures/plot_3.pdf", width = 12, height = 9)

ATE.var_adj.no_adj_2nd_line %>%
  ggplot() +
  geom_abline(aes(intercept = 0, slope = 1), colour = "red") +
  geom_hline(aes(yintercept = 0), colour = "grey") +
  geom_vline(aes(xintercept = 0), colour = "grey") +
  geom_point(aes(y = obs, x = diff.pred)) +
  geom_errorbar(aes(y = obs, ymin = lci, ymax = uci, x = diff.pred)) +
  scale_x_continuous("Predicted Conditional Average Treatment Discontinuation Difference", labels = scales::percent, breaks = seq(-1, 1, 0.02)) +
  scale_y_continuous("Decile Average Treatment Discontinuation Difference (no adjustment)", labels = scales::percent, breaks = seq(-1, 1, 0.02)) +
  ggtitle("2nd line Therapy discontinuation heterogeneity (no adjustment)") +
  facet_wrap(~Type, scales = "free") +
  theme_bw()

ATE.var_adj.overlap_match_2nd_line %>%
  ggplot() +
  geom_abline(aes(intercept = 0, slope = 1), colour = "red") +
  geom_hline(aes(yintercept = 0), colour = "grey") +
  geom_vline(aes(xintercept = 0), colour = "grey") +
  geom_point(aes(y = obs, x = diff.pred)) +
  geom_errorbar(aes(y = obs, ymin = lci, ymax = uci, x = diff.pred)) +
  scale_x_continuous("Predicted Conditional Average Treatment Discontinuation Difference", labels = scales::percent, breaks = seq(-1, 1, 0.02)) +
  scale_y_continuous("Decile Average Treatment Discontinuation Difference", labels = scales::percent, breaks = seq(-1, 1, 0.02)) +
  ggtitle("2nd line Therapy discontinuation heterogeneity (Overlap weight matching)") +
  facet_wrap(~Type, scales = "free") +
  theme_bw()

ATE.var_adj.IPW_match_2nd_line %>%
  ggplot() +
  geom_abline(aes(intercept = 0, slope = 1), colour = "red") +
  geom_hline(aes(yintercept = 0), colour = "grey") +
  geom_vline(aes(xintercept = 0), colour = "grey") +
  geom_point(aes(y = obs, x = diff.pred)) +
  geom_errorbar(aes(y = obs, ymin = lci, ymax = uci, x = diff.pred)) +
  scale_x_continuous("Predicted Conditional Average Treatment Discontinuation Difference", labels = scales::percent, breaks = seq(-1, 1, 0.02)) +
  scale_y_continuous("Decile Average Treatment Discontinuation Difference", labels = scales::percent, breaks = seq(-1, 1, 0.02)) +
  ggtitle("2nd line Therapy discontinuation heterogeneity (Overlap weight matching)") +
  facet_wrap(~Type, scales = "free") +
  theme_bw()

dev.off()





#:-----------------------------------------------------------------
# Combining pdfs

qpdf::pdf_combine(input = c("results/figures/plot_1.pdf", "results/figures/plot_2.pdf", "results/figures/plot_3.pdf"),
                  output = "results/figures/heterogeneity_differences_analysis.pdf")


file.remove(c("results/figures/plot_1.pdf", "results/figures/plot_2.pdf", "results/figures/plot_3.pdf"))
