####################
## Description: 
##  - In this file we:
##    - Fit a linear PS model to all variables.
#################### 



# load working directory
setwd("Samples/T2D_Discontinuation")

# load functions
source("code/00.set_up_data_and_functions.R")

# load libraries
library(tidyverse)
library(PSweight)

###############################################################################
###############################################################################
############################### Read Data In ##################################
###############################################################################
###############################################################################

# load dataset
cprd_dataset <- set_up_data(
  raw_data = "20240308_t2d_1stinstance",
  diagnosis = FALSE,
  therapies = c("DPP4", "GLP1", "SGLT2", "SU", "TZD"),
  dataset = "3m.disc.dataset"
) %>%
  drop_na(-stopdrug_6m_6mFU, -stopdrug_12m_6mFU)


########
### Features for a propensity score
########

drug.pscores <- SumStat(
  ps.formula = as.formula(paste0("drugclass ~ ", paste0(c(
    # Extra info
    "dstartdate_dm_dur", "dstartdate_age", "drugline", "numdrugs", 
    "smoking_cat", "imd2015_10", "gender",
    # Biomarkers
    "prehba1c", "preegfr", "prebmi", "prealt", 
    "prehdl"
  ), collapse = "+"))),
  data = cprd_dataset %>% mutate(drugclass = factor(drugclass, levels = c("GLP1", "DPP4", "SGLT2", "SU", "TZD"))),
  weight = c("overlap", "IPW")
)

# summary(drug.pscores, metric = "ASD")

saveRDS(drug.pscores, "results/PS_model/drug.pscores_model.rds")


##########################################
#  Save Propensity scores
#

ps.only_dataset <- data.frame(
  patid = cprd_dataset$patid,
  dstartdate = cprd_dataset$dstartdate,
  prop.score.GLP1 = drug.pscores$propensity[,1],
  prop.score.DPP4 = drug.pscores$propensity[,2],
  prop.score.SGLT2 = drug.pscores$propensity[,3],
  prop.score.SU = drug.pscores$propensity[,4],
  prop.score.TZD = drug.pscores$propensity[,5],
  prop.score.overlap = drug.pscores$`ps.weights` %>% as.data.frame() %>% select(overlap) %>% unlist(),
  prop.score.IPW = drug.pscores$`ps.weights` %>% as.data.frame() %>% select(IPW) %>% unlist()
)

saveRDS(ps.only_dataset, "results/PS_model/ps.dataset_lm_all.rds")



##########################################
#  Check balance
#

# drug.pscores <- readRDS("results/PS_model/drug.pscores_model.rds")

pdf("results/figures/02.covariate_balance.pdf", width = 12, height = 10)
plot(drug.pscores)
dev.off()


