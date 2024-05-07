####################
## Description: 
##  - In this file we:
##    - investigate discontinuation of GLP1RA
#################### 



# load working directory
setwd("Samples/T2D_Discontinuation")

# load functions
source("code/00.set_up_data.R")

# load libraries
library(tidyverse)
library(pROC)
library(caret)

###############################################################################
###############################################################################
############################### Read Data In ##################################
###############################################################################
###############################################################################

# load dataset
cprd_dataset <- set_up_data(
  raw_data = "20240308_t2d_1stinstance",
  diagnosis = FALSE,
  therapies = c("GLP1", "SGLT2"),
  dataset = "full.dataset"
) %>%
  drop_na(stopdrug_3m_6mFU)




cprd_dataset_model <- cprd_dataset %>%
  select(
    patid,
    dstartdate,
    stopdrug_3m_6mFU,
    drugclass,
    gender,
    dstartdate_age,
    prebmi,
    prehba1c,
    preegfr,
    dstartdate_dm_dur,
    prealt,
    prehdl,
    pretotalcholesterol,
    # preldl,
    ethnicity_5cat,
    smoking_cat,
    imd2015_10,
    drugline,
    numdrugs,
    stopdrug_3m_3mFU_MFN_hist,
    predrug_statins,
    CCI_index
  ) %>%
  mutate(
    drugclass = factor(drugclass, levels = c("GLP1", "SGLT2")),
    gender = factor(gender, levels = c(2, 1), labels = c("Female", "Male"))
  ) %>%
  # rename(
  #   "sex" = "gender",
  #   "agetx" = "dstartdate_age",
  #   "t2dmduration" = "dstartdate_dm_dur",
  #   "ncurrtx" = "numdrugs"
  # ) %>%
  as.data.frame() %>%
  drop_na() %>%
  filter(drugline != "1") %>%
  mutate(
    drugline = factor(drugline)
  )




##################################################################################
##################################################################################


## model 1 - can we predict discontinuation 

model_1 <- glm(stopdrug_3m_6mFU ~ gender + dstartdate_age + prebmi + prehba1c + preegfr + dstartdate_dm_dur + prealt + prehdl + pretotalcholesterol + ethnicity_5cat + smoking_cat + imd2015_10 + drugline + numdrugs + stopdrug_3m_3mFU_MFN_hist + CCI_index,
             family=binomial(link='logit'),
             data=cprd_dataset_model %>%
               filter(drugclass == "GLP1"))


## roc curve
pROC::roc(model_1$model$stopdrug_3m_6mFU, model_1$fitted.values)


## var importance
V = caret::varImp(model_1)

ggplot2::ggplot(V, aes(x=reorder(rownames(V),Overall), y=Overall)) +
  geom_point( color="blue", size=4, alpha=0.6)+
  geom_segment( aes(x=rownames(V), xend=rownames(V), y=0, yend=Overall), 
                color='skyblue') +
  xlab('Variable')+
  ylab('Overall Importance')+
  theme_light() +
  coord_flip() 




## model 2 - can we predict differential discontinuation on GLP1vsSGLT2

model_2 <- glm(stopdrug_3m_6mFU ~ gender*drugclass + 
                 dstartdate_age*drugclass + 
                 prebmi*drugclass + 
                 prehba1c*drugclass + 
                 preegfr*drugclass + 
                 dstartdate_dm_dur*drugclass + 
                 prealt*drugclass + 
                 prehdl*drugclass + 
                 pretotalcholesterol*drugclass + 
                 stopdrug_3m_3mFU_MFN_hist*drugclass +
                 CCI_index*drugclass +
                 ethnicity_5cat + smoking_cat + imd2015_10 + drugline + numdrugs,
               family=binomial(link='logit'),
               data=cprd_dataset_model)

## roc curve
pROC::roc(model_2$model$stopdrug_3m_6mFU, model_2$fitted.values)


## var importance
V = caret::varImp(model_2)

ggplot2::ggplot(V, aes(x=reorder(rownames(V),Overall), y=Overall)) +
  geom_point( color="blue", size=4, alpha=0.6)+
  geom_segment( aes(x=rownames(V), xend=rownames(V), y=0, yend=Overall), 
                color='skyblue') +
  xlab('Variable')+
  ylab('Overall Importance')+
  theme_light() +
  coord_flip() 





dev_data$pred <- predict.glm(model_2, type = 'response')
dev_data$y <- as.numeric(cprd_dataset_model$stopdrug_3m_6mFU) -1
dev_data$group <- cprd_dataset_model$stopdrug_3m_3mFU_MFN_hist
dev_data <- dev_data %>% as.data.frame()

calibration_plot(data = dev_data, obs = "y", pred = "pred", title = "Calibration plot for development data", y_lim = c(0, 0.4), x_lim=c(0, 0.4), nTiles = 30, group = "group")
#> $calibration_plot

## abstract about basic features do not predict discontinuation in observational study







######################################################################################
######################################################################################

# Is the tolerability of GLP1-RA modified by routine clinical features (age, duration, sex, BMI, HbA1c) and comorbidity profile?” [with SGLT2 as comparator group]



cprd_dataset_model <- cprd_dataset %>%
  select(
    patid,
    dstartdate,
    stopdrug_3m_6mFU,
    drugclass,
    gender,
    dstartdate_age,
    prebmi,
    prehba1c,
    dstartdate_dm_dur,
    smoking_cat,
    imd2015_10,
    drugline,
    numdrugs,
    stopdrug_3m_3mFU_MFN_hist,
    predrug_statins
  ) %>%
  mutate(
    drugclass = factor(drugclass, levels = c("GLP1", "SGLT2")),
    gender = factor(gender, levels = c(2, 1), labels = c("Female", "Male"))
  ) %>%
  # rename(
  #   "sex" = "gender",
  #   "agetx" = "dstartdate_age",
  #   "t2dmduration" = "dstartdate_dm_dur",
  #   "ncurrtx" = "numdrugs"
  # ) %>%
  as.data.frame() %>%
  drop_na() %>%
  filter(drugline != "1") %>%
  mutate(
    drugline = factor(drugline)
  )  %>%
  mutate(
    prehba1c = prehba1c/10,
    dstartdate_age = dstartdate_age/10,
    dstartdate_dm_dur = dstartdate_dm_dur/10,
    prebmi = prebmi/5
  )





## model 1 - can we predict discontinuation 

model_1 <- glm(stopdrug_3m_6mFU ~ gender + dstartdate_age + dstartdate_dm_dur + prehba1c + prebmi + imd2015_10 + smoking_cat + drugline + numdrugs + stopdrug_3m_3mFU_MFN_hist + predrug_statins,
               family=binomial(link='logit'),
               data=cprd_dataset_model %>%
                 filter(drugclass == "GLP1"))



## model 2 - can we predict discontinuation 

model_2 <- glm(stopdrug_3m_6mFU ~ gender + dstartdate_age + dstartdate_dm_dur + prehba1c + prebmi + imd2015_10 + smoking_cat + drugline + numdrugs + stopdrug_3m_3mFU_MFN_hist + predrug_statins,
               family=binomial(link='logit'),
               data=cprd_dataset_model %>%
                 filter(drugclass == "SGLT2"))



## model 3 - can we predict discontinuation 

model_3 <- glm(stopdrug_3m_6mFU ~ gender*drugclass + dstartdate_age*drugclass + dstartdate_dm_dur*drugclass + prehba1c*drugclass + prebmi*drugclass + imd2015_10 + smoking_cat + drugline + numdrugs + stopdrug_3m_3mFU_MFN_hist + predrug_statins,
               family=binomial(link='logit'),
               data=cprd_dataset_model)



# cardio, renal and liver

# probability of subgroup a is pred on drug a vs drug b

cprd_dataset_model <- cprd_dataset_model %>%
  mutate(subgroup = ifelse(prebmi >=35/5 & prehba1c <75/10, "a", ifelse(prebmi <35/5 & prehba1c >=75/10, "b", NA)))


model_4 <- glm(stopdrug_3m_6mFU ~ gender*drugclass + dstartdate_age*drugclass + dstartdate_dm_dur*drugclass + prehba1c*drugclass + prebmi*drugclass + imd2015_10 + smoking_cat + drugline + numdrugs + stopdrug_3m_3mFU_MFN_hist + predrug_statins,
               family=binomial(link='logit'),
               data=cprd_dataset_model)


subgroup_a <- cprd_dataset_model %>%
  filter(subgroup == "a") %>%
  mutate(drugclass = "SGLT2")

mean(predict(model_3, subgroup_a, type = "response"))

subgroup_a <- cprd_dataset_model %>%
  filter(subgroup == "a") %>%
  mutate(drugclass = "GLP1")

mean(predict(model_3, subgroup_a, type = "response"))


subgroup_b <- cprd_dataset_model %>%
  filter(subgroup == "b") %>%
  mutate(drugclass = "SGLT2")

mean(predict(model_3, subgroup_b, type = "response"))

subgroup_b <- cprd_dataset_model %>%
  filter(subgroup == "b") %>%
  mutate(drugclass = "GLP1")

mean(predict(model_3, subgroup_b, type = "response"))






















