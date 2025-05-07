######################################################################
##
##  In this file, we create a function for extracting the dataset
##    required for analysis.
##
######################################################################



# Data setup function

set_up_data <- function(
    raw_data = "20240308_t2d_1stinstance",
    diagnosis = FALSE,
    therapies = c("DPP4", "GLP1", "MFN", "SGLT2", "SU", "TZD"),
    dataset = "full.dataset",
    follow_up = "6-months",
    full_prescribing_history = TRUE
) {
  
  # load functions
  require(tidyverse)
  
  # Checks for 'raw_data'
  if(missing(raw_data)) {stop("'raw_data' needs to be supplied.")}
  # Checks for 'diagnosis
  if(missing(diagnosis)) {stop("'diagnosis' needs to be supplied.")}
  if(!(diagnosis %in% c(TRUE, FALSE))) {stop("'diagnosis' must be TRUE or FALSE.")}
  # Check for 'dataset'
  if(missing(dataset)) {stop("'dataset' needs to be supplied.")}
  if(!(dataset %in% c("full.dataset",
                      "3m.disc.dataset", "3m.disc.dataset.dev", "3m.disc.dataset.val", 
                      "6m.disc.dataset", "6m.disc.dataset.dev", "6m.disc.dataset.val", 
                      "12m.disc.dataset", "12m.disc.dataset.dev", "12m.disc.dataset.val"
  ))) {
    stop("'dataset' needs to be: full.dataset / 3m.disc.dataset / 3m.disc.dataset.dev / 3m.disc.dataset.val / 6m.disc.dataset / 6m.disc.dataset.dev / 6m.disc.dataset.val / 12m.disc.dataset / 12m.disc.dataset.dev / 12m.disc.dataset.val")
  }
  # Check for 'follow_up'
  if(missing(follow_up)) {stop("'follow_up' needs to be supplied.")}
  if(!(follow_up %in% c("3-months", "6-months"))) {stop("'follow_up' needs to be: 3-months / 6-months")}
  # Check for 'full_prescribing_history'
  if(missing(full_prescribing_history)) {stop("'full_prescribing_history' needs to be supplied.")}
  if(!(full_prescribing_history %in% c(TRUE, FALSE))) {stop("'full_prescribing_history' must be TRUE or FALSE.")}
  
  
  ###############################################
  
  
  # load dataset - 1st instance
  load(paste0("/slade/CPRD_data/2020_dataset/Mastermind/", raw_data, ".Rda"))   # name: t2d_1stinstance
  
  cprd_dataset <- t2d_1stinstance
  
  ###############################################
  
  # Select therapies of interest
  
  cprd_dataset <- cprd_dataset %>%
    filter(drugclass %in% therapies)
  
  
  ## Check therapies being initiated
  if (isTRUE(diagnosis)) {
    print("#####################################")
    print(paste("Patients initiating therapies:", nrow(cprd_dataset)))
    print(table(cprd_dataset$drugclass))
    print("#####################################")
  }
  
  
  ###############################################
  
  # If TRUE, remove missing diagnosis date
  if (isTRUE(full_prescribing_history)) {
    
    # Keep patients with full prescribing data
    cprd_dataset <- cprd_dataset %>%
      filter(!is.na(dm_diag_date))
    
    
    ## Check therapies being initiated
    if (isTRUE(diagnosis)) {
      print("#####################################")
      print(paste("Patients with full prescribing data:", nrow(cprd_dataset)))
      print(table(cprd_dataset$drugclass))
      print("#####################################")
    }
    
  }
  
  
  ###############################################
  
  # Keep patients above 2014-01-01
  
  cprd_dataset <- cprd_dataset %>%
    mutate(dstartdate = as.Date(dstartdate)) %>%
    filter(dstartdate > "2014-01-01")
  
  
  ## Check patients after data
  if (isTRUE(diagnosis)) {
    print("#####################################")
    print(paste("Patients after 2014-01-01:", nrow(cprd_dataset)))
    print(table(cprd_dataset$drugclass, useNA = "ifany"))
    print("#####################################")
  }
  
  
  ###############################################
  
  # Remove patients under 18
  
  cprd_dataset <- cprd_dataset %>%
    filter(dstartdate_age >= 18)
  
  
  ## Check patients after data
  if (isTRUE(diagnosis)) {
    print("#####################################")
    print(paste("Only adults 2014-01-01:", nrow(cprd_dataset)))
    print("#####################################")
  }
  
  ###############################################
  
  # Remove patients prescribed other glucose-lowering therapies as 1st line
  if (!("MFN" %in% therapies)) {
    
    cprd_dataset <- cprd_dataset %>%
      filter(drugline != 1)
    
    ## Check patients after data
    if (isTRUE(diagnosis)) {
      print("#####################################")
      print(paste("Only 2nd-line therapies:", nrow(cprd_dataset)))
      print("#####################################")
    }
    
  }
  
  
  ###############################################
  
  # Remove patients starting multiple drugs
  
  cprd_dataset <- cprd_dataset %>%
    filter(multi_drug_start == 0)
  
  
  ## Check patients after data
  if (isTRUE(diagnosis)) {
    print("#####################################")
    print(paste("Patients initiating only one therapy:", nrow(cprd_dataset)))
    print(table(cprd_dataset$drugclass, useNA = "ifany"))
    print("#####################################")
  }
  
  ###############################################
  
  # Remove patients on insulin
  
  cprd_dataset <- cprd_dataset %>%
    filter(INS == 0)
  
  ## Check patients after data
  if (isTRUE(diagnosis)) {
    print("#####################################")
    print(paste("Patients not on INS (insulin):", nrow(cprd_dataset)))
    print(table(cprd_dataset$drugclass, useNA = "ifany"))
    print("#####################################")
  }
  
  #####################################################################################
  #####################################################################################
  
  # Remove patients without HbA1c
  
  cprd_dataset <- cprd_dataset %>%
    filter(!is.na(prehba1c))
  
  
  ## Check patients after data
  if (isTRUE(diagnosis)) {
    print("#####################################")
    print(paste("Patients with HbA1c:", nrow(cprd_dataset)))
    print("#####################################")
  }
  
  #####################################################################################
  #####################################################################################
  
  # Remove patients with HbA1c above 53
  
  cprd_dataset <- cprd_dataset %>%
    filter(prehba1c > 53)
  
  
  ## Check patients after data
  if (isTRUE(diagnosis)) {
    print("#####################################")
    print(paste("Patients HbA1c above 53:", nrow(cprd_dataset)))
    print("#####################################")
  }
  
  
  #####################################################################################
  #####################################################################################
  
  # Remove patients with 3 month discontinuation data
  
  if (follow_up %in% c("3-months")) {
    
    cprd_dataset <- cprd_dataset %>%
      drop_na(stopdrug_3m_3mFU) 
    
    
    ## Check patients after data
    if (isTRUE(diagnosis)) {
      print("#####################################")
      print(paste("Patients with 3 month discontinuation data (3-month follow-up):", nrow(cprd_dataset)))
      print("#####################################")
    }
    
  } else {
    
    cprd_dataset <- cprd_dataset %>%
      drop_na(stopdrug_3m_6mFU) 
    
    
    ## Check patients after data
    if (isTRUE(diagnosis)) {
      print("#####################################")
      print(paste("Patients with 3 month discontinuation data (6-month follow-up):", nrow(cprd_dataset)))
      print("#####################################")
    }
    
  }
  
  
  # Create necessary variables
  
  ### selecting variables: Missing - number of chronic illnesses q: just sum over the vars listed below?
  cprd_dataset <- cprd_dataset %>%
    mutate(
      
      # Check if patient only had one prescription
      only_one_prescription = ifelse(as.Date(dstartdate) == as.Date(dstopdate), 1, 0),
      
      # Statins use
      predrug_statins = ifelse(!is.na(predrug_latest_statins), 1, 0),
      
      # Blood pressure medication in the last 6 months
      predrug_bloodmed = ifelse(
        (!is.na(predrug_latest_ace_inhibitors) & as.numeric(difftime(dstartdate, predrug_latest_ace_inhibitors)/30) < 6) |
          (!is.na(predrug_latest_beta_blockers) & as.numeric(difftime(dstartdate, predrug_latest_beta_blockers)/30) < 6) |
          (!is.na(predrug_latest_calcium_channel_blockers) & as.numeric(difftime(dstartdate, predrug_latest_calcium_channel_blockers)/30) < 6) |
          (!is.na(predrug_latest_thiazide_diuretics) & as.numeric(difftime(dstartdate, predrug_latest_thiazide_diuretics)/30) < 6) |
          (!is.na(predrug_latest_loop_diuretics) & as.numeric(difftime(dstartdate, predrug_latest_loop_diuretics)/30) < 6) |
          (!is.na(predrug_latest_ksparing_diuretics) & as.numeric(difftime(dstartdate, predrug_latest_ksparing_diuretics)/30) < 6) |
          (!is.na(predrug_latest_arb) & as.numeric(difftime(dstartdate, predrug_latest_arb)/30) < 6),
        1,
        0
      ),
      
      # Cardiovascular event
      predrug_cardio_event = ifelse(
        predrug_angina == 1 | predrug_myocardialinfarction == 1 | predrug_ihd == 1 | predrug_pad == 1 | predrug_revasc == 1 | predrug_stroke == 1,
        1,
        0
      ),
      
      # Heart problem event
      predrug_heart_event = ifelse(
        predrug_heartfailure == 1 | predrug_hypertension == 1,
        1,
        0
      ),
      
      # Microvascular event
      predrug_micro_event = ifelse(
        predrug_retinopathy == 1 | predrug_diabeticnephropathy == 1 | predrug_neuropathy == 1,
        1,
        0
      ),
      
      # Frailty proxy event
      predrug_frailty_proxy = ifelse(
        predrug_falls == 1 | predrug_lowerlimbfracture == 1,
        1,
        0
      )
      
    )
  
  
  
  #####################################################################################
  #####################################################################################
  
  ## If the follow-up is 3-months:
  if (follow_up %in% c("3-months")) {
    
    
    # Select variables
    
    cprd_dataset <- cprd_dataset %>%
      select(
        # patient info
        patid, dstartdate,
        # Outcome
        stopdrug_3m_3mFU,
        stopdrug_6m_3mFU,
        stopdrug_12m_3mFU,
        # Drug taken
        drugclass, drugsubstances, drugcombo,
        # Extra info
        only_one_prescription,
        dstartdate_dm_dur, dstartdate_age, drugline, numdrugs, smoking_cat, imd2015_10,
        predrug_statins, stopdrug_3m_3mFU_MFN_hist, ethnicity_5cat, gender, predrug_bloodmed,
        # Biomarkers
        prehba1c, preegfr, prebmi, 
        prehdl,
        # Comorbidities
        ## Hist of cardiovascular
        predrug_cardio_event,
        predrug_angina, predrug_myocardialinfarction, predrug_ihd, predrug_pad,
        predrug_revasc, predrug_stroke, 
        ## Heart problems
        predrug_heart_event,
        predrug_heartfailure, predrug_hypertension,
        ### Microvascular
        predrug_micro_event,
        predrug_retinopathy, predrug_diabeticnephropathy, predrug_neuropathy,
        ## CKD
        preckdstage, 
        ## CLD
        predrug_cld,
        ## Frailty proxy
        predrug_frailty_proxy
      ) %>%
      as.data.frame() %>%
      mutate(
        imd2015_10 = ifelse(imd2015_10 %in% c(1, 2), 1, ifelse(imd2015_10 %in% c(3, 4), 2, ifelse(imd2015_10 %in% c(5, 6), 3, ifelse(imd2015_10 %in% c(7, 8), 4, 5))))
      ) %>%
      mutate_at(
        c(
          # Outcome
          "stopdrug_3m_3mFU",
          "stopdrug_6m_3mFU",
          "stopdrug_12m_3mFU",
          # Extra info
          "only_one_prescription",
          "smoking_cat", "imd2015_10",
          "predrug_statins", "ethnicity_5cat", "gender", "predrug_bloodmed",
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
          # ## CKD
          # "preckdstage", 
          ## CLD
          "predrug_cld",
          ## Frailty proxy
          "predrug_frailty_proxy"
        ),
        as.factor
      ) %>%
      mutate(
        preckdstage = ifelse(is.na(preckdstage), "stage_0", preckdstage),
        preckdstage = factor(preckdstage, levels = c("stage_0", "stage_1", "stage_2", "stage_3a", "stage_3b", "stage_4", "stage_5"), labels = c("stage_0", "stage_1", "stage_2", "stage_3a", "stage_3b", "stage_4", "stage_5")),
        
        drugclass = factor(drugclass, levels = c("MFN", "GLP1", "DPP4", "SGLT2", "SU", "TZD")),
        stopdrug_3m_3mFU_MFN_hist = ifelse(is.na(stopdrug_3m_3mFU_MFN_hist), 0, ifelse(stopdrug_3m_3mFU_MFN_hist > 0, 1, 0)),
        stopdrug_3m_3mFU_MFN_hist = factor(stopdrug_3m_3mFU_MFN_hist),
        
        # make drugline have a limit 5+ above 4
        drugline = ifelse(drugline > 4, 5, drugline),
        drugline = factor(drugline, levels = c(1, 2, 3, 4, 5), labels = c("1", "2", "3", "4", "5+")),
        
        # make numdrugs have a limit 3+ above 2
        numdrugs = ifelse(numdrugs > 2, 3, numdrugs),
        numdrugs = factor(numdrugs, levels = c(1, 2, 3), labels = c("1", "2", "3+"))
      )
    
    
    if (dataset == "full.dataset") {return(cprd_dataset)}
    
    
    
    
    #####################################################################################
    #####################################################################################
    
    # Remove patients without discontinuation data
    
    if (dataset %in% c("3m.disc.dataset", "3m.disc.dataset.dev", "3m.disc.dataset.val")) {
      
      cprd_dataset <- cprd_dataset %>%
        drop_na(stopdrug_3m_3mFU) %>%
        mutate(row = 1:n())
      
      if (dataset == "3m.disc.dataset") {return(cprd_dataset %>% select(-row))}
      
      set.seed(123)
      cprd_dataset.dev <- cprd_dataset %>%
        sample_frac(.7)
      
      cprd_dataset.val <- cprd_dataset %>%
        filter(!(row %in% cprd_dataset.dev$row))
      
      
      ## Check patients after data
      if (isTRUE(diagnosis)) {
        print("#####################################")
        print(paste("Patients with 3m discontinuation data: DEV", nrow(cprd_dataset.dev), ", VAL", nrow(cprd_dataset.val)))
        print("#####################################")
      }
      
      if (dataset == "3m.disc.dataset.dev") {return(cprd_dataset.dev)}
      if (dataset == "3m.disc.dataset.val") {return(cprd_dataset.val)}
      
      
    } else if (dataset %in% c("6m.disc.dataset", "6m.disc.dataset.dev", "6m.disc.dataset.val")) {
      
      cprd_dataset <- cprd_dataset %>%
        drop_na(stopdrug_6m_3mFU) %>%
        mutate(row = 1:n())
      
      if (dataset == "6m.disc.dataset") {return(cprd_dataset %>% select(-row))}
      
      set.seed(123)
      cprd_dataset.dev <- cprd_dataset %>%
        sample_frac(.7)
      
      cprd_dataset.val <- cprd_dataset %>%
        filter(!(row %in% cprd_dataset.dev$row))
      
      
      ## Check patients after data
      if (isTRUE(diagnosis)) {
        print("#####################################")
        print(paste("Patients with 6m discontinuation data: DEV", nrow(cprd_dataset.dev), ", VAL", nrow(cprd_dataset.val)))
        print("#####################################")
      }
      
      if (dataset == "6m.disc.dataset.dev") {return(cprd_dataset.dev)}
      if (dataset == "6m.disc.dataset.val") {return(cprd_dataset.val)}
      
      
    } else if (dataset %in% c("12m.disc.dataset", "12m.disc.dataset.dev", "12m.disc.dataset.val")) {
      
      cprd_dataset <- cprd_dataset %>%
        drop_na(stopdrug_12m_3mFU) %>%
        mutate(row = 1:n())
      
      if (dataset == "12m.disc.dataset") {return(cprd_dataset %>% select(-row))}
      
      set.seed(123)
      cprd_dataset.dev <- cprd_dataset %>%
        sample_frac(.7)
      
      cprd_dataset.val <- cprd_dataset %>%
        filter(!(row %in% cprd_dataset.dev$row))
      
      
      ## Check patients after data
      if (isTRUE(diagnosis)) {
        print("#####################################")
        print(paste("Patients with 12m discontinuation data: DEV", nrow(cprd_dataset.dev), ", VAL", nrow(cprd_dataset.val)))
        print("#####################################")
      }
      
      if (dataset == "12m.disc.dataset.dev") {return(cprd_dataset.dev)}
      if (dataset == "12m.disc.dataset.val") {return(cprd_dataset.val)}
      
      
    }
    
    
    
  }
  
  
  #####################################################################################
  #####################################################################################
  
  ## If the follow-up is 6-months:
  if (follow_up %in% c("6-months")) {
    
    # Select variables
    
    cprd_dataset <- cprd_dataset %>%
      select(
        # patient info
        patid, dstartdate,
        # Outcome
        stopdrug_3m_6mFU,
        stopdrug_6m_6mFU,
        stopdrug_12m_6mFU,
        # Drug taken
        drugclass, drugsubstances, drugcombo,
        # Extra info
        only_one_prescription,
        dstartdate_dm_dur, dstartdate_age, drugline, numdrugs, smoking_cat, imd2015_10,
        predrug_statins, stopdrug_3m_3mFU_MFN_hist, ethnicity_5cat, gender, predrug_bloodmed,
        # Biomarkers
        prehba1c, preegfr, prebmi, 
        prehdl,
        # Comorbidities
        ## Hist of cardiovascular
        predrug_cardio_event,
        predrug_angina, predrug_myocardialinfarction, predrug_ihd, predrug_pad,
        predrug_revasc, predrug_stroke, 
        ## Heart problems
        predrug_heart_event,
        predrug_heartfailure, predrug_hypertension,
        ### Microvascular
        predrug_micro_event,
        predrug_retinopathy, predrug_diabeticnephropathy, predrug_neuropathy,
        ## CKD
        preckdstage, 
        ## CLD
        predrug_cld,
        ## Frailty proxy
        predrug_frailty_proxy
      ) %>%
      as.data.frame() %>%
      mutate(
        imd2015_10 = ifelse(imd2015_10 %in% c(1, 2), 1, ifelse(imd2015_10 %in% c(3, 4), 2, ifelse(imd2015_10 %in% c(5, 6), 3, ifelse(imd2015_10 %in% c(7, 8), 4, 5))))
      ) %>%
      mutate_at(
        c(
          # Outcome
          "stopdrug_3m_6mFU",
          "stopdrug_6m_6mFU",
          "stopdrug_12m_6mFU",
          # Extra info
          "only_one_prescription",
          "smoking_cat", "imd2015_10",
          "predrug_statins", "ethnicity_5cat", "gender", "predrug_bloodmed",
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
          # ## CKD
          # "preckdstage", 
          ## CLD
          "predrug_cld",
          ## Frailty proxy
          "predrug_frailty_proxy"
        ),
        as.factor
      ) %>%
      mutate(
        preckdstage = ifelse(is.na(preckdstage), "stage_0", preckdstage),
        preckdstage = factor(preckdstage, levels = c("stage_0", "stage_1", "stage_2", "stage_3a", "stage_3b", "stage_4", "stage_5"), labels = c("stage_0", "stage_1", "stage_2", "stage_3a", "stage_3b", "stage_4", "stage_5")),
        
        drugclass = factor(drugclass, levels = c("MFN", "GLP1", "DPP4", "SGLT2", "SU", "TZD")),
        stopdrug_3m_3mFU_MFN_hist = ifelse(is.na(stopdrug_3m_3mFU_MFN_hist), 0, ifelse(stopdrug_3m_3mFU_MFN_hist > 0, 1, 0)),
        stopdrug_3m_3mFU_MFN_hist = factor(stopdrug_3m_3mFU_MFN_hist),
        
        # make drugline have a limit 5+ above 4
        drugline = ifelse(drugline > 4, 5, drugline),
        drugline = factor(drugline, levels = c(1, 2, 3, 4, 5), labels = c("1", "2", "3", "4", "5+")),
        
        # make numdrugs have a limit 3+ above 2
        numdrugs = ifelse(numdrugs > 2, 3, numdrugs),
        numdrugs = factor(numdrugs, levels = c(1, 2, 3), labels = c("1", "2", "3+"))
      )
    
    
    if (dataset == "full.dataset") {return(cprd_dataset)}
    
    
    
    
    #####################################################################################
    #####################################################################################
    
    # Remove patients without discontinuation data
    
    if (dataset %in% c("3m.disc.dataset", "3m.disc.dataset.dev", "3m.disc.dataset.val")) {
      
      cprd_dataset <- cprd_dataset %>%
        drop_na(stopdrug_3m_6mFU) %>%
        mutate(row = 1:n())
      
      if (dataset == "3m.disc.dataset") {return(cprd_dataset %>% select(-row))}
      
      set.seed(123)
      cprd_dataset.dev <- cprd_dataset %>%
        sample_frac(.7)
      
      cprd_dataset.val <- cprd_dataset %>%
        filter(!(row %in% cprd_dataset.dev$row))
      
      
      ## Check patients after data
      if (isTRUE(diagnosis)) {
        print("#####################################")
        print(paste("Patients with 3m discontinuation data: DEV", nrow(cprd_dataset.dev), ", VAL", nrow(cprd_dataset.val)))
        print("#####################################")
      }
      
      if (dataset == "3m.disc.dataset.dev") {return(cprd_dataset.dev)}
      if (dataset == "3m.disc.dataset.val") {return(cprd_dataset.val)}
      
      
    } else if (dataset %in% c("6m.disc.dataset", "6m.disc.dataset.dev", "6m.disc.dataset.val")) {
      
      cprd_dataset <- cprd_dataset %>%
        drop_na(stopdrug_6m_6mFU) %>%
        mutate(row = 1:n())
      
      if (dataset == "6m.disc.dataset") {return(cprd_dataset %>% select(-row))}
      
      set.seed(123)
      cprd_dataset.dev <- cprd_dataset %>%
        sample_frac(.7)
      
      cprd_dataset.val <- cprd_dataset %>%
        filter(!(row %in% cprd_dataset.dev$row))
      
      
      ## Check patients after data
      if (isTRUE(diagnosis)) {
        print("#####################################")
        print(paste("Patients with 6m discontinuation data: DEV", nrow(cprd_dataset.dev), ", VAL", nrow(cprd_dataset.val)))
        print("#####################################")
      }
      
      if (dataset == "6m.disc.dataset.dev") {return(cprd_dataset.dev)}
      if (dataset == "6m.disc.dataset.val") {return(cprd_dataset.val)}
      
      
    } else if (dataset %in% c("12m.disc.dataset", "12m.disc.dataset.dev", "12m.disc.dataset.val")) {
      
      cprd_dataset <- cprd_dataset %>%
        drop_na(stopdrug_12m_6mFU) %>%
        mutate(row = 1:n())
      
      if (dataset == "12m.disc.dataset") {return(cprd_dataset %>% select(-row))}
      
      set.seed(123)
      cprd_dataset.dev <- cprd_dataset %>%
        sample_frac(.7)
      
      cprd_dataset.val <- cprd_dataset %>%
        filter(!(row %in% cprd_dataset.dev$row))
      
      
      ## Check patients after data
      if (isTRUE(diagnosis)) {
        print("#####################################")
        print(paste("Patients with 12m discontinuation data: DEV", nrow(cprd_dataset.dev), ", VAL", nrow(cprd_dataset.val)))
        print("#####################################")
      }
      
      if (dataset == "12m.disc.dataset.dev") {return(cprd_dataset.dev)}
      if (dataset == "12m.disc.dataset.val") {return(cprd_dataset.val)}
      
      
    }
    
  }
  
  
}



#:--------------------------------------------------


calc_ATE <- function(data, drugs, pred.variable = ".no_weight", weight.variable = NULL, break_points = NULL, breakdown = NULL, matching = FALSE, ntiles = 10, caliper = 0.05, replace = FALSE, order = "random", n_bootstrap = 10) {
  
  # output object
  final_dataset <- NULL
  
  # calculation of all drug combinations (from provided)
  combinations_matrix <- combn(drugs, 2)
  
  # number of combinations
  n_comb <- ncol(combinations_matrix)
  
  # break points
  if (!is.null(break_points)) {break_points <- c(-1, break_points, 1)}
  
  # iterate through each combinations
  for (comb_i in 1:n_comb) {
    
    # current drugs being analysed
    current_drugs <- combinations_matrix[,comb_i]
    
    if (is.null(break_points)) {
      
      # collect only the drugs we are interested
      data.initial <- data %>%
        filter(drugclass %in% current_drugs) %>%
        rename(
          "pred.drug_1" = paste0("pred", pred.variable, ".", current_drugs[1]),
          "pred.drug_2" = paste0("pred", pred.variable, ".", current_drugs[2]),
        ) %>%
        mutate(diff = pred.drug_1 - pred.drug_2) %>%
        mutate(group = ntile(diff, ntiles))
      
    } else {
      
      
      # collect only the drugs we are interested
      data.initial <- data %>%
        filter(drugclass %in% current_drugs) %>%
        rename(
          "pred.drug_1" = paste0("pred", pred.variable, ".", current_drugs[1]),
          "pred.drug_2" = paste0("pred", pred.variable, ".", current_drugs[2]),
        ) %>%
        mutate(diff = pred.drug_1 - pred.drug_2) %>%
        mutate(group = cut(diff, breaks = break_points)) %>%
        mutate(group = factor(group))
      
    }
    
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
      ) %>%
      left_join(
        data.initial %>%
          group_by(group) %>%
          filter(drugclass == current_drugs[1]) %>%
          mutate(
            N_drug1 = n()
          ) %>%
          ungroup() %>%
          select(group, N_drug1) %>%
          unique(),
        by = c("group")
      ) %>%
      left_join(
        data.initial %>%
          group_by(group) %>%
          filter(drugclass == current_drugs[2]) %>%
          mutate(
            N_drug2 = n()
          ) %>%
          ungroup() %>%
          select(group, N_drug2) %>%
          unique(),
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
      
      if (length(unique(data.new$drugclass)) != 2 | length(unique(data.new$stopdrug_3m_6mFU)) != 2) {
        
        # collect mean
        obs <- append(obs, NA)
        
        # collect lower bound CI
        lci <- append(lci, NA)
        
        # collect upper bound CI
        uci <- append(uci, NA)
        
      } else {
        
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
        
        
        # collect bootstrapped values
        bootstrap_differences <- NULL
        
        
        for (bootstrap_i in 1:n_bootstrap) {
          
          if (matching == FALSE) {
            
            # bootstrap dataset
            data.bootstrap <- data.new[sample(nrow(data.new), size = nrow(data.new), replace = TRUE),] 
            
            if (length(unique(data.bootstrap$drugclass)) != 2 | length(unique(data.bootstrap$stopdrug_3m_6mFU)) != 2) {
              
              # save differences
              bootstrap_differences <- append(bootstrap_differences, NA)
              
            } else {
              
              if (!is.null(weight.variable)) {
                
                # fit linear regression for decile in the matched dataset
                models[[i]] <- glm(as.formula("stopdrug_3m_6mFU ~ factor(drugclass)"), data=data.bootstrap, weights = propensity.score, family = quasibinomial())
                
              } else {
                
                # fit linear regression for decile in the matched dataset
                models[[i]] <- glm(as.formula(formula),data=data.bootstrap, family = binomial())
                
                # make predictions
                preds.drug1 <- predict(models[[i]], newdata = data.bootstrap %>% filter(drugclass == current_drugs[1]), type = "response")
                
                preds.drug2 <- predict(models[[i]], newdata = data.bootstrap %>% filter(drugclass == current_drugs[2]), type = "response")
                
                # save differences
                bootstrap_differences <- append(bootstrap_differences, mean(preds.drug1) - mean(preds.drug2))
                
              }
              
            }
            
            
          } else {
            
            if (!is.null(weight.variable)) {
              
              # model if propensity scores are provided
              matching_package_result <- MatchIt::matchit(
                formula = formula("drugclass ~ prehba1c"), # shouldn't be used since we are specifying 'distance' (propensity scores)
                data = data.new[which(data.new[,"group"] == i),], # select people in the quantile
                method = "nearest",
                distance = data.new[which(data.new[,"group"] == i),"propensity.score"],
                replace = replace,
                m.order = order,
                caliper = caliper,
                mahvars = NULL, estimand = "ATT", exact = NULL, antiexact = NULL, discard = "none", reestimate = FALSE, s.weights = NULL, std.caliper = TRUE, ratio = 1, verbose = FALSE, include.obj = FALSE,
              )
              
              data.match <- data.new %>%
                mutate(weights.matching = matching_package_result$weights) %>%
                filter(weights.matching > 0)
              
              # bootstrap dataset
              data.bootstrap <- data.match[sample(nrow(data.match), size = nrow(data.match), replace = TRUE),] 
              
              if (length(unique(data.bootstrap$drugclass)) != 2 | length(unique(data.bootstrap$stopdrug_3m_6mFU)) != 2) {
                
                # save differences
                bootstrap_differences <- append(bootstrap_differences, NA)
                
              } else {
                
                # fit linear regression for decile in the matched dataset
                models[[i]] <- glm(as.formula(formula),data=data.bootstrap, family = binomial())
                
                # make predictions
                preds.drug1 <- predict(models[[i]], newdata = data.bootstrap %>% filter(drugclass == current_drugs[1]), type = "response")
                
                preds.drug2 <- predict(models[[i]], newdata = data.bootstrap %>% filter(drugclass == current_drugs[2]), type = "response")
                
                # save differences
                bootstrap_differences <- append(bootstrap_differences, mean(preds.drug1) - mean(preds.drug2))
                
              }
              
            } else {
              
              stop("Please provide 'weight.variable' for matching.")
              
            }
            
          }
          
        }
        
        # collect mean
        obs <- append(obs, quantile(bootstrap_differences, probs = c(0.5), na.rm = TRUE))
        
        # collect lower bound CI
        lci <- append(lci, quantile(bootstrap_differences, probs = c(0.025), na.rm = TRUE))
        
        # collect upper bound CI
        uci <- append(uci, quantile(bootstrap_differences, probs = c(0.975), na.rm = TRUE))
        
      }
      
      
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




calc_odds_ratios <- function(data, drugs, pred.variable = ".no_weight", break_points = c(-0.1, -0.05, -0.025, 0, 0.025, 0.05, 0.1)) {
  
  # output object
  final_dataset <- NULL
  
  # calculation of all drug combinations (from provided)
  combinations_matrix <- combn(drugs, 2)
  
  # number of combinations
  n_comb <- ncol(combinations_matrix)
  
  break_points <- c(-1, break_points, 1)
  
  # iterate through each combinations
  for (comb_i in 1:n_comb) {
    
    # current drugs being analysed
    current_drugs <- combinations_matrix[,comb_i]
    
    # collect only the drugs we are interested
    data.initial <- data %>%
      filter(drugclass %in% current_drugs) %>%
      mutate(drugclass = factor(drugclass, levels = rev(current_drugs))) %>%
      rename(
        "pred.drug_1" = paste0("pred", pred.variable, ".", current_drugs[1]),
        "pred.drug_2" = paste0("pred", pred.variable, ".", current_drugs[2]),
      ) %>%
      mutate(diff = pred.drug_1 - pred.drug_2) %>%
      mutate(group = cut(diff, breaks = break_points)) %>%
      mutate(group = factor(group))
    
    # calculate the summary statistics for each ntile
    predicted_treatment_effect <- data.initial %>%
      plyr::ddply("group", dplyr::summarise,
                  N = length(diff)) %>%
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
    
    predicted_treatment_effect <- predicted_treatment_effect %>%
      mutate(
        group = factor(group, levels = rev(levels(predicted_treatment_effect$group)))
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
      
      if (length(unique(data.new$drugclass)) != 2 | length(unique(data.new$stopdrug_3m_6mFU)) != 2) {
        
        # collect mean
        obs <- append(obs, NA)
        
        # collect lower bound CI
        lci <- append(lci, NA)
        
        # collect upper bound CI
        uci <- append(uci, NA)
        
      } else {
        
        # formula
        formula <- "stopdrug_3m_6mFU ~ factor(drugclass)"
        
        # fit linear regression for decile in the matched dataset
        models[[i]] <- glm(as.formula(formula),data=data.new, family = binomial())
        
        
        # iterate through bootstrapped datasets to get CI
        
        # collect mean
        obs <- append(obs, models[[i]]$coefficients[[2]])
        
        # collect lower bound CI
        lci <- append(lci, confint(models[[i]])[2,1])
        
        # collect upper bound CI
        uci <- append(uci, confint(models[[i]])[2,2])
        
      }
      
    }
    
    # join treatment effects for deciles in a data.frame
    effects <- data.frame(predicted_treatment_effect,cbind(obs, lci, uci)) %>%
      mutate(
        obs = as.numeric(obs),
        lci = as.numeric(lci),
        uci = as.numeric(uci)
      )
    
    # returned list with fitted propensity model + decile treatment effects
    final_dataset <- rbind(
      final_dataset,
      effects %>% mutate(Type = paste(current_drugs[1], "vs", current_drugs[2]))
    )
    
  }
  
  return(final_dataset)
  
}



calc_predicted_risk <- function(data, drugs, pred.variable = ".no_weight", weight.variable = NULL, break_points = NULL, breakdown = NULL, matching = FALSE, ntiles = 10, caliper = 0.05, replace = FALSE, order = "random") {
  
  # output object
  final_dataset <- NULL
  
  # calculation of all drug combinations (from provided)
  combinations_matrix <- combn(drugs, 2)
  
  # number of combinations
  n_comb <- ncol(combinations_matrix)
  
  # break points
  if (!is.null(break_points)) {break_points <- c(-1, break_points, 1)}
  
  # iterate through each combinations
  for (comb_i in 1:n_comb) {
    
    # current drugs being analysed
    current_drugs <- combinations_matrix[,comb_i]
    
    if (is.null(break_points)) {
      
      # collect only the drugs we are interested
      data.initial <- data %>%
        filter(drugclass %in% current_drugs) %>%
        rename(
          "pred.drug_1" = paste0("pred", pred.variable, ".", current_drugs[1]),
          "pred.drug_2" = paste0("pred", pred.variable, ".", current_drugs[2]),
        ) %>%
        mutate(diff = pred.drug_1 - pred.drug_2) %>%
        mutate(group = ntile(diff, ntiles))
      
    } else {
      
      
      # collect only the drugs we are interested
      data.initial <- data %>%
        filter(drugclass %in% current_drugs) %>%
        rename(
          "pred.drug_1" = paste0("pred", pred.variable, ".", current_drugs[1]),
          "pred.drug_2" = paste0("pred", pred.variable, ".", current_drugs[2]),
        ) %>%
        mutate(diff = pred.drug_1 - pred.drug_2) %>%
        mutate(group = cut(diff, breaks = break_points)) %>%
        mutate(group = factor(group))
      
    }
    
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
      ) %>%
      left_join(
        data.initial %>%
          group_by(group) %>%
          filter(drugclass == current_drugs[1]) %>%
          mutate(
            N_drug1 = n()
          ) %>%
          ungroup() %>%
          select(group, N_drug1) %>%
          unique(),
        by = c("group")
      ) %>%
      left_join(
        data.initial %>%
          group_by(group) %>%
          filter(drugclass == current_drugs[2]) %>%
          mutate(
            N_drug2 = n()
          ) %>%
          ungroup() %>%
          select(group, N_drug2) %>%
          unique(),
        by = c("group")
      )
    
    
    # maximum number of deciles being tested
    quantiles <- length(unique(data.initial[,"group"]))
    
    # create lists with results
    mnumber = c(1:quantiles)
    models  <- as.list(1:quantiles)
    drug1_obs <- vector(); drug1_lci <- vector(); drug1_uci <- vector()
    drug2_obs <- vector(); drug2_lci <- vector(); drug2_uci <- vector()
    
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
      
      if (length(unique(data.new$drugclass)) != 2 | length(unique(data.new$stopdrug_3m_6mFU)) != 2) {
        
        # collect mean
        drug1_obs <- append(drug1_obs, NA)
        
        # collect lower bound CI
        drug1_lci <- append(drug1_lci, NA)
        
        # collect upper bound CI
        drug1_uci <- append(drug1_uci, NA)
        
        # collect mean
        drug2_obs <- append(drug2_obs, NA)
        
        # collect lower bound CI
        drug2_lci <- append(drug2_lci, NA)
        
        # collect upper bound CI
        drug2_uci <- append(drug2_uci, NA)
        
      } else {
        
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
        
        
        if (matching == FALSE) {
          
          if (!is.null(weight.variable)) {
            
            # fit linear regression for decile in the matched dataset
            models[[i]] <- glm(as.formula("stopdrug_3m_6mFU ~ factor(drugclass)"), data=data.new, weights = propensity.score, family = quasibinomial())
            
            
            stop("Not coded")
            
            
          } else {
            
            # fit linear regression for decile in the matched dataset
            models[[i]] <- glm(as.formula(formula),data=data.new, family = binomial())
            
            # make predictions
            preds.drug1 <- predict(models[[i]], newdata = data.new %>% filter(drugclass == current_drugs[1]) %>% slice(1), type = "response", se.fit = TRUE)
            
            preds.drug2 <- predict(models[[i]], newdata = data.new %>% filter(drugclass == current_drugs[2]) %>% slice(1), type = "response", se.fit = TRUE)
            
          }
          
          
          # collect mean
          drug1_obs <- append(drug1_obs, preds.drug1$fit)
          
          # collect lower bound CI
          if (preds.drug1$fit - 1.96*(preds.drug1$se.fit) < 0) {drug1_lci <- append(drug1_lci, 0)} else {drug1_lci <- append(drug1_lci, preds.drug1$fit - 1.96*(preds.drug1$se.fit))}
          
          # collect upper bound CI
          if (preds.drug1$fit + 1.96*(preds.drug1$se.fit) < 0) {drug1_uci <- append(drug1_uci, 0)} else {drug1_uci <- append(drug1_uci, preds.drug1$fit + 1.96*(preds.drug1$se.fit))}
          
          # collect mean
          drug2_obs <- append(drug2_obs, preds.drug2$fit)
          
          # collect lower bound CI
          if (preds.drug2$fit - 1.96*(preds.drug2$se.fit) < 0) {drug2_lci <- append(drug2_lci, 0)} else {drug2_lci <- append(drug2_lci, preds.drug2$fit - 1.96*(preds.drug2$se.fit))}
          
          # collect upper bound CI
          if (preds.drug2$fit + 1.96*(preds.drug2$se.fit) < 0) {drug2_uci <- append(drug2_uci, 0)} else {drug2_uci <- append(drug2_uci, preds.drug2$fit + 1.96*(preds.drug2$se.fit))}
          
          
          
        } else {
          
          stop("Not coded")
          
        }
        
      }
      
    }
    # join treatment effects for deciles in a data.frame
    effects <- data.frame(predicted_treatment_effect,cbind(obs = drug1_obs, lci = drug1_lci, uci = drug1_uci, drug_type = current_drugs[1])) %>%
      rbind(
        data.frame(predicted_treatment_effect,cbind(obs = drug2_obs, lci = drug2_lci, uci = drug2_uci, drug_type = current_drugs[2]))
      ) %>%
      mutate(
        obs = as.numeric(obs),
        lci = as.numeric(lci),
        uci = as.numeric(uci)
      )
    
    # returned list with fitted propensity model + decile treatment effects
    final_dataset <- rbind(
      final_dataset,
      effects %>% mutate(Type = paste(current_drugs[1], "vs", current_drugs[2]))
    )
    
  }
  
  return(final_dataset)
  
}