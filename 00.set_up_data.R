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
    dataset = "full.dataset"
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
  if(!(dataset %in% c("full.dataset", "ps.dataset"))) {stop("'dataset' needs to be: full.dataset / ps.dataset / ")}
  
  
  ###############################################
  
  
  # load dataset - 1st instance
  load(paste0("/slade/CPRD_data/mastermind_2022/", raw_data, ".Rda"))   # name: t2d_1stinstance
  
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
  
  # Keep patients above 2014-01-01
  
  cprd_dataset <- cprd_dataset %>%
    mutate(dstartdate = as.Date(dstartdate)) %>%
    filter(dstartdate > "2014-01-01")
  
  
  ## Check patients after data
  if (isTRUE(diagnosis)) {
    print("#####################################")
    print(paste("Patients after 2014-01-01:", nrow(cprd_dataset)))
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
  
  # Remove patients starting multiple drugs
  
  cprd_dataset <- cprd_dataset %>%
    filter(multi_drug_start == 0)
  
  
  ## Check patients after data
  if (isTRUE(diagnosis)) {
    print("#####################################")
    print(paste("Patients initiating therapies:", nrow(cprd_dataset)))
    print("#####################################")
  }
  
  
  #####################################################################################
  #####################################################################################

  # Remove patients with HbA1c above 53
  
  cprd_dataset <- cprd_dataset %>%
    filter(is.na(prehba1c) | prehba1c > 53)
  
  
  ## Check patients after data
  if (isTRUE(diagnosis)) {
    print("#####################################")
    print(paste("Patients HbA1c above 53:", nrow(cprd_dataset)))
    print("#####################################")
  }
  
  
  #####################################################################################
  #####################################################################################
  
  # Remove patients with CKD stage 3 or above
  
  cprd_dataset <- cprd_dataset %>%
    filter(preckdstage %in% c("stage_1", "stage_2"))
  
  ## Check patients after data
  if (isTRUE(diagnosis)) {
    print("#####################################")
    print(paste("Patients CKD stage 1 or 2:", nrow(cprd_dataset)))
    print("#####################################")
  }
  
  
  
  # Create necessary variables
  
  ### selecting variables: Missing - number of chronic illnesses q: just sum over the vars listed below?
  cprd_dataset <- cprd_dataset %>%
    mutate(
      
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
        !is.na(predrug_angina) | !is.na(predrug_myocardialinfarction) | !is.na(predrug_ihd) | !is.na(predrug_pad) | !is.na(predrug_revasc) | !is.na(predrug_stroke),
        1,
        0
      ),
      
      # Heart problem event
      predrug_heart_event = ifelse(
        !is.na(predrug_heartfailure) | !is.na(predrug_hypertension),
        1,
        0
      ),
      
      # Microvascular event
      predrug_micro_event = ifelse(
        !is.na(predrug_retinopathy) | !is.na(predrug_diabeticnephropathy) | !is.na(predrug_neuropathy),
        1,
        0
      )
      
    )
  
  
  
  #####################################################################################
  #####################################################################################
  
  
  # Select variables
  
  cprd_dataset <- cprd_dataset %>%
    select(
      # patient info
      patid, dstartdate,
      # Outcome
      stopdrug_3m_6mFU,
      # Drug taken
      drugclass, drugsubstances, drugcombo,
      # Extra info
      dstartdate_dm_dur, dstartdate_age, drugline, numdrugs, smoking_cat, imd2015_10,
      predrug_statins, stopdrug_3m_3mFU_MFN_hist, ethnicity_5cat, gender, predrug_bloodmed,
      # Biomarkers
      prehba1c, preegfr, prebmi, prealt, 
      pretotalcholesterol,
      # preldl, prehdl, pretriglyceride,
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
      
       
    ) %>%
    as.data.frame() %>%
    mutate(
      imd2015_10 = ifelse(imd2015_10 %in% c(1, 2), 1, ifelse(imd2015_10 %in% c(3, 4), 2, ifelse(imd2015_10 %in% c(5, 6), 3, ifelse(imd2015_10 %in% c(7, 8), 4, 5))))
    ) %>%
    mutate_at(
      c(
        # Outcome
        "stopdrug_3m_6mFU",
        # Extra info
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
        ## CKD
        "preckdstage", 
        ## CLD
        "predrug_cld"
        ),
      as.factor
    ) %>%
    mutate(
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
  
  
}








