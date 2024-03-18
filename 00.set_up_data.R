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
  if(!(dataset %in% c("full.dataset", "ps.dataset"))) {stop("'dataset' needs to be: full.dataset / ps.dataset")}
  
  
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
  
  
  # Create necessary variables
  
  ### selecting variables: Missing - number of chronic illnesses q: just sum over the vars listed below?
  cprd_dataset <- cprd_dataset %>%
    mutate(
      # Month of start
      dstartmonth = format(as.Date(dstartdate, format="%d/%m/%Y"),"%m"),
      
      
      # Charlson Comorbidity Index (CCI) https://www.mdcalc.com/calc/3917/charlson-comorbidity-index-cci
      ## Age limits <50 0, 50-59 1, 60-69 2, 70-79 3, >=80 4
      CCI_index = ifelse(dstartdate_age < 50, 0, ifelse(dstartdate_age < 60, 1, ifelse(dstartdate_age < 70, 2, ifelse(dstartdate_age < 80, 3, 4)))),
      ## Myocardial infarction Yes +1
      CCI_index = ifelse(predrug_incident_mi == 1, CCI_index + 1, CCI_index),
      ## Congestive Heart Failure -> using Heart failure Yes +1  (we have predrug_primary_hhf and predrug_heartfailure, will use predrug_heartfailure)
      CCI_index = ifelse(predrug_heartfailure == 1, CCI_index + 1, CCI_index),
      ## Peripheral vascular disease Yes +1
      # CCI_index = 
      ## Cerebrovascular accident or Transient ischemic attack -> using predrug_tia Yes +1
      CCI_index = ifelse(predrug_tia == 1, CCI_index + 1, CCI_index),
      ## Dementia Yes +1
      CCI_index = ifelse(predrug_dementia == 1, CCI_index + 1, CCI_index),
      ## CPOD - chronic obstructive pulmonary disease Yes +1
      CCI_index = ifelse(predrug_copd == 1, CCI_index + 1, CCI_index),
      ## Connective tissue disease Yes +1
      # CCI_index = 
      ## Peptic ulcer disease Yes +1
      # CCI_index = 
      ## Liver disease None 0, Mild +1, Moderate to severe +3
      # CCI_index = 
      ## Diabetes mellitus None or diet 0, Uncomplicated +1, End-organ damage +2
      CCI_index = CCI_index + 1
      ## Hemiplegia Yes +2
      # CCI_index = 
      ## Moderate to severe CKD (Severe = on dialysis, status post kidney transplant, uremia, moderate = creatinine >3 mg/dL (0.27 mmol/L)) Yes +2
      # CCI_index = 
      ## Solid tumor None 0, Localized +2, Metastatic +6
      # CCI_index = 
      ## Leukaemia Yes +2
      # CCI_index = 
      ## Lymphoma Yes +2
      # CCI_index = 
      ## AIDS Yes +2
      # CCI_index = 
    )
  
  
  
  #####################################################################################
  #####################################################################################
  
  
  # Select variables
  
  cprd_dataset <- cprd_dataset %>%
    select(
      patid, dstartdate,
      # Outcome
      stopdrug_3m_6mFU, stopdrug_6m_6mFU,
      # Drug taken
      drugclass, drugsubstances, drugcombo,
      # Biomarkers
      prefastingglucose, prehdl, pretriglyceride, precreatinine_blood, preldl, 
      prealt, preast, pretotalcholesterol, predbp, presbp, preacr, prehba1c, 
      preegfr, prealbumin_blood, prebilirubin, prehaematocrit, prehaemoglobin, 
      prepcr,
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
      ins_in_1_year, prebmi, smoking_cat, stopdrug_3m_3mFU_MFN_hist,
      alcohol_cat, fh_diabetes, dstartdate_age, dstartdate_dm_dur, dstartmonth,
      CCI_index, drugline
    ) %>%
    as.data.frame() %>%
    mutate_at(
      c(
        # Outcome
        "stopdrug_3m_6mFU", "stopdrug_6m_6mFU",
        # Commorbidities
        "preckdstage", "predrug_frailty_mild", "predrug_frailty_moderate", 
        "predrug_frailty_severe", "predrug_primary_hhf", "predrug_af", "predrug_angina", 
        "predrug_asthma", "predrug_bronchiectasis", "predrug_cld", "predrug_copd", 
        "predrug_cysticfibrosis", "predrug_dementia", "predrug_diabeticnephropathy", 
        "predrug_fh_premature_cvd", 
        "predrug_haem_cancer", "predrug_heartfailure", 
        "predrug_hypertension", "predrug_ihd", "predrug_myocardialinfarction", 
        "predrug_neuropathy", "predrug_otherneuroconditions", "predrug_pad", 
        "predrug_pulmonaryfibrosis", "predrug_pulmonaryhypertension", 
        "predrug_retinopathy", "predrug_revasc", "predrug_rheumatoidarthritis", 
        "predrug_solid_cancer", "predrug_solidorgantransplant", "predrug_stroke", 
        "predrug_tia", "predrug_anxiety_disorders", "predrug_medspecific_gi",
        "predrug_benignprostatehyperplasia", "predrug_micturition_control",
        "predrug_volume_depletion", "predrug_urinary_frequency", "predrug_falls",
        "predrug_lowerlimbfracture", "predrug_incident_mi", "predrug_incident_stroke",
        "predrug_dka", "predrug_osteoporosis", "predrug_unstableangina", 
        "predrug_amputation",
        "hosp_admission_prev_year",
        # Extra info
        "gender", "prac_region", "ethnicity_5cat", "imd2015_10",
        "ins_in_1_year", "smoking_cat",
        "alcohol_cat", "fh_diabetes"),
      as.factor
    ) %>%
    mutate(
      drugclass = factor(drugclass, levels = c("MFN", "GLP1", "DPP4", "SGLT2", "SU", "TZD")),
      stopdrug_3m_3mFU_MFN_hist = ifelse(is.na(stopdrug_3m_3mFU_MFN_hist), 0, ifelse(stopdrug_3m_3mFU_MFN_hist > 0, 1, 0)),
      stopdrug_3m_3mFU_MFN_hist = factor(stopdrug_3m_3mFU_MFN_hist)
    )
  
  
  if (dataset == "full.dataset") {return(cprd_dataset)}
  
  
  
  #####################################################################################
  #####################################################################################
  
  
  cprd_dataset <- cprd_dataset %>%
    select(
      patid, dstartdate,
      # Outcome
      stopdrug_3m_6mFU,
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
      ins_in_1_year, prebmi, smoking_cat, drugline, stopdrug_3m_3mFU_MFN_hist,
      alcohol_cat, dstartdate_age, dstartdate_dm_dur, dstartmonth,
      CCI_index,
      drugclass
    ) %>%
    drop_na() %>%
    as.data.frame()
  
  
  if (dataset == "ps.dataset") {return(cprd_dataset)}
  
}


