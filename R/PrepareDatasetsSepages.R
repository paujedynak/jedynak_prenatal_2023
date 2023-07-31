# quiets concerns of R CMD check when variables appear in pipelines
utils::globalVariables(c("ident", "t2", "t3", "mother_bmi", "mother_edu",
                         "delivery_mode", "parity", "ch_sex", "corr_lmp",
                         "cy1hao1_q10", "date_delivery", "date_lmp", "ddg_echo_T1_1",
                         "ddg_echo_T1_2", "ddg_echo_T1_3", "echo_T1_1", "echo_T1_2",
                         "echo_T1_3", "mean_ddg_echo", "mo_age", "mo_he", "mo_par",
                         "mo_we_bepr", "mother_height", "mother_weight", "mt1saa1_q01_1",
                         "po_datedel", "po_datelmp", "t2_log2", "t3_log2", "mother_edu",
                         "season_conc", "child_sex", "mother_age", "gest_age",
                         "puce", "plaque", "batch"))

#' Create covariates data
#'
#' @param covariate_data A .dta file containing bdd_grossesse_v4 data
#' @param exposure_data A .csv file containing phenol exposure data
#' @param exposures A character vector defining original names of exposures (from t2 and t3)
#'
#' @return A data.frame containing covariates data
#' @export
#' @import dplyr
#' @importFrom lubridate month
#' @importFrom tidyselect all_of
#' @importFrom tidyr replace_na

PrepareDatasetsSepages <- function(covariate_data,
                                   technical_covariates_data,
                                   exposure_data,
                                   exposures) {
  
  # Process phenol exposures
  # select analysis exposures and log trans exposures
  expo_log <- exposure_data %>%
    
    # Select standardized exposures
    dplyr::select(ident, tidyselect::all_of(exposures)) %>%
    dplyr::mutate(ident = as.character(ident)) %>% 
    
    # Log2 standardized exposures
    dplyr::mutate_if(is.numeric, funs("log2" = log2))
  
  colnames(expo_log) <- c("id", sub('.*cor_', '', colnames(expo_log)[-1]))
  
  expo_log_av <- expo_log %>%
    
    # add mean value of exposure for T2 and T3
    dplyr::rowwise() %>%
    dplyr::mutate(mean_expo = mean(c(t2, t3), na.rm = TRUE),
                  TCS_log2 = mean(c(t2_log2, t3_log2), na.rm = TRUE),
                  mo_urine_sample_date_t2 = lubridate::ymd(mo_urine_sample_date_t2),
                  mo_urine_sample_date_t3 = lubridate::ymd(mo_urine_sample_date_t3))
  
  # Process covariates
  
  # Note for deleted ID 18147, email sarah 18/06/2020:
  # >"cette femme a eu deux grossesses et la grossesse Sepages n'a pas aboutit
  # elle a eu un 2eme enfant ensuite et a accouché en décembre 2016 (mais cet
  # enfant n'est pas un enfant Sepages).  Quand nosu sommes allées dans les
  # maternités récupérer les poids de naissance / taille manquants, l'equipe
  # Sepages a du récupérer les info de cette deuxième grossesse...
  # >Donc Matthieu : supprime cette volontaire de ton analyse."
  
  
  # create covariates data
  cov <- covariate_data %>%
    dplyr::mutate(
      id = as.character(ident),
      
      # Birth place (clinic)
      birth_place = factor(cy1hao1_q03,
                           
                           # Merge two last levels (Other)
                           labels = c("HCE", "Mutual", "Belled", "Other", "Other")),
      
      # Maternal age
      mother_age = as.numeric(mo_age),
      
      # Maternal active smoking
      maternal_smoke_bef_pregn = dplyr::case_when(mo_tob_avgr == 0 ~ 0,
                                                  mo_tob_avgr > 0 ~ 1),
      
      maternal_smoke_gest_yn = mo_tob_grstt1_yn,
      maternal_smoke_anytime_yn = mo_tob_gr_anyt_yn_n2,
      
      # Maternal education
      mother_edu =
        dplyr::case_when(
          mo_dipl %in% c(1, 2, 3) ~ "bac+2",
          mo_dipl == 4 ~ "bac+3,4",
          mo_dipl == 5 ~ "bac>5"
        ),
      mother_edu = factor(mother_edu, 
                          levels = c("bac+2", "bac+3,4", "bac>5")),
      
      # Maternal vitamin use
      mother_vitamin = dplyr::na_if(mo_vit_preg_c, "NA"),
      mother_vitamin = factor(mother_vitamin, 
                              labels = c("Never", "T1", "T2", "T3")),
      
      mother_folicA = dplyr::na_if(mo_vitfol_preg_c, "NA"),
      mother_folicA = factor(mother_folicA, 
                             labels = c("Never", "T1", "T2", "T3")),
      
      # Parity
      parity = factor(mo_par,
                      
                      # Combine two last levels One and more
                      labels = c("Nulliparous", "One_and_more", "One_and_more")),
      
      # Child sex
      child_sex = factor(ch_sex,
                         labels = c("Male", "Female")),
      
      # Echo based conception date
      echo_T1_1 = as.Date(ddg_echo_T1_1),
      echo_T1_2 = as.Date(ddg_echo_T1_2),
      echo_T1_3 = as.Date(ddg_echo_T1_3),
      
      # LMP date
      date_lmp = as.Date(po_datelmp),
      
      # Delivery date
      date_delivery = as.Date(po_datedel),
      
      # Maternal weight
      mother_weight = mo_we_bepr,
      
      # Maternal height
      mother_height = mo_he,
      
      # Maternal country of origin
      # mt1saa1_q01p2, mt1saa1_q01p3, mt1saa1_q01p4 are empty
      mother_origin_france = dplyr::case_when(
        mt1saa1_q01_1 == 1 ~ "France",
        is.na(mt1saa1_q01_1) & mt1saa1_q01p1 != "" ~ "Not_France",
        is.na(mt1saa1_q01_1) & mt1saa1_q01p5 != "" ~ "Not_France"),
      
      mother_origin_france = factor(mother_origin_france),
      
      # Paternal country of origin
      # ft2sac1_q01_4, ft2sac1_q01_5 are empty
      father_origin_france = dplyr::case_when(
        ft2sac1_q01_1 == 1 ~ "France",
        is.na(ft2sac1_q01_1) & !is.na(ft2sac1_q01_2) ~ "Not_France",
        is.na(ft2sac1_q01_1) & !is.na(ft2sac1_q01_3) ~ "Not_France"),
      
      father_origin_france = factor(father_origin_france),

      # Delivery mode general
      delivery_mode = factor(po_delmod,
                             labels = c("Vaginal", "Cesarean")),
      
      # Pre-pregnancy HTA
      pre_pregn_HTA = factor(mt1haa1_q18,
                             labels = c("No", "Yes")),
      
      # Gestational HTA
      gest_HTA = dplyr::case_when(
        mt1haa1_q40p1 == 0 ~ "No",
        mt1haa1_q40p1 == 1 ~ "Yes",
        mt3haf1_q01 == 0 ~ "No",
        mt3haf1_q01 == 1 ~ "Yes"),
      
      gest_HTA = factor(gest_HTA),
      
      # Pre-pregnancy and gestational HTA
      HTA = dplyr::case_when(
        pre_pregn_HTA == "Yes" | gest_HTA == "Yes" ~ "Yes",
        pre_pregn_HTA == "No" & gest_HTA == "No" ~ "No"),
      
      HTA = factor(HTA),
      
      # Gestational diabetes
      gest_diab = dplyr::case_when(
        mt1haa1_q40p7 == 0 ~ "No",
        mt1haa1_q40p7 == 1 ~ "Yes",
        mt3haf1_q04 == 0 ~ "No",
        mt3haf1_q04 == 1 ~ "Yes"),
      
      gest_diab = factor(gest_diab),
      
      # Maternal anxiety and depression score
      mother_anx_depr_tot = mo_hadtotscore_grt3_imp,
      
      .keep = "none") %>%
    
    dplyr::rowwise() %>%
    dplyr::mutate(mean_ddg_echo = mean.Date(c(echo_T1_1, echo_T1_2, echo_T1_3), na.rm = TRUE),
                  .keep = "unused") %>%

    dplyr::mutate(diff = as.numeric(mean_ddg_echo - (date_lmp + 14)),
                  
                  # Based on the date of the LMP or gestational duration assessed by the obstetrician if it differed
                  # from the LMP-based estimate by more than 2 weeks
                  corr_lmp = dplyr::case_when(abs(diff) > 14 ~ mean_ddg_echo,
                                              TRUE ~ date_lmp),
                  corr_lmp = lubridate::ymd(corr_lmp),
                  # one woman 17668 had the difference of >14 days. 
                  # I am not using abs(diff) <= 14 ~ date_lmp but TRUE ~ date_lmp to avoid situation that for missing diff 
                  # (and this happens when all three ddg_echo are missing) the corr_lmp would be NA. In my approach corr_lmp
                  # will be date_lmp
                  
                  # Calculate gestational duration based on the LMP corrected with US data
                  gest_age = date_delivery - corr_lmp,
                  gest_age = as.numeric(gest_age / 7),
                  
                  # Calculate pre-pregnancy BMI [weight (kg) / height (cm) / height (cm)] x 10,000
                  mother_bmi = (mother_weight / mother_height / mother_height) * 10000,
                  mother_bmi = dplyr::case_when(mother_bmi < 18.5 ~ "Underweight",
                                                mother_bmi >= 18.5 & mother_bmi < 25 ~ "Normal_weight",
                                                mother_bmi >= 25 & mother_bmi < 30 ~ "Overweight_obesity",
                                                mother_bmi >= 30 ~ "Overweight_obesity"),
                  mother_bmi = as.factor(mother_bmi),
                  mother_bmi = stats::relevel(mother_bmi, ref = "Normal_weight"),
                  
                  # Calculate season of conception
                  season_conc = case_when(
                    lubridate::month(corr_lmp + 14) %in% c(1:3) ~ "Jan_March",
                    lubridate::month(corr_lmp + 14) %in% c(4:6) ~ "April_June",
                    lubridate::month(corr_lmp + 14) %in% c(7:9) ~ "July_Sept",
                    lubridate::month(corr_lmp + 14) %in% c(10:12) ~ "Oct_Dec"),
                  season_conc = factor(season_conc, levels = c("Jan_March", "April_June", "July_Sept", "Oct_Dec")),
                  
                  maternal_smoke_sum = sum(maternal_smoke_bef_pregn, maternal_smoke_gest_yn, maternal_smoke_anytime_yn, na.rm = TRUE),
                  maternal_smoke_sum = dplyr::case_when(is.na(maternal_smoke_bef_pregn) & is.na(maternal_smoke_gest_yn) &
                                                          is.na(maternal_smoke_anytime_yn) ~ NA_real_,
                                                        TRUE ~ maternal_smoke_sum),
                  mother_act_smoke = case_when(maternal_smoke_sum == 0 ~ "Didn't_smoke",
                                               maternal_smoke_sum != 0 ~ "Smoked_before_OrAnd_in_pregn"), 
                  mother_act_smoke = factor(mother_act_smoke)
    )
  
  cov <- select(cov,
                -dplyr::contains("echo"),
                -dplyr::contains("date"),
                -dplyr::contains("maternal"),
                -mother_weight,
                -mother_height,
                -diff
  )
  
  # Merge exposures with covariates
  expo_cov <- merge(expo_log_av, cov, by = "id") %>%  # dim 479
    dplyr::mutate(gest_age_urine_sampl_t2 = round(as.double((mo_urine_sample_date_t2 - corr_lmp) / 7), 1),
                  gest_age_urine_sampl_t3 = round(as.double((mo_urine_sample_date_t3 - corr_lmp) / 7), 1))
  
  # Process technical confounders
  technical_conf <- technical_covariates_data %>% # dim 395   6
    
    # Select variables of interest. "Plaque" is the correct variable coding the plate number
    # Change chip, plate, and batch to factors
    dplyr::transmute(id = as.character(ident),
                     chip = factor(chip),
                     plate = factor(plate),
                     batch = as.factor(batch))
  
  # Merge exposures-covariates with technical factors
  expo_cov_tech_conf <- merge(expo_cov, technical_conf, by = "id") # %>% # dim 395  27
  
  return(expo_cov_tech_conf)
}

#' Impute missing data in categorical covariates (with mode)
#'
#' @param data_w_missing A data.frame with covariates containing missing values
#'
#' @return A data.frame with covariates with missing values imputed
#' @export

ImputeMissingCov <- function(data_w_missing) {
  
  # df %>% mutate_if(is.numeric, funs(replace(.,is.na(.), mean(., na.rm = TRUE)))) %>%
  #   mutate_if(is.factor, funs(replace(.,is.na(.), Mode(na.omit(.)))))
  
  Mode <- function(x) {
    ux <- unique(x)
    ux[which.max(tabulate(match(x, ux)))]
  }
  
  data_complete <- data_w_missing %>%
    dplyr::mutate_if(is.factor, funs(tidyr::replace_na(., Mode(na.omit(.)))))
  
  return(data_complete)
}
