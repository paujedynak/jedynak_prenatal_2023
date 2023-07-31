#' Main mediation analysis
#'
#' @param CpG A numeric vector containing methylation data from 1 CpG
#' @param expo_cov A data.frame containing all variables needed for regression (exposures and confounders)
#' @param mediator A string defining mediator name
#' @param exposure A character string with the exposure name
#' @param seed An integer defining a seed
#' @param sex A logical defining if a sex-stratififed analysis is run, default to false
#'
#' @import dplyr
#' @importFrom stats as.formula
#' @importFrom stats lm
#' @importFrom stats na.omit
#' @importFrom tibble column_to_rownames
#' @importFrom tibble rownames_to_column
#' @importFrom ccmm ccmm
#'
#' @return A data.frame containing a p value for the total mediation effect for each CpG

.Mediation <- function(CpG,
                       expo_cov,
                       mediator,
                       exposure,
                       seed, 
                       sex = FALSE) {
  
  # Keeping CpG as a matrix will preserve the row names (apply coerces it to numeric vector)
  CpG_mat <- as.matrix(CpG)
  
  # Change the CpG colname name to "y" as it will be used as such in the regression formula
  colnames(CpG_mat) <- "y"
  
  # Change ID to rownames so covariates can be merged with methylation data
  expo_cov_id <- expo_cov %>%
    tibble::column_to_rownames(var = "id")
  
  # Create a subset of data containing methylation data of one CpG (y) and exposure-covariates data (x)
  dataset <- merge(x = expo_cov_id,
                   y = CpG_mat,
                   by = "row.names") %>%
    stats::na.omit()
  
  if (isFALSE(sex)) {
    
    # mediation package does not work if the lm model contains expression 'formula' instead of variables in a form
    # of Y ~ X + C1 + C2... That is why generic function passing the names of variables to lm will not work
    # Yuan 2021: gestational age, sex, and ancestry were not significantly associated with cell composition.
    # But the sample size supporting these findings was small and future studies with more appropriate power are
    # needed to answer how much these factors play in contributing to placental cell composition variability
    Mreg <- stats::lm(pc1 ~ TCS_log2 +
                        mother_bmi +
                        child_sex +
                        gest_age +
                        mother_act_smoke +
                        mother_age +
                        parity +
                        mother_edu +
                        season_conc,
                      data = dataset)
    
    Yreg <- stats::lm(y ~ TCS_log2 +
                        mother_bmi +
                        child_sex +
                        gest_age +
                        mother_act_smoke +
                        mother_age +
                        parity +
                        mother_edu +
                        season_conc +
                        batch +
                        plate +
                        chip +
                        pc1,
                      data = dataset)
    
  } else {
    
    Mreg <- stats::lm(pc1 ~ TCS_log2 +
                        mother_bmi +
                        gest_age +
                        mother_act_smoke +
                        mother_age +
                        parity +
                        mother_edu +
                        season_conc,
                      data = dataset)
    
    Yreg <- stats::lm(y ~ TCS_log2 +
                        mother_bmi +
                        gest_age +
                        mother_act_smoke +
                        mother_age +
                        parity +
                        mother_edu +
                        season_conc +
                        batch +
                        plate +
                        chip +
                        pc1,
                      data = dataset)
  }
  
  set.seed(seed)
  
  # Run mediation analysis using results of linear regressions
  med_analysis <- mediation::mediate(Mreg,
                                     Yreg,
                                     treat = exposure,
                                     mediator = mediator,
                                     boot = FALSE)
  
  # Extract the p value for the average causal mediation effect (ACME)
  acme <- data.frame("acme_CIs" = stringr::str_c(round(med_analysis$d0, 3),
                                                 " (", 
                                                 round(med_analysis$d0.ci[1], 3),
                                                 "; ",
                                                 round(med_analysis$d0.ci[2], 3),
                                                 ")"),
                     "acme_p_value" = round(med_analysis$d0.p, 2))
  
  return(acme)
}


#' Main mediation analysis using ccmm package
#'
#' @param outcome A numeric vector containing methylation data from 1 CpG
#' @param expo_cov A data.frame containing all variables needed for regression (exposures, covariates, mediators)
#' @param mediators A character string defining mediators names
#' @param treatment A character string with the exposure name
#' @param covariates A character string with covariates names
#' @param seed An integer defining a seed
#' @param sex A logical defining if a sex-stratififed analysis is run, default to false
#'
#' @import dplyr
#' @import ccmm
#' @importFrom stats as.formula
#' @importFrom stats lm
#' @importFrom stats na.omit
#' @importFrom tibble column_to_rownames
#' @importFrom tibble rownames_to_column
#'
#' @return A data.frame containing a p value for the total mediation effect for each CpG

.MediationCCMM <- function(CpG,
                           expo_cov,
                           exposure,
                           confounders,
                           confounders_dum,
                           cells, 
                           seed) {
  
  options(scipen = 999)
  
  CpG <- as.matrix(CpG)
  
  meth_expo_cov_sel <- expo_cov %>% 
    dplyr::select(id, 
                  tidyselect::all_of(exposure), 
                  tidyselect::all_of(confounders), 
                  tidyselect::all_of(cells)) %>% 
    tibble::column_to_rownames("id") %>% 
    
    # While merging, the CpG name changes to "V1"
    merge(CpG, by = 0) %>% 
    fastDummies::dummy_cols(remove_first_dummy = TRUE,
                            remove_selected_columns = TRUE) %>% 
    stats::na.omit()
  
  outcome <- meth_expo_cov_sel[, "V1"]
  mediators <- as.matrix(meth_expo_cov_sel[, cells])
  treatment <- meth_expo_cov_sel[, exposure]
  covariates <- as.matrix(meth_expo_cov_sel[, confounders_dum])
  
  set.seed(seed)
  
  ccmm_analysis <- ccmm::ccmm(y = outcome, 
                              M = mediators, 
                              tr = treatment, 
                              x = covariates)
  
  ccmm <- data.frame("ccmm_CIs" = stringr::str_c(round(ccmm_analysis$TIDE, 3),
                                                 " (", 
                                                 round(ccmm_analysis$TIDE.CI[1], 3),
                                                 "; ",
                                                 round(ccmm_analysis$TIDE.CI[2], 3),
                                                 ")")) %>% 
    cbind(t(ccmm_analysis$IDEs))
  
  return(ccmm)
}



#' Parallelized mediation analysis
#'
#' @param meth_data A dataframe containing methylation data
#' @param expo_cov A data.frame containing all variables needed for regression (exposures and confounders)
#' @param mediator A string defining mediator name
#' @param ncores The number of cores used for parallel computing, by default all available cores
#' @param exposure A character string with the exposure name
#' @param seed An integer defining a seed
#' @param sex A logical defining if a sex-stratififed analysis is run, default to false
#'
#' @import parallel
#' @import doParallel
#' @import foreach
#'
#' @return A data.frame containing a p value for the total mediation effect for each CpG
#' @export
#'
MediationAnalysisSepages <- function(meth_data,
                                     expo_cov,
                                     mediator = NULL,
                                     ncores = bigstatsr::nb_cores(),
                                     exposure,
                                     confounders = NULL,
                                     confounders_dum = NULL,
                                     cells = NULL,
                                     seed = 111666, 
                                     sex = FALSE,
                                     ccmm = FALSE) {
  
  # Set the parallel computing conditions
  
  # Create a set of copies of R running in parallel and communicating over sockets
  cl <- parallel::makeCluster(getOption("cl.cores", ncores))
  
  # Export used functions and objects to be defined before creation of the clusters
  parallel::clusterEvalQ(cl, c(library(magrittr)))
  
  # Register the parallel backend
  doParallel::registerDoParallel(cl)
  
  # Explicitly register a sequential parallel backend. This will prevent a warning message
  # from being issued if the %dopar% function is called and no parallel backend has been registered.
  foreach::registerDoSEQ()
  
  # Stop the cluster
  on.exit(parallel::stopCluster(cl))
  
  if (isFALSE(ccmm)) {
    
    res_med <- parallel::parApply(cl = cl,
                                  X = meth_data,
                                  MARGIN = 2,
                                  FUN = .Mediation,
                                  expo_cov = expo_cov,
                                  mediator = mediator,
                                  exposure = exposure,
                                  seed,
                                  sex)
    
  } else {
    
    res_med <- parallel::parApply(cl = cl,
                                  X = meth_data,
                                  MARGIN = 2,
                                  FUN = .MediationCCMM,
                                  expo_cov = expo_cov,
                                  confounders = confounders,
                                  exposure = exposure,
                                  confounders_dum = confounders_dum,
                                  cells = cells,
                                  seed)
  }
  
  
  tidy_res <- do.call(rbind, res_med) %>% 
    tibble::rownames_to_column("CpG")
  
  if (isTRUE(ccmm)) {
    
    colnames(tidy_res)[-c(1:2)] <- cells
    tidy_res <- dplyr::mutate_if(tidy_res, is.numeric, round, 5)
  }
  
  return(tidy_res)
}
