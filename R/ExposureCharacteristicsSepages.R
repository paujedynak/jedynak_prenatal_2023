# quiets concerns of R CMD check re: the .'s that appear in pipelines
if (getRversion() >= "2.15.1") {
  utils::globalVariables(c("."))
}

# quiets concerns of R CMD check when variables appear in pipelines
utils::globalVariables(c("mean_expo", "t2", "t3"))


#' Calculates exposure characteristics
#'
#' @param input_data A data.frame containing information on exposure concentrations and <LOD
#' @param lod_values A double defining LOD value for a compound
#' @param path A string defining the path to save the output
#' @param file_name A string defining specific part of the file name
#'
#' @return A data.frame with exposure descriptive statistics
#' @export
#'
#' @import dplyr
#' @import tidyr
#' @importFrom naniar replace_with_na_all
#' @importFrom stats quantile
#' @importFrom utils write.csv
#'
ExposureCharacteristicsSepages <- function(input_data,
                                           lod_values,
                                           path,
                                           file_name) {
  
  expo_summary <- input_data %>%
    dplyr::mutate(LOD = lod_values) %>% 
    tidyr::pivot_longer(cols = c("t2", "t3", "mean_expo"),
                        names_to = "period",
                        values_to = "val") %>% 
    dplyr::group_by(period) %>% 
    
    # Calculate percentage of exposure samples >LOD
    dplyr::summarise(n = sum(!is.na(val)),
                     "pct_det" = round((sum(val >= LOD, na.rm = TRUE) * 100) / n, 1),
                     "p5" = round(stats::quantile(val, na.rm = TRUE, probs = 0.05), 1),
                     "median" = round(stats::quantile(val, na.rm = TRUE, probs = 0.5), 1),
                     "p95" = round(stats::quantile(val, na.rm = TRUE, probs = 0.95), 1)
    ) %>% 
    
    tidyr::pivot_wider(
      names_from = period,
      values_from = n:p95) %>%
    dplyr::select(contains("t2"), contains("t3"), dplyr::everything(), -pct_det_mean_expo) %>% 
  
  dplyr::mutate(LOD = lod_values) %>% 
    dplyr::select(LOD, dplyr::everything())
  
  write.csv(expo_summary, here::here(path, paste0(file_name, ".csv")), row.names = FALSE)
  
  return(expo_summary)
}
