resid_check <- function(CpG,
                        meth_data,
                        expo_cov,
                        exposure,
                        confounders,
                        technical_confounders) {
  
  # Keeping CpG as a matrix will preserve the row names (apply coerces it to numeric vector)
  CpG_data <- as.matrix(meth_data[, CpG])
  
  # Change the CpG colname name to "y" as it will be used as such in the regression formula
  
  colnames(CpG_data) <- "y"
  
  # Change ID to rownames so covariates can be merged with methylation dataset
  expo_cov_id <- expo_cov %>%
    column_to_rownames(var = "id")
  
  # Create a subset of data containing methylation data of one CpG (y) and exposure-covariates data (x)
  data <- merge(x = expo_cov_id, y = CpG_data, by = "row.names")
  
  formula <-
    as.formula(paste(
      "y ~",
      exposure,
      "+",
      paste(confounders, collapse = " + "),
      "+",
      paste(technical_confounders, collapse = " + ")
    ))
  
  fit <- rlm(formula = formula,
             data = data,
             maxit = 400)
  
  plot <- plot(fitted(fit), 
               resid(fit),
               main = CpG)
  abline(0, 0)
  
  return(plot)
}