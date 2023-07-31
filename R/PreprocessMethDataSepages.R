#' Prepare methylation file before further processing (remove cross-reactive and sex probes), add IDs
#'
#' @param Bvalues A matrix of BMIQ data provided by Johanna, contains duplicates
#' @param path A string defining path for saving the methylation design file
#' @param file_name A string defining file name for saving
#'
#' @return
#' @export

PreprocessMethDataSepages <- function(Bvalues,
                                      path, 
                                      file_name) {
  
  # ################################################################################
  # # 1. Removal of CpGs close to SNPs (dist<2 bp) and ChrXY
  # #------------------------------------------------------------------------------#
  # 
  # Bvalues <- readRDS(bvalFile)
  
  Bvalues <- DMRcate::rmSNPandCH(Bvalues, 
                                 dist = 2, 
                                 mafcut = 0.05,
                                 rmcrosshyb = FALSE, 
                                 rmXY = TRUE)
  # dim(Bvalues)
  
  # 784 494 CpGs 
  
  ################################################################################
  # 3. Removal of cross-reactive probes
  #------------------------------------------------------------------------------#
  
  xloci <- maxprobes::xreactive_probes(array_type = "EPIC")
  xloci <- unlist(unique(xloci))
  # length(xloci)
  # 87 464 known cross-reactive probes
  
  x <- which(rownames(Bvalues) %in% xloci)
  # length(x)
  # 31 917 known cross-reactive probes still in the dataset
  # round(length(x)/nrow(Bvalues),2)
  
  Bvalues <- Bvalues[-x,]
  # dim(Bvalues)
  # 752 577  CpGs 
  
  # ################################################################################
  # # Save
  # #------------------------------------------------------------------------------#
  # 
  saveRDS(Bvalues, here::here(path, file_name))
  
  return(Bvalues) 
  
}