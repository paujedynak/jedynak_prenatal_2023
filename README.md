# Prenatal exposure to triclosan assessed in multiple urine samples and placental DNA methylation

## Article info

Repository contains reproducible analyses for the paper: [Jedynak et al., 2023](https://www.sciencedirect.com/science/article/pii/S0269749123011995).

## Analysis overview

Analyses depend on R (>= 4.3.0).

Code used to produce the files is found in the corresponding file: `jedynak_prenatal_2023.Rmd` (executable file).

This analysis was performed under Windows 11 x64 (build 22621) using:    
* [R 4.3.1](https://cran.r-project.org/bin/windows/base) (2023-06-16 ucrt)    
* [renv 1.0.0](https://cran.r-project.org/web/packages/renv/index.html) dependency management package


### Packages

All packages used in the analyses are saved in the `renv/library/R-4.3/x86_64-w64-mingw32` folder at the version used to produce these results, under control of the `renv` package manager. Re-running the analysis requires execution of `renv::restore()` which will upgrade/ downgrade user's packages to the versions used in the present study. This operation will modify the packages only locally (for the project), so it will not affect user's package library. All "in-house" packages are available in the R/ folder.


### To run

Re-running the analysis requires an additional `data/raw_data` folder that is not shared here. These data can only be provided upon request and after approval by the SEPAGES consortium (contact: [Sarah Lyon-Caen](sarah.lyon-caen@univ-grenoble-alpes.fr)). Running the script: `jedynak_prenatal_2023.Rmd` will allow to fully reproduce the analyses, figures and tables.


## Repo organization

### data/raw_data/ folder

Analysis input data-files are not made available as they contain sensitive information. The analysis input data files would be:

* `data_sepages_211022.sas7bdat` = covariates data, including potential confounders, n = 484
* `phenols_phthalates_pregnancy_mr_2020-03-26_448.csv` = phenols exposure data, n = 479
* `bmiq_processed.RDS` = methylation data (no replicates, 814,481 CpGs, not filtered, normalized using BMIQ method), n = 395
* `SEPAGES_SampleSheet_nodup_20211208.csv` = DNA methylation measurement technical factors (batch, plate, chip), n = 395
* `SEPAGES_CC.planet_rpc.csv` = placental cell mix proportions
* `data_nsample_230314.Rdata` = data on number of urine samples per pool

Additional file `EWAS_rep_nl_FDR.RDS` is not available in the repo due to its large size.


### R/ folder

This folder contains all the in-house functions and packages used for the analyses.


## Session info

```
R version 4.3.1 (2023-06-16 ucrt)
Platform: x86_64-w64-mingw32/x64 (64-bit)
Running under: Windows 11 x64 (build 22621)

Matrix products: default


locale:
[1] LC_COLLATE=English_Europe.utf8  LC_CTYPE=English_Europe.utf8   
[3] LC_MONETARY=English_Europe.utf8 LC_NUMERIC=C                   
[5] LC_TIME=English_Europe.utf8    

time zone: Europe/Madrid
tzcode source: internal

attached base packages:
[1] parallel  stats4    stats     graphics  grDevices datasets  utils    
[8] methods   base     

other attached packages:
 [1] lubridate_1.9.2                                    
 [2] forcats_1.0.0                                      
 [3] stringr_1.5.0                                      
 [4] dplyr_1.1.2                                        
 [5] purrr_1.0.1                                        
 [6] readr_2.1.4                                        
 [7] tidyr_1.3.0                                        
 [8] tibble_3.2.1                                       
 [9] tidyverse_2.0.0                                    
[10] scales_1.2.1                                       
[11] RobustRegressions_1.5.1                            
[12] robCompositions_2.3.1                              
[13] data.table_1.14.8                                  
[14] pls_2.8-2                                          
[15] rio_0.5.29                                         
[16] Plots_1.3                                          
[17] IlluminaHumanMethylationEPICanno.ilm10b4.hg19_0.6.0
[18] minfi_1.46.0                                       
[19] bumphunter_1.42.0                                  
[20] locfit_1.5-9.8                                     
[21] iterators_1.0.14                                   
[22] foreach_1.5.2                                      
[23] Biostrings_2.68.1                                  
[24] XVector_0.40.0                                     
[25] SummarizedExperiment_1.30.2                        
[26] Biobase_2.60.0                                     
[27] MatrixGenerics_1.12.2                              
[28] matrixStats_1.0.0                                  
[29] GenomicRanges_1.52.0                               
[30] GenomeInfoDb_1.36.1                                
[31] IRanges_2.34.1                                     
[32] S4Vectors_0.38.1                                   
[33] BiocGenerics_0.46.0                                
[34] magrittr_2.0.3                                     
[35] janitor_2.2.0                                      
[36] Hmisc_5.1-0                                        
[37] here_1.0.1                                         
[38] Helpers_0.1.0                                      
[39] ggrepel_0.9.3                                      
[40] ggpubr_0.6.0                                       
[41] DMR_0.1.3                                          
[42] ccmm_1.0                                           
[43] bacon_1.28.0                                       
[44] ellipse_0.5.0                                      
[45] BiocParallel_1.34.2                                
[46] ggplot2_3.4.2                                      

loaded via a namespace (and not attached):
  [1] dichromat_2.0-0.1         progress_1.2.2           
  [3] pacman_0.5.1              nnet_7.3-19              
  [5] HDF5Array_1.28.1          TH.data_1.1-2            
  [7] vctrs_0.6.3               digest_0.6.33            
  [9] png_0.1-8                 proxy_0.4-27             
 [11] pcaPP_2.0-3               corrplot_0.92            
 [13] deldir_1.0-9              renv_1.0.0               
 [15] hdrcde_3.4                MASS_7.3-60              
 [17] reshape_0.8.9             withr_2.5.0              
 [19] psych_2.3.6               xfun_0.39                
 [21] survival_3.5-5            doRNG_1.8.6              
 [23] memoise_2.0.1             diptest_0.76-0           
 [25] MatrixModels_0.5-2        zoo_1.8-12               
 [27] DEoptimR_1.1-0            Formula_1.2-5            
 [29] prettyunits_1.1.1         GGally_2.1.2             
 [31] prabclus_2.3-2            KEGGREST_1.40.0          
 [33] httr_1.4.6                bigparallelr_0.3.2       
 [35] rstatix_0.7.2             restfulr_0.0.15          
 [37] hash_2.2.6.2              cvTools_0.3.2            
 [39] rhdf5filters_1.12.1       fpc_2.2-10               
 [41] rhdf5_2.44.0              rstudioapi_0.15.0        
 [43] generics_0.1.3            base64enc_0.1-3          
 [45] curl_5.0.1                sfsmisc_1.1-15           
 [47] mitools_2.4               zlibbioc_1.46.0          
 [49] GenomeInfoDbData_1.2.10   quadprog_1.5-8           
 [51] rms_6.7-0                 doParallel_1.0.17        
 [53] xtable_1.8-4              pracma_2.4.2             
 [55] evaluate_0.21             S4Arrays_1.0.5           
 [57] BiocFileCache_2.8.0       preprocessCore_1.62.1    
 [59] hms_1.1.3                 colorspace_2.1-0         
 [61] filelock_1.0.2            readxl_1.4.3             
 [63] flexmix_2.3-19            lmtest_0.9-40            
 [65] snakecase_0.11.0          modeltools_0.2-23        
 [67] lattice_0.21-8            genefilter_1.82.1        
 [69] robustbase_0.99-0         SparseM_1.81             
 [71] cowplot_1.1.1             survey_4.2-1             
 [73] XML_3.99-0.14             class_7.3-22             
 [75] pillar_1.9.0              nlme_3.1-162             
 [77] compiler_4.3.1            stringi_1.7.12           
 [79] GenomicAlignments_1.36.0  plyr_1.8.8               
 [81] fda_6.1.4                 crayon_1.5.2             
 [83] abind_1.4-5               BiocIO_1.10.0            
 [85] truncnorm_1.0-9           flock_0.7                
 [87] haven_2.5.3               sp_2.0-0                 
 [89] bit_4.0.5                 sandwich_3.0-2           
 [91] multcomp_1.4-25           bigassertr_0.1.6         
 [93] codetools_0.2-19          perry_0.3.1              
 [95] openssl_2.1.0             ggfortify_0.4.16         
 [97] e1071_1.7-13              biovizBase_1.48.0        
 [99] fds_1.8                   multtest_2.56.0          
[101] splines_4.3.1             Rcpp_1.0.11              
[103] fastDummies_1.7.3         quantreg_5.96            
[105] dbplyr_2.3.3              sparseMatrixStats_1.12.2 
[107] cellranger_1.1.0          interp_1.1-4             
[109] knitr_1.43                blob_1.2.4               
[111] utf8_1.2.3                AnnotationFilter_1.24.0  
[113] checkmate_2.2.0           DelayedMatrixStats_1.22.1
[115] Gviz_1.44.0               openxlsx_4.2.5.2         
[117] ggsignif_0.6.4            Matrix_1.6-0             
[119] tzdb_0.4.0                pkgconfig_2.0.3          
[121] tools_4.3.1               cachem_1.0.8             
[123] RSQLite_2.3.1             NADA_1.6-1.1             
[125] DBI_1.1.3                 fastmap_1.1.1            
[127] rmarkdown_2.23            grid_4.3.1               
[129] Rsamtools_2.16.0          broom_1.0.5              
[131] xlsx_0.6.5                BiocManager_1.30.21.1    
[133] VariantAnnotation_1.46.0  carData_3.0-5            
[135] scrime_1.3.5              rpart_4.1.19             
[137] yaml_2.3.7                deSolve_1.36             
[139] latticeExtra_0.6-30       foreign_0.8-84           
[141] rtracklayer_1.60.0        illuminaio_0.42.0        
[143] cli_3.6.1                 siggenes_1.74.0          
[145] GEOquery_2.68.0           lifecycle_1.0.3          
[147] askpass_1.1               rainbow_3.7              
[149] mvtnorm_1.2-2             kernlab_0.9-32           
[151] backports_1.4.1           annotate_1.78.0          
[153] timechange_0.2.0          gtable_0.3.3             
[155] rjson_0.2.21              ggridges_0.5.4           
[157] limma_3.56.2              bitops_1.0-7             
[159] bit64_4.0.5               base64_2.0.1             
[161] zip_2.3.0                 ranger_0.15.1            
[163] polspline_1.1.23          rrcov_1.7-4              
[165] coMET_1.32.0              lazyeval_0.2.2           
[167] bigstatsr_1.5.12          zCompositions_1.4.0-1    
[169] htmltools_0.5.5           rJava_1.0-6              
[171] rappdirs_0.3.3            ensembldb_2.24.0         
[173] glue_1.6.2                VIM_6.2.2                
[175] RCurl_1.98-1.12           rprojroot_2.0.3          
[177] mclust_6.0.0              ks_1.14.0                
[179] mnormt_2.1.1              BSgenome_1.68.0          
[181] jpeg_0.1-10               gridExtra_2.3            
[183] boot_1.3-28.1             R6_2.5.1                 
[185] vcd_1.4-11                xlsxjars_0.6.1           
[187] GenomicFeatures_1.52.1    cluster_2.1.4            
[189] rngtools_1.5.2            Rhdf5lib_1.22.0          
[191] beanplot_1.3.1            robustHD_0.7.4           
[193] DelayedArray_0.26.6       tidyselect_1.2.0         
[195] ProtGenerics_1.32.0       htmlTable_2.4.1          
[197] xml2_1.3.5                QCEWAS_1.2-3             
[199] car_3.1-2                 AnnotationDbi_1.62.2     
[201] munsell_0.5.0             KernSmooth_2.23-22       
[203] laeken_0.5.2              nor1mix_1.3-0            
[205] htmlwidgets_1.6.2         RColorBrewer_1.1-3       
[207] biomaRt_2.56.1            rlang_1.1.1              
[209] fansi_1.0.4
```
