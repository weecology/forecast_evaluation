source("functions.R")

load_dependencies()
data_set <- prep_data_set()
models <- run_models(data_set, run = FALSE)
make_figures(data_set, models)


# sessionInfo(): 
# 
# R version 4.0.3 (2020-10-10)
# Platform: x86_64-w64-mingw32/x64 (64-bit)
# Running under: Windows 10 x64 (build 19041)
# 
# Matrix products: default
# 
# locale:
#  LC_COLLATE=English_United States.1252 
#  LC_CTYPE=English_United States.1252   
#  LC_MONETARY=English_United States.1252
#  LC_NUMERIC=C                          
#  LC_TIME=English_United States.1252    
# 
# attached base packages:
#  stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
#  e1071_1.7-4        viridis_0.5.1      viridisLite_0.3.0  scoringRules_1.0.1
#  runjags_2.0.4-6    portalr_0.3.6      png_0.1-7          coda_0.19-4       
# 
# loaded via a namespace (and not attached):
#  Rcpp_1.0.6        knitr_1.25        magrittr_2.0.1    MASS_7.3-53      
#  tidyselect_1.1.0  munsell_0.5.0     colorspace_2.0-0  lattice_0.20-41  
#  R6_2.5.0          rlang_0.4.10      dplyr_1.0.3       tools_4.0.3      
#  parallel_4.0.3    grid_4.0.3        gtable_0.3.0      xfun_0.10        
#  class_7.3-17      ellipsis_0.3.1    tibble_3.0.5      lifecycle_0.2.0  
#  crayon_1.3.4      gridExtra_2.3     tidyr_1.1.2       purrr_0.3.4      
#  ggplot2_3.3.3     vctrs_0.3.6       clisymbols_1.2.0  glue_1.4.2       
#  pillar_1.4.7      compiler_4.0.3    generics_0.1.0    scales_1.1.1     
#  lubridate_1.7.9.2 pkgconfig_2.0.3 
# 