---
title: "Mandibular shape in the Meso-Neolithic transition: the Zvejnieki case study"
author: "Maria Ana Correia"
date: '01/09/2025'
subtitle: "Analyses Scripts"
output:
  html_document:
    code_folding: hide
    code_link: true
    keep_md: true
    toc: true
    toc_float: true
    df_print: paged
  pdf_document:
    latex_engine: xelatex
    citation_package: default
    includes:
      in_header: mypackages.tex
bibliography: references.bib
link-citations: true
csl: apa-single-spaced
---

This file documents the geometric morphometrics of Zvejnieki mandibles. 3D landmarks were collected using 3DSlicer and statistical analyses used geomorph and SlicerMorph/SllicerMorphR [@adams2025; @rolfe2021].









TODO: figure out missinglandmarks() TRY to add more jsonsstr





# References

::: {#refs}
:::


```
## R version 4.5.1 (2025-06-13 ucrt)
## Platform: x86_64-w64-mingw32/x64
## Running under: Windows 11 x64 (build 26100)
## 
## Matrix products: default
##   LAPACK version 3.12.1
## 
## locale:
## [1] LC_COLLATE=Portuguese_Portugal.utf8  LC_CTYPE=Portuguese_Portugal.utf8   
## [3] LC_MONETARY=Portuguese_Portugal.utf8 LC_NUMERIC=C                        
## [5] LC_TIME=Portuguese_Portugal.utf8    
## 
## time zone: Europe/Lisbon
## tzcode source: internal
## 
## attached base packages:
## [1] stats     graphics  grDevices utils     datasets  methods   base     
## 
## other attached packages:
##  [1] SlicerMorphR_0.0.2.0000 jsonlite_2.0.0          geomorph_4.0.10        
##  [4] Matrix_1.7-3            rgl_1.3.24              RRPP_2.1.2             
##  [7] osfr_0.2.9              usethis_3.2.0           arrow_21.0.0.1         
## [10] styler_1.10.3           lubridate_1.9.4         forcats_1.0.0          
## [13] stringr_1.5.1           dplyr_1.1.4             purrr_1.1.0            
## [16] readr_2.1.5             tidyr_1.3.1             tibble_3.3.0           
## [19] ggplot2_3.5.2           tidyverse_2.0.0        
## 
## loaded via a namespace (and not attached):
##  [1] gtable_0.3.6       xfun_0.53          bslib_0.9.0        htmlwidgets_1.6.4 
##  [5] lattice_0.22-7     tzdb_0.5.0         vctrs_0.6.5        tools_4.5.1       
##  [9] generics_0.1.4     curl_7.0.0         parallel_4.5.1     pkgconfig_2.0.3   
## [13] R.oo_1.27.1        RColorBrewer_1.1-3 assertthat_0.2.1   lifecycle_1.0.4   
## [17] R.cache_0.17.0     compiler_4.5.1     farver_2.1.2       htmltools_0.5.8.1 
## [21] sass_0.4.10        yaml_2.3.10        crayon_1.5.3       pillar_1.11.0     
## [25] jquerylib_0.1.4    R.utils_2.13.0     cachem_1.1.0       mime_0.13         
## [29] nlme_3.1-168       tidyselect_1.2.1   digest_0.6.37      stringi_1.8.7     
## [33] rprojroot_2.1.1    fastmap_1.2.0      grid_4.5.1         cli_3.6.5         
## [37] magrittr_2.0.3     base64enc_0.1-3    triebeard_0.4.1    crul_1.6.0        
## [41] ape_5.8-1          withr_3.0.2        scales_1.4.0       bit64_4.6.0-1     
## [45] timechange_0.3.0   rmarkdown_2.29     httr_1.4.7         jpeg_0.1-11       
## [49] bit_4.6.0          R.methodsS3_1.8.2  hms_1.1.3          memoise_2.0.1     
## [53] evaluate_1.0.5     knitr_1.50         urltools_1.7.3.1   rlang_1.1.6       
## [57] Rcpp_1.1.0         glue_1.8.0         httpcode_0.3.0     vroom_1.6.5       
## [61] rstudioapi_0.17.1  R6_2.6.1           fs_1.6.6
```
