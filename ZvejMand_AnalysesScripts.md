---
title: "Mandibular shape in the Meso-Neolithic transition: the Zvejnieki case study"
author: "Maria Ana Correia"
date: '02/09/2025'
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


``` r
#this makes images save in folder in directory
knitr::opts_chunk$set(
  echo = TRUE, #shows code
  warning = FALSE, message = FALSE, #stops warning messages
  fig.path = "images/",
  dev = c("svg", "png", "tiff"), #saves figures as svg, tiff, and png in images folder
  dpi = 500, #publishing quality for combination art (Elsevier)
  tidy.opts=list(width.cutoff=60), # stops code from running off page
  tidy=TRUE
)

# function to load or install from CRAN or GitHub
# installs missing packages when sharing with others
load_or_install <- function(pkg, github = NULL) {
  # If package not available, install it
  if (!requireNamespace(pkg, quietly = TRUE)) {
    if (is.null(github)) {
      # Install from CRAN
      install.packages(pkg, dependencies = TRUE)
    } else {
      # Install from GitHub (needs remotes package)
      if (!requireNamespace("remotes", quietly = TRUE)) {
        install.packages("remotes")
      }
      remotes::install_github(github, dependencies = TRUE)
    }
  }
  # Finally load it
  suppressPackageStartupMessages(
    library(pkg, character.only = TRUE)
  )
}
```


``` r
# CRAN packages
cran_pkgs <- c(
  "tidyverse", #everyday data analyses
  "styler", #source code formatter
  "arrow", #cross-language development platform to export .parquet
  "usethis", #automates repetitive tasks that arise during project setup
   "osfr", #interface for OSF
  "geomorph" #geometric morphometrics
  )

purrr::walk(cran_pkgs, load_or_install)

# GitHub packages (supply the repo as github = "user/repo")
load_or_install("SlicerMorphR", github = "SlicerMorph/SlicerMorphR")
#import SlicerMorph dataset into R
```


``` r
#!setting decimals
fmt_decimals <- function(decimals = 0) {
  function(x) format(x, nsmall = decimals, scientific = FALSE)
}
#graphical settings for ggplot
my_theme <- theme(
  axis.text = element_text(size = 8, colour = "black"),
  # makes numbers smaller and black (consider final display)
  axis.ticks = element_line(linewidth = 0.5,
                            colour = "black"),
  # same for ticks
  axis.title = element_text(size = 10),
  # and for axis titles
  panel.grid.minor = element_blank(),
  panel.background = element_blank(),
  panel.grid.major = element_blank(),
  panel.border = element_rect(
    colour = "black", fill = NA, size = 0.5)
)

#calculate outliers
is_outlier <- function(x) {
  return(
    x < quantile(x, 0.25, na.rm = TRUE) - 1.5 * IQR(x, na.rm = TRUE) |
      x > quantile(x, 0.75, na.rm = TRUE) + 1.5 * IQR(x, na.rm = TRUE)
  )
}
```

Fragmented specimens were virtually pieced together in 3D slicer using the Fiducial Registration Wizard [@godinho2020]. Then, in the Markups Module of 3DSlicer, coordinates were extracted from a total of 21 anatomical landmarks from the most complete hemi-mandible of each specimen to capture mandibular morphology [@godinho2022]. The use of left hemimandibles was favoured, but because that was also the favoured side when sampling, the sample ended up even.

# Preparing Landmark Data for `geomorph`

We implemented two parallel workflows for preparing landmark data into the format required by the `geomorph` package. The goal in both cases is to produce a 3D numeric array with dimensions *(p landmarks × k coordinates × n specimens)*, suitable for downstream analyses.

1.  **From `.mark.json` files**\
    Individual specimen landmark files produced by PAT were parsed from JSON into matrices of coordinates. These were then stacked into a 3D array, with each slice representing one specimen.

2.  **From `.csv` files**\
    Landmark coordinates saved as CSVs were read into matrices and similarly stacked into a 3D array with the same dimensional structure.

Both approaches produce the object `array3d` for use in `geomorph`, and use depends on whether user has access to raw (`.mark.json`) or derived (`.csv`) data on [OSF](https://osf.io/vkat9/) [@foster2017][^1]. Some safety checks were also included to ensure the produced dataset adheres to `geomorph` criteria.

[^1]: in the package `osfr`, PATs are required to upload files, create projects/components, access information about your private projects, or download files in your private projects. PATs are not required for accessing information about public projects or downloading public files, but authentication with a PAT will increase the rate limit on the API


``` r
# --- Dual-mode OSF pipeline for landmarks ---
# Project vkat9 | Raw Data: kwafd | Derived Data: 9fnsp | Analyses: a75wu

# --- Directories ---
raw_dir <- "data/raw"
derived_dir <- "data/derived"
usethis::use_directory(raw_dir)
usethis::use_directory(derived_dir)

# --- Check PAT ---
if (nzchar(Sys.getenv("OSF_PAT"))) {
  message("Authenticated with OSF PAT → building array from JSON")

  # --- PAT branch ---
  # 1) Retrieve JSON files from Raw Data component
  json_files <- osf_retrieve_node("kwafd") |>
    osf_ls_files() |>
    dplyr::filter(str_detect(name, "\\.json$"))

  # 2) Download JSON files to raw_dir
  osf_download(json_files, path = raw_dir, conflicts = "overwrite")

  # 3) Specimen IDs from filenames
  json_paths <- list.files(raw_dir, pattern = "\\.json$", full.names = TRUE)
  specimen_ids <- tools::file_path_sans_ext(basename(json_paths))

  # 4) Read JSON into numeric matrices
  read_lmk_matrix <- function(path) {
    m <- SlicerMorphR::read.markups.json(path)
    m <- as.matrix(m)
    m <- apply(m, 2, as.numeric)
    if (is.null(colnames(m))) colnames(m) <- c("X", "Y", "Z")
    m
  }
  landmark_list <- map(json_paths, read_lmk_matrix) |>
    set_names(specimen_ids)

  # 5) Build array3d
  p <- nrow(landmark_list[[1]])
  k <- 3
  n <- length(landmark_list)

  array3d <- array(
    NA_real_,
    dim = c(p, k, n),
    dimnames = list(
      landmark = seq_len(p),
      coord    = c("X", "Y", "Z"),
      specimen = names(landmark_list)
    )
  )
  for (i in seq_along(landmark_list)) array3d[, , i] <- landmark_list[[i]]
  storage.mode(array3d) <- "double"

  # 6) Save array3d locally as RDS
  saveRDS(array3d, file.path(derived_dir, "array3d.RDS"))

  # 7) Create CSV + Parquet from array3d
  df2d <- data.frame(specimen = dimnames(array3d)$specimen)
  for (i in 1:p) {
    df2d[[paste0("X", i)]] <- array3d[i, 1, ]
    df2d[[paste0("Y", i)]] <- array3d[i, 2, ]
    df2d[[paste0("Z", i)]] <- array3d[i, 3, ]
  }

  # 8) Upload CSV + Parquet to OSF Derived Data
  derived_node <- osf_retrieve_node("9fnsp")
  tmp_csv <- file.path(tempdir(), "landmarks.csv")
  tmp_parquet <- file.path(tempdir(), "landmarks.parquet")

  write_csv(df2d, tmp_csv)
  write_parquet(df2d, tmp_parquet)

  osf_upload(derived_node, path = tmp_csv, conflicts = "overwrite")
  osf_upload(derived_node, path = tmp_parquet, conflicts = "overwrite")
} else {
  # --- No PAT branch ---
  message("No OSF PAT → rebuilding array from Derived Data CSV")

  # 1) Download CSV from Derived Data
  csv_file <- derived_node |>
    osf_ls_files() |>
    dplyr::filter(name == "landmarks.csv")

  osf_download(csv_file, path = derived_dir, conflicts = "error")

  # 2) Rebuild array3d from CSV
  df2d <- read_csv(file.path(derived_dir, "landmarks.csv"))
  df_matrix <- as.matrix(df2d %>% select(-specimen))

  p <- ncol(df_matrix) / 3
  k <- 3
  n <- nrow(df_matrix)

  array3d <- geomorph::arrayspecs(df_matrix, p = p, k = k)

  # 3) Set dimnames identical to PAT branch
  dimnames(array3d)[[1]] <- seq_len(p) # landmarks
  dimnames(array3d)[[2]] <- c("X", "Y", "Z") # coords
  dimnames(array3d)[[3]] <- df2d$specimen # specimen names
}
```

<div data-pagedtable="false">
  <script data-pagedtable-source type="application/json">
{"columns":[{"label":["name"],"name":[1],"type":["chr"],"align":["left"]},{"label":["id"],"name":[2],"type":["chr"],"align":["left"]},{"label":["meta"],"name":[3],"type":["list"],"align":["right"]}],"data":[{"1":"landmarks.parquet","2":"68b4a0ea3560139ce496af8a","3":"<named list [3]>"}],"options":{"columns":{"min":{},"max":[10]},"rows":{"min":[10],"max":[10]},"pages":{}}}
  </script>
</div>

``` r
# --- Safety checks ---
# Ensure array has correct structure
is.numeric(array3d)
```

```
## [1] TRUE
```

``` r
length(dim(array3d)) == 3
```

```
## [1] TRUE
```

``` r
dim(array3d)  # should be (landmarks, 3, specimens)
```

```
## [1] 21  3 10
```

``` r
# Confirm specimen IDs are in the 3rd dimension
dimnames(array3d)[[3]]
```

```
##  [1] "ILH_Zvejnieki_34.115.228.mrk"     "ILH_Zvejnieki_34.118.234.mrk"    
##  [3] "ILH_Zvejnieki_34.22.63_miss.mrk"  "ILH_Zvejnieki_34.3.16.mrk"       
##  [5] "ILH_Zvejnieki_34.43.112_miss.mrk" "ILH_Zvejnieki_34.46.115.mrk"     
##  [7] "ILH_Zvejnieki_34.58.134.mrk"      "ILH_Zvejnieki_34.66.153.mrk"     
##  [9] "ILH_Zvejnieki_34.74.164.mrk"      "ILH_Zvejnieki_34.85.181.mrk"
```

``` r
# Check requirement: number of specimens > number of coordinates
dim(array3d)[3] > dim(array3d)[2]
```

```
## [1] TRUE
```

When specimens were incomplete, the location of the missing LMs was estimated using the thin plate spline (TPS) function (`estimate.missing()`) of the `geomorph` package.


``` r
#estimate.missing() WARNING: n > k
estimateLMs <- estimate.missing(array3d, method = "TPS")
```


``` r
# 1) GPA with fixed landmarks (ProcD = TRUE to compute/return Procrustes distances)
gpa <- geomorph::gpagen(A = estimateLMs, ProcD = TRUE)
```

```
## 
## Performing GPA
##   |                                                                              |                                                                      |   0%  |                                                                              |==================                                                    |  25%  |                                                                              |===================================                                   |  50%  |                                                                              |======================================================================| 100%
## 
## Making projections... Finished!
```

``` r
# 2) PCA on aligned coordinates
pca <- geomorph::gm.prcomp(gpa$coords)

# 3) Capture existing contents before Slicer export
pre_export <- dir(derived_dir, full.names = TRUE, all.files = FALSE)

# 4) Export SlicerMorph bundle to derived_dir
# Note: argument renamed to 'output_folder' per your current SlicerMorphR version
SlicerMorphR::geomorph2slicermorph2(
  gpa = gpa,
  pca = pca,
  output.folder = derived_dir
)

# 5) Determine which files/folders were newly created by the export
post_export <- dir(derived_dir, full.names = TRUE, all.files = FALSE)
slicer_exports <- setdiff(post_export, pre_export)

# 6) --- Check PAT ---
if (nzchar(Sys.getenv("OSF_PAT"))) {
  cat("Authenticated with OSF PAT\n")

  # Upload all non-RDS artifacts from derived_dir (recursively), overwriting on OSF
  files_to_upload <- list.files(
    derived_dir,
    recursive = TRUE,
    full.names = TRUE,
    include.dirs = FALSE
  )
  files_to_upload <- files_to_upload[
    !grepl("\\.rds$", files_to_upload, ignore.case = TRUE)
  ]

  if (length(files_to_upload) == 0) {
    cat("No non-RDS files found to upload from derived_dir.\n")
  } else {
    cat("Uploading ", length(files_to_upload), " file(s) to OSF (overwrite on conflict)...\n", sep = "")
    for (f in files_to_upload) {
      osfr::osf_upload(derived_node, path = f, conflicts = "overwrite")
    }
    cat("Upload complete to OSF node 9fnsp.\n")
  }
} else {
  # --- No PAT branch ---
  cat("No OSF PAT. Outputs only saved locally in: ", derived_dir, "\n", sep = "")
}
```

```
## Authenticated with OSF PAT
## Uploading 6 file(s) to OSF (overwrite on conflict)...
## Upload complete to OSF node 9fnsp.
```

TODO: upload to SlicerMorph again to see errors in placing landmarks and see if I can see shapes (video); check out on paper stuff to produce graphs on ggplot using `make_ggplot`; code not displaying on html and removed derived data folder from github



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
