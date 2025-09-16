---
title: "Mandibular shape in the Meso-Neolithic transition: the Zvejnieki case study"
subtitle: "Analyses Scripts"
author: "Maria Ana Correia"
date: "15/09/2025"
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
bibliography: references.bib
link-citations: true
csl: "apa-single-spaced"
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
  "tidyverse", # everyday data analyses
  "styler", # source code formatter
  "formatR", # format output
  "arrow", # cross-language development platform to export .parquet
  "usethis", # automates repetitive tasks that arise during project setup
  "osfr", # interface for OSF
  "geomorph", # geometric morphometrics
  "jsonlite", # JSON parser and generator
  "httr"
)

purrr::walk(cran_pkgs, load_or_install)

# GitHub packages (supply the repo as github = "user/repo")
load_or_install("SlicerMorphR", github = "SlicerMorph/SlicerMorphR")
# import SlicerMorph dataset into R
```


``` r
#!setting decimals
fmt_decimals <- function(decimals = 0) {
  function(x) format(x, nsmall = decimals, scientific = FALSE)
}
# graphical settings for ggplot
my_theme <- theme(
  axis.text = element_text(size = 8, colour = "black"),
  # makes numbers smaller and black (consider final display)
  axis.ticks = element_line(
    linewidth = 0.5,
    colour = "black"
  ),
  # same for ticks
  axis.title = element_text(size = 10),
  # and for axis titles
  panel.grid.minor = element_blank(),
  panel.background = element_blank(),
  panel.grid.major = element_blank(),
  panel.border = element_rect(
    colour = "black", fill = NA, size = 0.5
  )
)

# calculate outliers
is_outlier <- function(x) {
  return(
    x < quantile(x, 0.25, na.rm = TRUE) - 1.5 * IQR(x, na.rm = TRUE) |
      x > quantile(x, 0.75, na.rm = TRUE) + 1.5 * IQR(x, na.rm = TRUE)
  )
}
```

Fragmented specimens were virtually pieced together in 3D slicer using the Fiducial Registration Wizard [@godinho2020]. Then, in the Markups Module of 3DSlicer, coordinates were extracted from a total of 21 anatomical landmarks from the most complete hemi-mandible of each specimen to capture mandibular morphology [@godinho2022]. The use of left hemimandibles was favoured, but because that was also the favoured side when sampling, theere were less left mandibles by the end of the process. Landmark coordinates of left mandibles were reflected to look like right (see Obsidian 3D Data)

# Preparing Landmark Data for `geomorph`

We implemented two parallel workflows for preparing landmark data into the format required by the `geomorph` package. The goal in both cases is to produce a 3D numeric array with dimensions *(p landmarks × k coordinates × n specimens)*, suitable for downstream analyses.

1.  **From `.mark.json` files**\
    Individual specimen landmark files produced by PAT were parsed from JSON into matrices of coordinates. These were then stacked into a 3D array, with each slice representing one specimen.

2.  **From `.csv` files**\
    Landmark coordinates saved as CSVs were read into matrices and similarly stacked into a 3D array with the same dimensional structure.

Both approaches produce the object `array3d` for use in `geomorph`, and use depends on whether user has access to raw (`.mark.json`) or derived (`.csv`) data on [OSF](https://osf.io/vkat9/) [@foster2017][^1]. Some safety checks were also included to ensure that the produced dataset adheres to `geomorph` criteria. THIS FAILED BECAUSE IT ONLY DOWNLAODS 84 files but there are 87.

[^1]: in the package `osfr`, PATs are required to upload files, create projects/components, access information about your private projects, or download files in your private projects. PATs are not required for accessing information about public projects or downloading public files, but authentication with a PAT will increase the rate limit on the API


``` r
# --- Dual-mode OSF pipeline for landmarks --- Project
# vkat9 | Raw Data: kwafd | Derived Data: 9fnsp | Analyses:
# a75wu

# --- Directories ---
raw_dir <- "data/raw"
derived_dir <- "data/derived"
usethis::use_directory(raw_dir)
usethis::use_directory(derived_dir)

derived_node <- osfr::osf_retrieve_node("9fnsp")
raw_node <- osfr::osf_retrieve_node("kwafd")

# 5) Gate OSF uploads so they only run on my machine
allow_flag <- tolower(Sys.getenv("ALLOW_OSF_UPLOAD"))
allow_upload <- nzchar(Sys.getenv("OSF_PAT")) && allow_flag %in%
    c("1", "true", "yes")

if (allow_upload) {
    cat("Authenticated with OSF PAT → building array from JSON\n")
    osfr::osf_auth(Sys.getenv("OSF_PAT"))

    # 1) Retrieve JSON files from Raw Data component
    json_files <- raw_node |>
        osfr::osf_ls_files(n = Inf) |>
        dplyr::filter(grepl("\\.json$", name, ignore.case = TRUE)) |>
        dplyr::distinct(id, .keep_all = TRUE)

    # 2) Download JSON files to raw_dir
    osfr::osf_download(json_files, path = raw_dir, conflicts = "overwrite")

    # 3) Specimen IDs from filenames
    json_paths <- list.files(raw_dir, pattern = "\\.json$", full.names = TRUE)
    specimen_ids <- tools::file_path_sans_ext(basename(json_paths))

    # 4) Read JSON into numeric matrices
    read_lmk_matrix <- function(path) {
        m <- SlicerMorphR::read.markups.json(path)
        m <- as.matrix(m)
        m <- apply(m, 2, as.numeric)
        if (is.null(colnames(m)))
            colnames(m) <- c("X", "Y", "Z")
        m
    }
    landmark_list <- purrr::map(json_paths, read_lmk_matrix) |>
        purrr::set_names(specimen_ids)

    # 5) Build array3d
    p <- nrow(landmark_list[[1]])
    k <- 3
    n <- length(landmark_list)

    array3d <- array(NA_real_, dim = c(p, k, n), dimnames = list(landmark = seq_len(p),
        coord = c("X", "Y", "Z"), specimen = names(landmark_list)))
    for (i in seq_along(landmark_list)) array3d[, , i] <- landmark_list[[i]]
    storage.mode(array3d) <- "double"

    # 6) Save array3d locally as RDS (only .RDS goes to
    # derived_dir)
    saveRDS(array3d, file.path(derived_dir, "array3d.RDS"))

    # 7) Create CSV + Parquet from array3d in-memory, write
    # to temp, upload to OSF, then remove temp files
    df2d <- data.frame(specimen = dimnames(array3d)$specimen,
        stringsAsFactors = FALSE)
    for (i in seq_len(p)) {
        df2d[[paste0("X", i)]] <- array3d[i, 1, ]
        df2d[[paste0("Y", i)]] <- array3d[i, 2, ]
        df2d[[paste0("Z", i)]] <- array3d[i, 3, ]
    }

    tmp_csv <- file.path(tempdir(), "landmarks.csv")
    tmp_parquet <- file.path(tempdir(), "landmarks.parquet")

    readr::write_csv(df2d, tmp_csv)
    arrow::write_parquet(df2d, tmp_parquet)

    osfr::osf_upload(derived_node, path = tmp_csv, conflicts = "overwrite")
    osfr::osf_upload(derived_node, path = tmp_parquet, conflicts = "overwrite")

    unlink(c(tmp_csv, tmp_parquet), force = TRUE)
} else {
    # --- No PAT branch ---
    cat("No OSF PAT → rebuilding array from Derived Data CSV on OSF\n")

    # 1) Locate CSV on OSF (do not save locally in project)
    csv_file <- derived_node |>
        osfr::osf_ls_files() |>
        dplyr::filter(name == "landmarks.csv")

    # 2) Download CSV to tempdir and read from there
    dl <- osfr::osf_download(csv_file, path = tempdir(), conflicts = "overwrite")
    csv_path <- dl$local_path
    df2d <- readr::read_csv(csv_path, show_col_types = FALSE)

    # 3) Rebuild array3d from CSV
    df_matrix <- as.matrix(dplyr::select(df2d, -specimen))

    p <- ncol(df_matrix)/3
    k <- 3
    n <- nrow(df_matrix)

    array3d <- geomorph::arrayspecs(df_matrix, p = p, k = k)

    # 4) Set dimnames identical to PAT branch
    dimnames(array3d)[[1]] <- seq_len(p)  # landmarks
    dimnames(array3d)[[2]] <- c("X", "Y", "Z")  # coords
    dimnames(array3d)[[3]] <- df2d$specimen  # specimen names

    # 5) Save array3d locally as RDS (only .RDS goes to
    # derived_dir)
    saveRDS(array3d, file.path(derived_dir, "array3d.RDS"))
}

# --- Safety checks ---
array_report <- function(array3d) {
    # Basic facts
    is_num <- is.numeric(array3d)
    dims <- dim(array3d)
    ndims <- if (is.null(dims))
        0L else length(dims)
    p <- if (ndims >= 1)
        dims[1] else NA_integer_
    k <- if (ndims >= 2)
        dims[2] else NA_integer_
    n <- if (ndims >= 3)
        dims[3] else NA_integer_

    has_3_dims <- ndims == 3
    coords_ok <- isTRUE(k %in% c(2, 3))
    spec_gt_coord <- isTRUE(n > k)

    # Specimen names in 3rd dim
    dn <- dimnames(array3d)
    spec_names <- if (!is.null(dn) && length(dn) >= 3)
        dn[[3]] else NULL
    has_spec_names <- !is.null(spec_names) && length(spec_names) ==
        n && all(!is.na(spec_names)) && all(nzchar(spec_names))

    # Non-finite counts (don’t fail if present; just
    # report)
    if (is_num) {
        na_count <- sum(is.na(array3d))
        nan_count <- sum(is.nan(array3d))
        inf_count <- sum(is.infinite(array3d))
    } else {
        # is.nan/is.infinite are numeric-only; report NA if
        # not numeric
        na_count <- sum(is.na(array3d))
        nan_count <- NA_integer_
        inf_count <- NA_integer_
    }

    # Print a simple, compact report
    cat("Array checks\n")
    cat("  Numeric: ", is_num, "\n", sep = "")
    cat("  3 dims: ", has_3_dims, "  dims: ", paste(ifelse(is.na(c(p,
        k, n)), "NA", c(p, k, n)), collapse = " x "), "\n", sep = "")
    cat("  Coords k in {2,3}: ", coords_ok, "\n", sep = "")
    cat("  Specimen names present: ", has_spec_names, "\n", sep = "")
    cat("  n specimens > n coords: ", spec_gt_coord, "\n", sep = "")
    cat("  Non-finite values: NA=", na_count, ", NaN=", nan_count,
        ", Inf=", inf_count, "\n", sep = "")

    invisible(list(is_numeric = is_num, has_3_dims = has_3_dims,
        dims = c(p = p, k = k, n = n), coords_ok = coords_ok,
        has_specimen_names = has_spec_names, specimens_gt_coords = spec_gt_coord,
        nonfinite = list(na = na_count, nan = nan_count, inf = inf_count)))
}

array_report(array3d)
```

<div data-pagedtable="false">
  <script data-pagedtable-source type="application/json">
{"columns":[{"label":["name"],"name":[1],"type":["chr"],"align":["left"]},{"label":["id"],"name":[2],"type":["chr"],"align":["left"]},{"label":["meta"],"name":[3],"type":["list"],"align":["right"]}],"data":[{"1":"landmarks.csv","2":"68b4a0e5fb038679a95e5f5a","3":"<named list [3]>"}],"options":{"columns":{"min":{},"max":[10]},"rows":{"min":[10],"max":[10]},"pages":{}}}
  </script>
</div><div data-pagedtable="false">
  <script data-pagedtable-source type="application/json">
{"columns":[{"label":["name"],"name":[1],"type":["chr"],"align":["left"]},{"label":["id"],"name":[2],"type":["chr"],"align":["left"]},{"label":["meta"],"name":[3],"type":["list"],"align":["right"]}],"data":[{"1":"landmarks.parquet","2":"68b4a0ea3560139ce496af8a","3":"<named list [3]>"}],"options":{"columns":{"min":{},"max":[10]},"rows":{"min":[10],"max":[10]},"pages":{}}}
  </script>
</div>


```
## Array checks
##   Numeric: TRUE
##   3 dims: TRUE  dims: 21 x 3 x 87
##   Coords k in {2,3}: TRUE
##   Specimen names present: TRUE
##   n specimens > n coords: TRUE
##   Non-finite values: NA=414, NaN=0, Inf=0
```


``` r
# Build specimen_metadata for exactly the specimens in
# array3d (no specimen-name changes), write
# data/derived/specimen_metadata.csv, list specimens with
# no sex AND no period AND no culture, and print per-period
# counts: total, number with missing == 0, and male/female
# counts.  Reads raw CSVs once. No library() calls; uses
# package namespaces.

# paths
scan_path <- "data/raw/scan.csv"
period_path <- "data/raw/period.csv"
demo_path <- "data/raw/demographics.csv"
out_meta_path <- "data/derived/specimen_metadata.csv"

# read raw inputs exactly once
scan_raw <- readr::read_csv(scan_path, show_col_types = FALSE)
period_raw <- readr::read_csv(period_path, show_col_types = FALSE)
demo_raw <- readr::read_csv(demo_path, show_col_types = FALSE)

# prepare tables (use specimen values exactly as provided;
# do not rename specimens)
scan <- scan_raw %>%
    dplyr::mutate(specimen = as.character(specimen), side = as.character(side),
        missing = as.numeric(missing)) %>%
    dplyr::select(specimen, side, missing)

period_tbl <- period_raw %>%
    dplyr::mutate(specimen = as.character(specimen)) %>%
    dplyr::filter(!is.na(specimen)) %>%
    dplyr::select(specimen, culture, period) %>%
    dplyr::distinct(specimen, .keep_all = TRUE) %>%
    dplyr::mutate(culture = as.character(culture), period = as.character(period))

demo_tbl <- demo_raw %>%
    dplyr::mutate(specimen = as.character(specimen)) %>%
    dplyr::filter(!is.na(specimen)) %>%
    dplyr::select(specimen, sex) %>%
    dplyr::distinct(specimen, .keep_all = TRUE) %>%
    dplyr::mutate(sex = as.character(sex))

# combined metadata (scan is master; exact matching by
# specimen)
combined <- scan %>%
    dplyr::left_join(period_tbl, by = "specimen") %>%
    dplyr::left_join(demo_tbl, by = "specimen") %>%
    dplyr::select(specimen, side, missing, culture, period, sex)

specimens_in_array <- as.character(dimnames(arr)[[3]])

# scaffold preserves array order and ensures exactly those
# specimens appear
scaffold <- data.frame(specimen = specimens_in_array, stringsAsFactors = FALSE)
meta_filtered <- scaffold %>%
    dplyr::left_join(combined, by = "specimen")

# write final combined metadata for the specimens in
# array3d (exact specimen names preserved)
readr::write_csv(meta_filtered, out_meta_path)
message("Wrote specimen metadata for specimens present in array3d -> ",
    out_meta_path)

# 1) specimens with no sex AND no period AND no culture
# (report their exact specimen IDs)
no_meta_mask <- is.na(meta_filtered$sex) & is.na(meta_filtered$period) &
    is.na(meta_filtered$culture)
no_meta_ids <- meta_filtered$specimen[which(no_meta_mask)]
if (length(no_meta_ids) == 0) {
    message("No specimens have all of sex, period and culture missing.")
} else {
    message("Specimens with no sex, no period, and no culture (exact specimen names):")
    cat(paste0(no_meta_ids, collapse = "\n"), "\n")
}

# 2) per-period summary: - total specimens per period
# (period kept as-is; NA shown as '<missing>') - number
# with missing == 0 (complete landmarks) - male and female
# counts per period (simple normalization by leading
# letter; does not rename specimens)
meta2 <- meta_filtered %>%
    dplyr::mutate(period_key = ifelse(is.na(period) | period ==
        "", "<missing>", as.character(period)), complete0_missing = ifelse(is.na(missing),
        FALSE, missing == 0), sex_norm = dplyr::case_when(is.na(sex) ~
        NA_character_, grepl("^[Ff]", sex) ~ "female", grepl("^[Mm]",
        sex) ~ "male", TRUE ~ "other"))

period_summary <- meta2 %>%
    dplyr::group_by(period_key) %>%
    dplyr::summarise(total = dplyr::n(), complete0_missing = sum(complete0_missing,
        na.rm = TRUE), n_male = sum(sex_norm == "male", na.rm = TRUE),
        n_female = sum(sex_norm == "female", na.rm = TRUE)) %>%
    dplyr::arrange(dplyr::desc(total))

cat("\nPer-period summary (period, total, complete0_missing, n_male, n_female):\n")
print(period_summary)

invisible(list(specimen_metadata = meta_filtered, no_meta_ids = no_meta_ids,
    period_summary = period_summary))
```


```
## Rows scan vs combined: 87 87 
## 
## No duplicates in scan_raw 
## 
## No duplicates in period_raw 
## 
## No duplicates in demo_raw 
## 
## Unmatched period rows: 22 
## # A tibble: 6 × 19
##   specimen    box   burial culture period labcode `14C`   `±` status  d13C  d15N
##   <chr>       <chr> <chr>  <chr>   <chr>  <chr>   <dbl> <dbl> <chr>  <dbl> <dbl>
## 1 ILH_Zvejni… 34-8  37     Narva   LM/EN  OxA-40…  6901    27 UNPUB… -23.8  12.7
## 2 ILH_Zvejni… 34-17 56     Kunda   MM     OxA-40…  7161    27 UNPUB… -21.4   9.9
## 3 ILH_Zvejni… 34-23 64     Kunda/… LM     NO DATE    NA    NA <NA>    NA    NA  
## 4 ILH_Zvejni… 34-33 92     Narva   LM/EN  Hela-1…  6510    50 PUBLI… -22.9  12  
## 5 ILH_Zvejni… 34-50 121    Narva   LM/EN  Ua-198…  6145    80 PUBLI… -23.1  10.9
## 6 ILH_Zvejni… 34-51 122    Narva   LM/EN  OxA-59…  6395    75 PUBLI… -23.3  11.9
## # ℹ 8 more variables: `C:N` <dbl>, ...13 <lgl>, ...14 <lgl>, ...15 <lgl>,
## #   ...16 <lgl>, ...17 <lgl>, ...18 <lgl>, ...19 <chr>
## Unmatched demographics rows: 20 
## # A tibble: 6 × 2
##   specimen                sex   
##   <chr>                   <chr> 
## 1 ILH_Zvejnieki_34.8.37   Male  
## 2 ILH_Zvejnieki_34.17.56  Male  
## 3 ILH_Zvejnieki_34.23.64  Male  
## 4 ILH_Zvejnieki_34.33.92  Male  
## 5 ILH_Zvejnieki_34.50.121 Female
## 6 ILH_Zvejnieki_34.52.123 <NA>  
## 
## Coverage counts:
## # A tibble: 1 × 4
##   n_specimens have_culture have_period have_sex
##         <int>        <int>       <int>    <int>
## 1          87           74          77       75
## 
## Coverage percentages:
## # A tibble: 3 × 4
##   n_specimens metric       count percent
##         <int> <chr>        <int>   <dbl>
## 1          87 have_culture    74    85.1
## 2          87 have_period     77    88.5
## 3          87 have_sex        75    86.2
## 
## Specimens with zero added metadata: 1 
## # A tibble: 1 × 1
##   specimen                
##   <chr>                   
## 1 ILH_Zvejnieki_34.176.105
## 
## NA counts per column:
## specimen     side  missing  culture   period      sex 
##        0        0        0       13       10       12 
## 
## Factor levels (NULL means column not factor yet):
## $side
## [1] "Left"  "Right"
## 
## $culture
## [1] "Comb Ware"   "Corded Ware" "Kunda"       "Kunda/Narva" "Narva"      
## 
## $period
## [1] "EIA"    "EN/MN"  "LBA/IA" "LM"     "LM/EN"  "LN"     "MM"     "MN"    
## 
## $sex
## [1] "Female" "Male"  
## 
## 
## Distribution of 'missing' values:
## # A tibble: 6 × 2
##   missing     n
##     <dbl> <int>
## 1       0    27
## 2       1    17
## 3       2    21
## 4       3     9
## 5       4    10
## 6       5     3
```

# Estimating `NA` using `estimate.missing()`

When specimens were incomplete, the location of the missing LMs was estimated using the thin plate spline (TPS) option of the function `estimate.missing()` of the `geomorph` package. Up to a maximum of 5 landmarks per specimen were estimated, since reconstruction of mandibles with more than 5 missing landmarks has been shown to result in reconstruction error and bias [@godinho2020a]. **TODO** Population specific reference mean specimens were used to reconstruct incomplete specimens as using inadequate references may also lead to large estimation errors [@neeser2009].


``` r
# estimate.missing() WARNING: n > k
estimateLMs <- estimate.missing(array3d, method = "TPS")
```

# Conduct GPA using `geomorph` and export using `SlicerMorphR`

Using the `gpapgen` function of the `geomorh` package, Generalized Procrustes Analysis (GPA) was used to superimpose all landmark configurations and remove the effects of location, size, and orientation on the raw coordinates. Using the `gm.prcomp` function of the `geomorph` package, Principal Component Analysis (PCA) was used to reduce dimensionality and examine shape differences between specimens.

The resulting shape variables were then used to examine morphological variance and hypothetical similarities and/or differences between groups. Using the `geomorph2slicermorph2` function of the `SlicerMorphR` package, the variables were exported so they could be visualised in 3D slicer by warping a surface along the relevant PCs. There's a bug in this, so that I have to go into `analysis.json` and remove `""`, such that change `“ExcludedLM”: "[]"` is `“ExcludedLM”: []` and `“SemiLandmarks”: "[]"` is `“SemiLandmarks”: []` , before uploading the data into 3DSlicer. Question on 3DSlicer Community on [this](https://discourse.slicer.org/t/slicer-morph-gpa-interactive-3d-vizualization-cannot-warp-around-pcs/44416/2).


``` r
# 1) GPA with fixed landmarks (ProcD = TRUE to
# compute/return Procrustes distances)
gpa <- geomorph::gpagen(A = estimateLMs, ProcD = TRUE)
```

```
## 
## Performing GPA
##   |                                                                              |                                                                      |   0%  |                                                                              |==================                                                    |  25%  |                                                                              |===================================                                   |  50%  |                                                                              |====================================================                  |  75%  |                                                                              |======================================================================| 100%
## 
## Making projections... Finished!
```

``` r
# 2) PCA on aligned coordinates
pca <- geomorph::gm.prcomp(gpa$coords)

# 3) Export SlicerMorph bundle to derived_dir Note:
# argument renamed to 'output_folder' per your current
# SlicerMorphR version
SlicerMorphR::geomorph2slicermorph2(gpa = gpa, pca = pca, output.folder = derived_dir)
```


``` r
# 4) upload to OSF if my OSF_PAT is present
if (allow_upload) {
    cat("OSF PAT branch: uploading to OSF node 9fnsp\n")

    # Upload all non-RDS artifacts from derived_dir
    # (recursively), overwrite on conflict
    files_to_upload <- list.files(derived_dir, recursive = TRUE,
        full.names = TRUE, include.dirs = FALSE)
    files_to_upload <- files_to_upload[!grepl("\\.rds$", files_to_upload,
        ignore.case = TRUE)]

    for (f in files_to_upload) {
        osfr::osf_upload(derived_node, path = f, conflicts = "overwrite")
    }
} else {
    cat("No OSF PAT branch: saved outputs locally in: ", derived_dir,
        "\n", sep = "")
}
```

# Conduct statistical analyses

Specimens were grouped based on chronology and region (i.e., Mesolithic Iberia, Levantine Chalcolithic, etc.). Potential differences in size (i.e., centroid size) were examined using the Kruskal–Wallis test. Shape differences were first tested using a nonparametric test Permutational Multivariate ANOVA (PERMANOVA) to assess potential multivariate shape differences in the different groups120. This was based on the first 29 PCs, which explain \~ 95% of the total variance and was implemented in Past121 using 10,000 permutations. The Kruskal–Wallis test (followed by post-hoc tests) was used to test for differences in the two first PCs, which were used together with surface warping to visualize shape differences across the groups. Kruskal–Wallis tests were implemented using the R package ggstatsplot122. The use of non-parametric testing was necessary for centroid size and PCs 1 and 2 after the Shapiro–Wilk’s test revealed that the ANOVA assumption of normality of residuals was not met for size, and Levene’s test of homogeneity of variances was not met for shape. These tests, along with the Durbin Watson test for independence of residuals, were carried out in the car R package123. Lastly, the morphol.disparity function of the geomorph package111 was used to examine hypothetical differences in shape variance across groups (i.e., disparity using Procrustes distances).

TODO: check out on paper stuff to produce graphs on ggplot using `make_ggplot`;


``` r
Y.gpa <- gpagen(plethodon$land, PrinAxes = FALSE)
gdf <- geomorph.data.frame(Y.gpa)
attributes(gdf) 
gdf <- geomorph.data.frame(Y.gpa, 
                           species = plethodon$species, 
                           site = plethodon$site) attributes(gdf) 
# Using geomorph.data.frame to facilitate analysis
anova(procD.lm(coords ~ Csize + species * site, data = gdf))
```

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
## time zone: Europe/Berlin
## tzcode source: internal
## 
## attached base packages:
## [1] stats     graphics  grDevices utils     datasets  methods   base     
## 
## other attached packages:
##  [1] SlicerMorphR_0.0.2.0000 httr_1.4.7              jsonlite_2.0.0         
##  [4] geomorph_4.0.10.999     Matrix_1.7-3            rgl_1.3.24             
##  [7] RRPP_2.1.2              osfr_0.2.9              usethis_3.2.0          
## [10] arrow_21.0.0.1          formatR_1.14            styler_1.10.3          
## [13] lubridate_1.9.4         forcats_1.0.0           stringr_1.5.1          
## [16] dplyr_1.1.4             purrr_1.1.0             readr_2.1.5            
## [19] tidyr_1.3.1             tibble_3.3.0            ggplot2_4.0.0          
## [22] tidyverse_2.0.0        
## 
## loaded via a namespace (and not attached):
##  [1] gtable_0.3.6       xfun_0.53          bslib_0.9.0        htmlwidgets_1.6.4 
##  [5] lattice_0.22-7     tzdb_0.5.0         vctrs_0.6.5        tools_4.5.1       
##  [9] generics_0.1.4     parallel_4.5.1     curl_7.0.0         pkgconfig_2.0.3   
## [13] R.oo_1.27.1        RColorBrewer_1.1-3 S7_0.2.0           assertthat_0.2.1  
## [17] lifecycle_1.0.4    R.cache_0.17.0     compiler_4.5.1     farver_2.1.2      
## [21] htmltools_0.5.8.1  sass_0.4.10        yaml_2.3.10        crayon_1.5.3      
## [25] pillar_1.11.0      jquerylib_0.1.4    R.utils_2.13.0     cachem_1.1.0      
## [29] mime_0.13          nlme_3.1-168       tidyselect_1.2.1   digest_0.6.37     
## [33] stringi_1.8.7      rprojroot_2.1.1    fastmap_1.2.0      grid_4.5.1        
## [37] cli_3.6.5          magrittr_2.0.3     base64enc_0.1-3    utf8_1.2.6        
## [41] triebeard_0.4.1    crul_1.6.0         ape_5.8-1          withr_3.0.2       
## [45] scales_1.4.0       bit64_4.6.0-1      timechange_0.3.0   rmarkdown_2.29    
## [49] jpeg_0.1-11        bit_4.6.0          R.methodsS3_1.8.2  hms_1.1.3         
## [53] memoise_2.0.1      evaluate_1.0.5     knitr_1.50         urltools_1.7.3.1  
## [57] rlang_1.1.6        Rcpp_1.1.0         glue_1.8.0         httpcode_0.3.0    
## [61] vroom_1.6.5        rstudioapi_0.17.1  R6_2.6.1           fs_1.6.6
```
