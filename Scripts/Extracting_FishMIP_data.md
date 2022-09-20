Extracting FishMIP Data
================
Denisse Fierro Arcos
2022-09-19

-   [Introduction](#introduction)
-   [Loading R libraries](#loading-r-libraries)
-   [Defining relevant directories](#defining-relevant-directories)
    -   [Obtaining list of files needed to create
        summaries](#obtaining-list-of-files-needed-to-create-summaries)
-   [Loading LME mask to extract
    data](#loading-lme-mask-to-extract-data)
-   [Defining function to merge files and extract
    data](#defining-function-to-merge-files-and-extract-data)

## Introduction

This notebook will guide you through the steps of how to extract data
using the `csv` masks created using the
[`Creating_Your_Own_Mask_From_Shapefiles`](%22Creating_Your_Own_Mask_From_Shapefiles.md%22).

## Loading R libraries

``` r
library(tidyverse)
# library(dtplyr)
# library(data.table)
library(here)
library(parallel)
#library(tictoc)
```

## Defining relevant directories

``` r
yannick_dir <- "/rd/gem/private/users/yannickr"
original_effort_dir <- "/rd/gem/private/users/yannickr/effort_mapped_bycountry"
aggregated_files_dir <- "/rd/gem/private/users/ldfierro/effort_mapped_by_country_aggregated/"
```

### Obtaining list of files needed to create summaries

``` r
#Files located under the original effort directory
original_effort_files <- list.files(original_effort_dir, pattern = ".csv")

#Checking we got the correct files
head(original_effort_files)
```

    ## [1] "mapped_1950_APW_100.csv" "mapped_1950_APW_104.csv"
    ## [3] "mapped_1950_APW_12.csv"  "mapped_1950_APW_120.csv"
    ## [5] "mapped_1950_APW_124.csv" "mapped_1950_APW_132.csv"

## Loading LME mask to extract data

For this example, we will be using the 0.5 degree mask for LMEs, which
matches the resolution of the outputs that we will be used for these
data extractions.

``` r
mask_df <- read_csv("../Data/Masks/fishMIP_regional_05deg_ISIMIP3a.csv")
```

    ## Rows: 16452 Columns: 3── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: ","
    ## chr (1): region
    ## dbl (2): Lon, Lat
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

``` r
#Check contents of mask
head(mask_df)
```

    ## # A tibble: 6 × 3
    ##     Lon   Lat region    
    ##   <dbl> <dbl> <chr>     
    ## 1  21.2  64.2 Baltic.Sea
    ## 2  21.8  64.2 Baltic.Sea
    ## 3  22.2  64.2 Baltic.Sea
    ## 4  22.8  64.2 Baltic.Sea
    ## 5  23.2  64.2 Baltic.Sea
    ## 6  20.8  63.8 Baltic.Sea

## Defining function to merge files and extract data

Given that we have 29852 files to process, we will create a function
before looping over all files.

<!-- #Set select from dplyr as default -->
<!-- select <- dplyr::select -->
<!-- # Defining function to merge files ---------------------------------------- -->
<!-- join_effort_data <- function(this_file_name, df, prefix_name){ -->
<!--   this_source_path <- file.path(original_effort_dir, this_file_name) -->
<!--   this_destination_path <- paste0(aggregated_files_dir, prefix_name, "_aggregated_", this_file_name) -->
<!--   if(!file.exists(this_destination_path)){ -->
<!--     Year <- as.numeric(str_extract(this_file_name, pattern =  "([[:digit:]])+")) -->
<!--     these_data <- read_csv(this_source_path) -->
<!--     #dtplyr approach -->
<!--     these_aggregated_data <-  these_data %>% -->
<!--       left_join(df, by=c("Lat", "Lon")) %>% -->
<!--       group_by(region, SAUP, Gear, FGroup, Sector) %>% -->
<!--       summarise(NomActive = sum(NomActive, na.rm = TRUE), -->
<!--                 EffActive = sum(EffActive, na.rm = TRUE), -->
<!--                 NV= sum(NV, na.rm = TRUE), -->
<!--                 P= sum(P, na.rm = TRUE), -->
<!--                 GT= sum(GT, na.rm = TRUE)) %>% -->
<!--       mutate(Year = Year) %>%  -->
<!--       as.data.table() %>%  -->
<!--       lazy_dt() -->
<!--     fwrite(x = these_aggregated_data, file = this_destination_path) -->
<!--   } -->
<!-- } -->
<!-- # Parallelising work ---------------------------------------------------- -->
<!-- chunk_size <- 500 #chunk size for processing -->
<!-- effort_list_split <- split(original_effort_files, ceiling(seq_along(original_effort_files)/chunk_size)) -->
<!-- length(effort_list_split) #60 chunks -->
<!-- tic() -->
<!-- for(i in 1:length(effort_list_split)){ -->
<!--   file_chunk <- effort_list_split[[i]] -->
<!--   message("Processing chunk #", i, " of ", length(effort_list_split)) -->
<!--   mclapply(X = file_chunk, FUN = join_effort_data, deg1_df, "1deg", mc.cores = 40) -->
<!--   mclapply(X = file_chunk, FUN = join_effort_data, deg025_df, "025deg", mc.cores = 40) -->
<!-- } -->
<!-- toc() -->
<!-- #check the files -->
<!-- newly_written_files_1deg <- list.files(aggregated_files_dir, pattern = "1deg", full.names = TRUE) -->
<!-- newly_written_files_025deg <- list.files(aggregated_files_dir, pattern = "025deg", full.names = TRUE) -->
<!-- #pick one randomly  -->
<!-- map(newly_written_files[[500]], fread) -->
<!-- combined_aggregated_effort_1deg <- rbindlist(mclapply(X = newly_written_files_1deg, FUN = fread, mc.cores = 40)) -->
<!-- fwrite(x = combined_aggregated_effort_1deg, file.path(aggregated_files_dir, "1deg_all_effort_aggregated.csv")) -->
<!-- combined_aggregated_effort_025deg <- rbindlist(mclapply(X = newly_written_files_025deg, FUN = fread, mc.cores = 40)) -->
<!-- fwrite(x = combined_aggregated_effort_025deg, file.path(aggregated_files_dir, "025deg_all_effort_aggregated.csv")) -->
