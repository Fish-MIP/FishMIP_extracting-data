Creating raster mask from FAO regions shapefile
================
Denisse Fierro Arcos
2022-11-07

- <a href="#introduction" id="toc-introduction">Introduction</a>
- <a href="#loading-r-libraries" id="toc-loading-r-libraries">Loading R
  libraries</a>
- <a href="#loading-shapefiles" id="toc-loading-shapefiles">Loading
  shapefiles</a>
- <a href="#plotting-fao-regions-shapefile"
  id="toc-plotting-fao-regions-shapefile">Plotting FAO regions
  shapefile</a>
- <a href="#extracting-names-and-codes-for-fao-regions"
  id="toc-extracting-names-and-codes-for-fao-regions">Extracting names and
  codes for FAO regions</a>
- <a href="#creating-a-multilayer-raster-mask-based-on-merged-shapefile"
  id="toc-creating-a-multilayer-raster-mask-based-on-merged-shapefile">Creating
  a multilayer raster mask based on merged shapefile</a>
  - <a href="#loading-input-rasters" id="toc-loading-input-rasters">Loading
    input rasters</a>
  - <a href="#defining-function-to-create-rasters-from-shapefiles"
    id="toc-defining-function-to-create-rasters-from-shapefiles">Defining
    function to create rasters from shapefiles</a>
  - <a href="#applying-function-to-list-containing-all-shapefiles"
    id="toc-applying-function-to-list-containing-all-shapefiles">Applying
    function to list containing all shapefiles</a>
- <a href="#python-based-code"
  id="toc-python-based-code"><code>Python</code>-based code</a>
  - <a href="#loading-libraries" id="toc-loading-libraries">Loading
    libraries</a>
  - <a href="#loading-raster-using-xarray"
    id="toc-loading-raster-using-xarray">Loading raster using
    <code>xarray</code></a>
  - <a href="#plotting-results" id="toc-plotting-results">Plotting
    results</a>

## Introduction

This notebook will go through the process of creating a raster mask of
the FAO major fishing areas available through its [WFS
service](https://www.fao.org/fishery/geoserver/fifao/ows?service=WFS&request=GetFeature&version=1.0.0&typeName=fifao:FAO_AREAS_CWP).
Rasters are created in `R` and the final details are done with `xarray`
library for `Python`.

## Loading R libraries

``` r
library(sf)
library(raster)
library(tidyverse)
library(reticulate)
```

## Loading shapefiles

A copy of the FAO Major Fishing Areas was saved locally from the [FAO’s
WFS
service](https://www.fao.org/fishery/geoserver/fifao/ows?service=WFS&request=GetFeature&version=1.0.0&typeName=fifao:FAO_AREAS_CWP).

``` r
#Load shapefile with FAO regions
fao_reg <- read_sf("../Spatial_Data/FAO_shapefiles/FAO_MajorAreas.shp") %>% 
  #Subset of columns
  select(-c(F_LEVEL, F_STATUS, SUBOCEAN:F_SUBUNIT, NAME_FR:SURFACE)) %>% 
  #Turning character columns into factors
  mutate_if(is.character, as.factor)

#We can check the results of the first two rows
head(fao_reg, 2)
```

    ## Simple feature collection with 2 features and 3 fields
    ## Geometry type: POLYGON
    ## Dimension:     XY
    ## Bounding box:  xmin: -180 ymin: 35.7829 xmax: 126.5826 ymax: 89.99
    ## Geodetic CRS:  WGS 84
    ## # A tibble: 2 × 4
    ##   F_CODE OCEAN   NAME_EN                                                geometry
    ##   <fct>  <fct>   <fct>                                             <POLYGON [°]>
    ## 1 18     Arctic  Arctic Sea         ((-40 89.99, -40 83.2997, -40.0028 83.2998,…
    ## 2 61     Pacific Pacific, Northwest ((126.5826 35.7829, 126.5825 35.7829, 126.5…

## Plotting FAO regions shapefile

``` r
#Loading land shapefile to include in plot
land <- rnaturalearth::ne_countries(type = "countries", returnclass = "sf")
#Plotting FAO shapefile
fao_reg %>% 
  ggplot()+
  geom_sf(aes(fill = NAME_EN))+
  geom_sf(data = land, inherit.aes = F, color = "gray")+
  theme_bw()
```

![](Creating_Mask_FAO_regions_files/figure-gfm/plot_shapefile-1.png)<!-- -->

## Extracting names and codes for FAO regions

We will use this information to save correct names for regions in the
raster

``` r
#Create a data frame of unique FAO regions
fao_names_codes <- fao_reg %>% 
  #Tranforming into data frame
  st_drop_geometry() %>% 
  #Extracting unique regions
  distinct(F_CODE, NAME_EN) %>% 
  #Ordering by region code
  arrange(F_CODE) %>% 
  #Fixing up names to use as raster layer name
  mutate_all(as.character) %>% 
  #Removing any spaces and commas from names and replacing with an underscore "_"
  mutate(NAME_EN = str_remove(string = NAME_EN, pattern = ","), 
         NAME_EN = str_replace_all(NAME_EN, " ", "_"))

#We can check some of the results
head(fao_names_codes, 2)
```

    ## # A tibble: 2 × 2
    ##   F_CODE NAME_EN           
    ##   <chr>  <chr>             
    ## 1 18     Arctic_Sea        
    ## 2 21     Atlantic_Northwest

``` r
#Saving FAO region keys
fao_names_codes %>% 
  write_csv("../Spatial_Data/FAO_shapefiles/FAO-regions_keys.csv")
```

## Creating a multilayer raster mask based on merged shapefile

We will now create a multilayer mask that matches the grid used in the
physical model forcings. Most models use the same 1 degree grid, with
the exception of the DBPM ecosystem model, which uses a different 1
degree grid, and the DBEM model, which uses a 0.5 degree grid.

Note that masks *must* match the grid of the model from which data is
being extracted. This means that you will need to create a new mask for
each grid that is different. In the chunk below, you will find the three
different grids identified in the ecosystem models.

### Loading input rasters

``` r
#Loading sample raster to be used as target for rasterising FAO regions
#Most ecosystem models use this one degree grid
# ras <- raster("../Spatial_Data/InputRasters/gfdl-mom6-cobalt2_obsclim_deptho_onedeg_global_fixed.nc")

#Sample from DBPM model
ras <- raster("../Spatial_Data/InputRasters/dbpm_ipsl-cm6a-lr_nobasd_historical_nat_default_tcb_global_monthly_1850_2014.nc")[[1]]
```

    ## Loading required namespace: ncdf4

``` r
#Sample from DBEM model
# ras <- raster("../Spatial_Data/InputRasters/dbem_ipsl-cm6a-lr_nobasd_historical_nat_default_tcb_global_annual_1951_2014.nc")[[1]]

#We will define a few extra variables to automate creation of file names for each mask
#Model resolution
res <- "1deg"

#Model associated with grid. Leave blank if multiple models use same grid
mod_name <- "_DBPM"

#Plotting raster
plot(ras)
```

![](Creating_Mask_FAO_regions_files/figure-gfm/load_rasters_input-1.png)<!-- -->

### Defining function to create rasters from shapefiles

We will define our own function that will use the shapefiles above to
create rasters.

``` r
#Defining function which needs a shapefile and a raster as input
shp_to_raster <- function(shp, nc_raster){
  #The final raster will have ones where within the shapefile boundaries
  rasterize(shp, nc_raster, field = 1)
}
```

### Applying function to list containing all shapefiles

``` r
#Split shapefile into regions prior to transforming into raster
fao_reg_list <- fao_reg %>% 
  group_by(F_CODE) %>% 
  group_split()

#Applying function to raster list
grid_raster <- map(fao_reg_list, shp_to_raster, ras) %>% 
  #Stacking rasters to create a single multilayer raster
  stack()

#Checking results of stacked raster (first six regions)...
plot(grid_raster[[1:6]])
```

![](Creating_Mask_FAO_regions_files/figure-gfm/mapping_function-1.png)<!-- -->
We will plot the shapefile to compare results.

``` r
#...against shapefile, to make sure they match
fao_reg %>% 
  #Creating a new column with the code and name to identify them easily in graph
  unite("reg_code", F_CODE, NAME_EN, remove = F) %>% 
  ggplot()+
  geom_sf(aes(fill = reg_code))+
  facet_wrap(~F_CODE)
```

![](Creating_Mask_FAO_regions_files/figure-gfm/unnamed-chunk-1-1.png)<!-- -->
The results match! We can now save the results locally and move onto
`Python` for the final touches.

``` r
#Define file name
filename <- paste0("FAO-regions_", res, "mask", mod_name, ".nc")

#Saving raster to disk
writeRaster(grid_raster, file.path("../Spatial_Data/Masks", filename), format = "CDF", overwrite = T,
            varname = "FAO_regions", zname = "FAO_reg_name")
```

# `Python`-based code

We will now start `Python` and save the correct names for the FAO
regions in the `netcdf` file we created in `R`.

``` r
#Activating conda
use_condaenv("CMIP6_data")
```

## Loading libraries

``` python
import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import os
```

## Loading raster using `xarray`

``` python
#Getting filename
fn = f'FAO-regions_{r.res}mask{r.mod_name}.nc'

#Loading multilayer raster as dataset
mask = xr.open_dataset(os.path.join('../Spatial_Data/Masks', fn))
#Checking saved file in R
mask
```

<div><svg style="position: absolute; width: 0; height: 0; overflow: hidden">
<defs>
<symbol id="icon-database" viewBox="0 0 32 32">
<path d="M16 0c-8.837 0-16 2.239-16 5v4c0 2.761 7.163 5 16 5s16-2.239 16-5v-4c0-2.761-7.163-5-16-5z"></path>
<path d="M16 17c-8.837 0-16-2.239-16-5v6c0 2.761 7.163 5 16 5s16-2.239 16-5v-6c0 2.761-7.163 5-16 5z"></path>
<path d="M16 26c-8.837 0-16-2.239-16-5v6c0 2.761 7.163 5 16 5s16-2.239 16-5v-6c0 2.761-7.163 5-16 5z"></path>
</symbol>
<symbol id="icon-file-text2" viewBox="0 0 32 32">
<path d="M28.681 7.159c-0.694-0.947-1.662-2.053-2.724-3.116s-2.169-2.030-3.116-2.724c-1.612-1.182-2.393-1.319-2.841-1.319h-15.5c-1.378 0-2.5 1.121-2.5 2.5v27c0 1.378 1.122 2.5 2.5 2.5h23c1.378 0 2.5-1.122 2.5-2.5v-19.5c0-0.448-0.137-1.23-1.319-2.841zM24.543 5.457c0.959 0.959 1.712 1.825 2.268 2.543h-4.811v-4.811c0.718 0.556 1.584 1.309 2.543 2.268zM28 29.5c0 0.271-0.229 0.5-0.5 0.5h-23c-0.271 0-0.5-0.229-0.5-0.5v-27c0-0.271 0.229-0.5 0.5-0.5 0 0 15.499-0 15.5 0v7c0 0.552 0.448 1 1 1h7v19.5z"></path>
<path d="M23 26h-14c-0.552 0-1-0.448-1-1s0.448-1 1-1h14c0.552 0 1 0.448 1 1s-0.448 1-1 1z"></path>
<path d="M23 22h-14c-0.552 0-1-0.448-1-1s0.448-1 1-1h14c0.552 0 1 0.448 1 1s-0.448 1-1 1z"></path>
<path d="M23 18h-14c-0.552 0-1-0.448-1-1s0.448-1 1-1h14c0.552 0 1 0.448 1 1s-0.448 1-1 1z"></path>
</symbol>
</defs>
</svg>
<style>/* CSS stylesheet for displaying xarray objects in jupyterlab.
 *
 */

:root {
  --xr-font-color0: var(--jp-content-font-color0, rgba(0, 0, 0, 1));
  --xr-font-color2: var(--jp-content-font-color2, rgba(0, 0, 0, 0.54));
  --xr-font-color3: var(--jp-content-font-color3, rgba(0, 0, 0, 0.38));
  --xr-border-color: var(--jp-border-color2, #e0e0e0);
  --xr-disabled-color: var(--jp-layout-color3, #bdbdbd);
  --xr-background-color: var(--jp-layout-color0, white);
  --xr-background-color-row-even: var(--jp-layout-color1, white);
  --xr-background-color-row-odd: var(--jp-layout-color2, #eeeeee);
}

html[theme=dark],
body[data-theme=dark],
body.vscode-dark {
  --xr-font-color0: rgba(255, 255, 255, 1);
  --xr-font-color2: rgba(255, 255, 255, 0.54);
  --xr-font-color3: rgba(255, 255, 255, 0.38);
  --xr-border-color: #1F1F1F;
  --xr-disabled-color: #515151;
  --xr-background-color: #111111;
  --xr-background-color-row-even: #111111;
  --xr-background-color-row-odd: #313131;
}

.xr-wrap {
  display: block !important;
  min-width: 300px;
  max-width: 700px;
}

.xr-text-repr-fallback {
  /* fallback to plain text repr when CSS is not injected (untrusted notebook) */
  display: none;
}

.xr-header {
  padding-top: 6px;
  padding-bottom: 6px;
  margin-bottom: 4px;
  border-bottom: solid 1px var(--xr-border-color);
}

.xr-header > div,
.xr-header > ul {
  display: inline;
  margin-top: 0;
  margin-bottom: 0;
}

.xr-obj-type,
.xr-array-name {
  margin-left: 2px;
  margin-right: 10px;
}

.xr-obj-type {
  color: var(--xr-font-color2);
}

.xr-sections {
  padding-left: 0 !important;
  display: grid;
  grid-template-columns: 150px auto auto 1fr 20px 20px;
}

.xr-section-item {
  display: contents;
}

.xr-section-item input {
  display: none;
}

.xr-section-item input + label {
  color: var(--xr-disabled-color);
}

.xr-section-item input:enabled + label {
  cursor: pointer;
  color: var(--xr-font-color2);
}

.xr-section-item input:enabled + label:hover {
  color: var(--xr-font-color0);
}

.xr-section-summary {
  grid-column: 1;
  color: var(--xr-font-color2);
  font-weight: 500;
}

.xr-section-summary > span {
  display: inline-block;
  padding-left: 0.5em;
}

.xr-section-summary-in:disabled + label {
  color: var(--xr-font-color2);
}

.xr-section-summary-in + label:before {
  display: inline-block;
  content: '►';
  font-size: 11px;
  width: 15px;
  text-align: center;
}

.xr-section-summary-in:disabled + label:before {
  color: var(--xr-disabled-color);
}

.xr-section-summary-in:checked + label:before {
  content: '▼';
}

.xr-section-summary-in:checked + label > span {
  display: none;
}

.xr-section-summary,
.xr-section-inline-details {
  padding-top: 4px;
  padding-bottom: 4px;
}

.xr-section-inline-details {
  grid-column: 2 / -1;
}

.xr-section-details {
  display: none;
  grid-column: 1 / -1;
  margin-bottom: 5px;
}

.xr-section-summary-in:checked ~ .xr-section-details {
  display: contents;
}

.xr-array-wrap {
  grid-column: 1 / -1;
  display: grid;
  grid-template-columns: 20px auto;
}

.xr-array-wrap > label {
  grid-column: 1;
  vertical-align: top;
}

.xr-preview {
  color: var(--xr-font-color3);
}

.xr-array-preview,
.xr-array-data {
  padding: 0 5px !important;
  grid-column: 2;
}

.xr-array-data,
.xr-array-in:checked ~ .xr-array-preview {
  display: none;
}

.xr-array-in:checked ~ .xr-array-data,
.xr-array-preview {
  display: inline-block;
}

.xr-dim-list {
  display: inline-block !important;
  list-style: none;
  padding: 0 !important;
  margin: 0;
}

.xr-dim-list li {
  display: inline-block;
  padding: 0;
  margin: 0;
}

.xr-dim-list:before {
  content: '(';
}

.xr-dim-list:after {
  content: ')';
}

.xr-dim-list li:not(:last-child):after {
  content: ',';
  padding-right: 5px;
}

.xr-has-index {
  font-weight: bold;
}

.xr-var-list,
.xr-var-item {
  display: contents;
}

.xr-var-item > div,
.xr-var-item label,
.xr-var-item > .xr-var-name span {
  background-color: var(--xr-background-color-row-even);
  margin-bottom: 0;
}

.xr-var-item > .xr-var-name:hover span {
  padding-right: 5px;
}

.xr-var-list > li:nth-child(odd) > div,
.xr-var-list > li:nth-child(odd) > label,
.xr-var-list > li:nth-child(odd) > .xr-var-name span {
  background-color: var(--xr-background-color-row-odd);
}

.xr-var-name {
  grid-column: 1;
}

.xr-var-dims {
  grid-column: 2;
}

.xr-var-dtype {
  grid-column: 3;
  text-align: right;
  color: var(--xr-font-color2);
}

.xr-var-preview {
  grid-column: 4;
}

.xr-var-name,
.xr-var-dims,
.xr-var-dtype,
.xr-preview,
.xr-attrs dt {
  white-space: nowrap;
  overflow: hidden;
  text-overflow: ellipsis;
  padding-right: 10px;
}

.xr-var-name:hover,
.xr-var-dims:hover,
.xr-var-dtype:hover,
.xr-attrs dt:hover {
  overflow: visible;
  width: auto;
  z-index: 1;
}

.xr-var-attrs,
.xr-var-data {
  display: none;
  background-color: var(--xr-background-color) !important;
  padding-bottom: 5px !important;
}

.xr-var-attrs-in:checked ~ .xr-var-attrs,
.xr-var-data-in:checked ~ .xr-var-data {
  display: block;
}

.xr-var-data > table {
  float: right;
}

.xr-var-name span,
.xr-var-data,
.xr-attrs {
  padding-left: 25px !important;
}

.xr-attrs,
.xr-var-attrs,
.xr-var-data {
  grid-column: 1 / -1;
}

dl.xr-attrs {
  padding: 0;
  margin: 0;
  display: grid;
  grid-template-columns: 125px auto;
}

.xr-attrs dt,
.xr-attrs dd {
  padding: 0;
  margin: 0;
  float: left;
  padding-right: 10px;
  width: auto;
}

.xr-attrs dt {
  font-weight: normal;
  grid-column: 1;
}

.xr-attrs dt:hover span {
  display: inline-block;
  background: var(--xr-background-color);
  padding-right: 10px;
}

.xr-attrs dd {
  grid-column: 2;
  white-space: pre-wrap;
  word-break: break-all;
}

.xr-icon-database,
.xr-icon-file-text2 {
  display: inline-block;
  vertical-align: middle;
  width: 1em;
  height: 1.5em !important;
  stroke-width: 0;
  stroke: currentColor;
  fill: currentColor;
}
</style><pre class='xr-text-repr-fallback'>&lt;xarray.Dataset&gt;
Dimensions:       (longitude: 360, latitude: 180, FAO_reg_name: 19)
Coordinates:
  * longitude     (longitude) float64 -180.0 -179.0 -178.0 ... 177.0 178.0 179.0
  * latitude      (latitude) float64 89.5 88.5 87.5 86.5 ... -87.5 -88.5 -89.5
  * FAO_reg_name  (FAO_reg_name) int32 1 2 3 4 5 6 7 8 ... 13 14 15 16 17 18 19
Data variables:
    crs           int32 -2147483647
    FAO_regions   (FAO_reg_name, latitude, longitude) float32 ...
Attributes:
    Conventions:  CF-1.4
    created_by:   R, packages ncdf4 and raster (version 3.5-29)
    date:         2022-12-15 19:36:13</pre><div class='xr-wrap' style='display:none'><div class='xr-header'><div class='xr-obj-type'>xarray.Dataset</div></div><ul class='xr-sections'><li class='xr-section-item'><input id='section-684540aa-01c2-4318-b3e3-e77d4ff3308e' class='xr-section-summary-in' type='checkbox' disabled ><label for='section-684540aa-01c2-4318-b3e3-e77d4ff3308e' class='xr-section-summary'  title='Expand/collapse section'>Dimensions:</label><div class='xr-section-inline-details'><ul class='xr-dim-list'><li><span class='xr-has-index'>longitude</span>: 360</li><li><span class='xr-has-index'>latitude</span>: 180</li><li><span class='xr-has-index'>FAO_reg_name</span>: 19</li></ul></div><div class='xr-section-details'></div></li><li class='xr-section-item'><input id='section-634b5874-a888-4b5b-ae05-fbc8c15d0a90' class='xr-section-summary-in' type='checkbox'  checked><label for='section-634b5874-a888-4b5b-ae05-fbc8c15d0a90' class='xr-section-summary' >Coordinates: <span>(3)</span></label><div class='xr-section-inline-details'></div><div class='xr-section-details'><ul class='xr-var-list'><li class='xr-var-item'><div class='xr-var-name'><span class='xr-has-index'>longitude</span></div><div class='xr-var-dims'>(longitude)</div><div class='xr-var-dtype'>float64</div><div class='xr-var-preview xr-preview'>-180.0 -179.0 ... 178.0 179.0</div><input id='attrs-f052a30f-684e-433b-b4d5-407ec1ca0860' class='xr-var-attrs-in' type='checkbox' ><label for='attrs-f052a30f-684e-433b-b4d5-407ec1ca0860' title='Show/Hide attributes'><svg class='icon xr-icon-file-text2'><use xlink:href='#icon-file-text2'></use></svg></label><input id='data-5250beaa-db80-4147-a762-c0c6230a2a00' class='xr-var-data-in' type='checkbox'><label for='data-5250beaa-db80-4147-a762-c0c6230a2a00' title='Show/Hide data repr'><svg class='icon xr-icon-database'><use xlink:href='#icon-database'></use></svg></label><div class='xr-var-attrs'><dl class='xr-attrs'><dt><span>units :</span></dt><dd>degrees_east</dd><dt><span>long_name :</span></dt><dd>longitude</dd></dl></div><div class='xr-var-data'><pre>array([-180., -179., -178., ...,  177.,  178.,  179.])</pre></div></li><li class='xr-var-item'><div class='xr-var-name'><span class='xr-has-index'>latitude</span></div><div class='xr-var-dims'>(latitude)</div><div class='xr-var-dtype'>float64</div><div class='xr-var-preview xr-preview'>89.5 88.5 87.5 ... -88.5 -89.5</div><input id='attrs-6f90d807-c5f1-4d0f-8319-a6f17fe09337' class='xr-var-attrs-in' type='checkbox' ><label for='attrs-6f90d807-c5f1-4d0f-8319-a6f17fe09337' title='Show/Hide attributes'><svg class='icon xr-icon-file-text2'><use xlink:href='#icon-file-text2'></use></svg></label><input id='data-e8e7bbde-d4e0-455b-894f-2a403b88af0e' class='xr-var-data-in' type='checkbox'><label for='data-e8e7bbde-d4e0-455b-894f-2a403b88af0e' title='Show/Hide data repr'><svg class='icon xr-icon-database'><use xlink:href='#icon-database'></use></svg></label><div class='xr-var-attrs'><dl class='xr-attrs'><dt><span>units :</span></dt><dd>degrees_north</dd><dt><span>long_name :</span></dt><dd>latitude</dd></dl></div><div class='xr-var-data'><pre>array([ 89.5,  88.5,  87.5,  86.5,  85.5,  84.5,  83.5,  82.5,  81.5,  80.5,
        79.5,  78.5,  77.5,  76.5,  75.5,  74.5,  73.5,  72.5,  71.5,  70.5,
        69.5,  68.5,  67.5,  66.5,  65.5,  64.5,  63.5,  62.5,  61.5,  60.5,
        59.5,  58.5,  57.5,  56.5,  55.5,  54.5,  53.5,  52.5,  51.5,  50.5,
        49.5,  48.5,  47.5,  46.5,  45.5,  44.5,  43.5,  42.5,  41.5,  40.5,
        39.5,  38.5,  37.5,  36.5,  35.5,  34.5,  33.5,  32.5,  31.5,  30.5,
        29.5,  28.5,  27.5,  26.5,  25.5,  24.5,  23.5,  22.5,  21.5,  20.5,
        19.5,  18.5,  17.5,  16.5,  15.5,  14.5,  13.5,  12.5,  11.5,  10.5,
         9.5,   8.5,   7.5,   6.5,   5.5,   4.5,   3.5,   2.5,   1.5,   0.5,
        -0.5,  -1.5,  -2.5,  -3.5,  -4.5,  -5.5,  -6.5,  -7.5,  -8.5,  -9.5,
       -10.5, -11.5, -12.5, -13.5, -14.5, -15.5, -16.5, -17.5, -18.5, -19.5,
       -20.5, -21.5, -22.5, -23.5, -24.5, -25.5, -26.5, -27.5, -28.5, -29.5,
       -30.5, -31.5, -32.5, -33.5, -34.5, -35.5, -36.5, -37.5, -38.5, -39.5,
       -40.5, -41.5, -42.5, -43.5, -44.5, -45.5, -46.5, -47.5, -48.5, -49.5,
       -50.5, -51.5, -52.5, -53.5, -54.5, -55.5, -56.5, -57.5, -58.5, -59.5,
       -60.5, -61.5, -62.5, -63.5, -64.5, -65.5, -66.5, -67.5, -68.5, -69.5,
       -70.5, -71.5, -72.5, -73.5, -74.5, -75.5, -76.5, -77.5, -78.5, -79.5,
       -80.5, -81.5, -82.5, -83.5, -84.5, -85.5, -86.5, -87.5, -88.5, -89.5])</pre></div></li><li class='xr-var-item'><div class='xr-var-name'><span class='xr-has-index'>FAO_reg_name</span></div><div class='xr-var-dims'>(FAO_reg_name)</div><div class='xr-var-dtype'>int32</div><div class='xr-var-preview xr-preview'>1 2 3 4 5 6 7 ... 14 15 16 17 18 19</div><input id='attrs-7a2fcc08-732c-491b-bb96-eef0872cac28' class='xr-var-attrs-in' type='checkbox' ><label for='attrs-7a2fcc08-732c-491b-bb96-eef0872cac28' title='Show/Hide attributes'><svg class='icon xr-icon-file-text2'><use xlink:href='#icon-file-text2'></use></svg></label><input id='data-f32e7593-66ef-434f-a74e-4de1fb116c71' class='xr-var-data-in' type='checkbox'><label for='data-f32e7593-66ef-434f-a74e-4de1fb116c71' title='Show/Hide data repr'><svg class='icon xr-icon-database'><use xlink:href='#icon-database'></use></svg></label><div class='xr-var-attrs'><dl class='xr-attrs'><dt><span>units :</span></dt><dd>unknown</dd><dt><span>long_name :</span></dt><dd>FAO_reg_name</dd></dl></div><div class='xr-var-data'><pre>array([ 1,  2,  3,  4,  5,  6,  7,  8,  9, 10, 11, 12, 13, 14, 15, 16, 17, 18,
       19], dtype=int32)</pre></div></li></ul></div></li><li class='xr-section-item'><input id='section-7941acf2-9d5e-4220-8add-eed9ecfc7ab3' class='xr-section-summary-in' type='checkbox'  checked><label for='section-7941acf2-9d5e-4220-8add-eed9ecfc7ab3' class='xr-section-summary' >Data variables: <span>(2)</span></label><div class='xr-section-inline-details'></div><div class='xr-section-details'><ul class='xr-var-list'><li class='xr-var-item'><div class='xr-var-name'><span>crs</span></div><div class='xr-var-dims'>()</div><div class='xr-var-dtype'>int32</div><div class='xr-var-preview xr-preview'>...</div><input id='attrs-2adc6267-9249-48fa-9cea-55c82a90575c' class='xr-var-attrs-in' type='checkbox' ><label for='attrs-2adc6267-9249-48fa-9cea-55c82a90575c' title='Show/Hide attributes'><svg class='icon xr-icon-file-text2'><use xlink:href='#icon-file-text2'></use></svg></label><input id='data-155451fd-9fe7-4f01-9d3d-7a331b48624b' class='xr-var-data-in' type='checkbox'><label for='data-155451fd-9fe7-4f01-9d3d-7a331b48624b' title='Show/Hide data repr'><svg class='icon xr-icon-database'><use xlink:href='#icon-database'></use></svg></label><div class='xr-var-attrs'><dl class='xr-attrs'><dt><span>proj4 :</span></dt><dd>+proj=longlat +datum=WGS84 +no_defs</dd></dl></div><div class='xr-var-data'><pre>array(-2147483647, dtype=int32)</pre></div></li><li class='xr-var-item'><div class='xr-var-name'><span>FAO_regions</span></div><div class='xr-var-dims'>(FAO_reg_name, latitude, longitude)</div><div class='xr-var-dtype'>float32</div><div class='xr-var-preview xr-preview'>...</div><input id='attrs-8971015b-a99f-4f91-a468-cff2cf7d87f2' class='xr-var-attrs-in' type='checkbox' ><label for='attrs-8971015b-a99f-4f91-a468-cff2cf7d87f2' title='Show/Hide attributes'><svg class='icon xr-icon-file-text2'><use xlink:href='#icon-file-text2'></use></svg></label><input id='data-3e3518d5-e6be-4e1d-9996-5ede90742adf' class='xr-var-data-in' type='checkbox'><label for='data-3e3518d5-e6be-4e1d-9996-5ede90742adf' title='Show/Hide data repr'><svg class='icon xr-icon-database'><use xlink:href='#icon-database'></use></svg></label><div class='xr-var-attrs'><dl class='xr-attrs'><dt><span>grid_mapping :</span></dt><dd>crs</dd><dt><span>proj4 :</span></dt><dd>+proj=longlat +datum=WGS84 +no_defs</dd><dt><span>min :</span></dt><dd>[1. 1. 1. 1. 1. 1. 1. 1. 1. 1. 1. 1. 1. 1. 1. 1. 1. 1. 1.]</dd><dt><span>max :</span></dt><dd>[1. 1. 1. 1. 1. 1. 1. 1. 1. 1. 1. 1. 1. 1. 1. 1. 1. 1. 1.]</dd></dl></div><div class='xr-var-data'><pre>[1231200 values with dtype=float32]</pre></div></li></ul></div></li><li class='xr-section-item'><input id='section-6b5df638-e854-42e1-b454-11590e18dbba' class='xr-section-summary-in' type='checkbox'  checked><label for='section-6b5df638-e854-42e1-b454-11590e18dbba' class='xr-section-summary' >Attributes: <span>(3)</span></label><div class='xr-section-inline-details'></div><div class='xr-section-details'><dl class='xr-attrs'><dt><span>Conventions :</span></dt><dd>CF-1.4</dd><dt><span>created_by :</span></dt><dd>R, packages ncdf4 and raster (version 3.5-29)</dd><dt><span>date :</span></dt><dd>2022-12-15 19:36:13</dd></dl></div></li></ul></div></div>

Here we can see that the FAO region names are not saved correctly. They
are numbered based on its location on the shapefile. We can update this
using the data frame with the unique FAO regions we previously created.

``` python
#We can simply load the data frame from the R environment.
FAO_keys = r.fao_names_codes

#Ensuring ID column is read as integer
FAO_keys['F_CODE'] = FAO_keys['F_CODE'].astype('int')

#We could also load it from our disk using the line below.
#FAO_keys = pd.read_csv("../Spatial_Data/FAO_shapefiles/FAO-regions_keys.csv")

#Checking list
FAO_keys
```

    ##     F_CODE                      NAME_EN
    ## 0       18                   Arctic_Sea
    ## 1       21           Atlantic_Northwest
    ## 2       27           Atlantic_Northeast
    ## 3       31     Atlantic_Western_Central
    ## 4       34     Atlantic_Eastern_Central
    ## 5       37  Mediterranean_and_Black_Sea
    ## 6       41           Atlantic_Southwest
    ## 7       47           Atlantic_Southeast
    ## 8       48           Atlantic_Antarctic
    ## 9       51         Indian_Ocean_Western
    ## 10      57         Indian_Ocean_Eastern
    ## 11      58       Indian_Ocean_Antarctic
    ## 12      61            Pacific_Northwest
    ## 13      67            Pacific_Northeast
    ## 14      71      Pacific_Western_Central
    ## 15      77      Pacific_Eastern_Central
    ## 16      81            Pacific_Southwest
    ## 17      87            Pacific_Southeast
    ## 18      88            Pacific_Antarctic

We can now update the names on the `netcdf` file.

``` python
#We will use the values in the data frame we loaded above
mask['FAO_reg_name'] = FAO_keys.F_CODE.tolist()

#Checking results
mask = mask.rename({'latitude': 'lat', 'longitude': 'lon'})
```

## Plotting results

We will plot all regions below to ensure we got them all correctly.

``` python
#We will loop through each layer
for reg in mask.FAO_regions:
  #Plotting results
  fig = plt.figure()
  ax = fig.add_subplot(111)
  reg.plot(ax = ax, levels = [1, 2])
  plt.title(reg.FAO_reg_name.values.tolist())
  plt.show()
```

<img src="Creating_Mask_FAO_regions_files/figure-gfm/unnamed-chunk-8-1.png" width="672" /><img src="Creating_Mask_FAO_regions_files/figure-gfm/unnamed-chunk-8-2.png" width="672" /><img src="Creating_Mask_FAO_regions_files/figure-gfm/unnamed-chunk-8-3.png" width="672" /><img src="Creating_Mask_FAO_regions_files/figure-gfm/unnamed-chunk-8-4.png" width="672" /><img src="Creating_Mask_FAO_regions_files/figure-gfm/unnamed-chunk-8-5.png" width="672" /><img src="Creating_Mask_FAO_regions_files/figure-gfm/unnamed-chunk-8-6.png" width="672" /><img src="Creating_Mask_FAO_regions_files/figure-gfm/unnamed-chunk-8-7.png" width="672" /><img src="Creating_Mask_FAO_regions_files/figure-gfm/unnamed-chunk-8-8.png" width="672" /><img src="Creating_Mask_FAO_regions_files/figure-gfm/unnamed-chunk-8-9.png" width="672" /><img src="Creating_Mask_FAO_regions_files/figure-gfm/unnamed-chunk-8-10.png" width="672" /><img src="Creating_Mask_FAO_regions_files/figure-gfm/unnamed-chunk-8-11.png" width="672" /><img src="Creating_Mask_FAO_regions_files/figure-gfm/unnamed-chunk-8-12.png" width="672" /><img src="Creating_Mask_FAO_regions_files/figure-gfm/unnamed-chunk-8-13.png" width="672" /><img src="Creating_Mask_FAO_regions_files/figure-gfm/unnamed-chunk-8-14.png" width="672" /><img src="Creating_Mask_FAO_regions_files/figure-gfm/unnamed-chunk-8-15.png" width="672" /><img src="Creating_Mask_FAO_regions_files/figure-gfm/unnamed-chunk-8-16.png" width="672" /><img src="Creating_Mask_FAO_regions_files/figure-gfm/unnamed-chunk-8-17.png" width="672" /><img src="Creating_Mask_FAO_regions_files/figure-gfm/unnamed-chunk-8-18.png" width="672" /><img src="Creating_Mask_FAO_regions_files/figure-gfm/unnamed-chunk-8-19.png" width="672" /><img src="Creating_Mask_FAO_regions_files/figure-gfm/unnamed-chunk-8-20.png" width="672" />
When comparing to the shapefile plots in the `R` section of this
notebook, we can see that the regions are named correctly. This means
that we can save our results now.

``` python
#Creating filename
fn = f'FAO-regions-corrected_{r.res}mask{r.mod_name}.nc'

#Saving result
mask.to_netcdf(os.path.join('../Spatial_Data/Masks', fn))
```
