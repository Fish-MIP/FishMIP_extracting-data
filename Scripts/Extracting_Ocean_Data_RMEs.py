#!/usr/bin/python

#Extracting model data for regional models
#Author: Denisse Fierro Arcos
#2023-03-30
#This script runs in the DKRZ server. Ensure the Python
#module is loaded into your DKRZ session before running
#this script. This is done in the command line as follows:
#module load python3
#python3 Extracting_Ocean_Data_RMEs.py
#Note that the regional boundaries used here were created using
#the "Creating_Your_Own_Mask_From_Shapefiles" script available
#in this repository.

### Loading libraries
import xarray as xr
import pandas as pd
from glob import glob
import os
import re

### Region masks
#Loading masks containing the boundaries of the regions of interest
#We will be using a 1 deg resolution mask and a 0.25 deg resolution mask
mask_1deg = xr.open_dataarray(r'Masks/fishMIP_regional_1degmask_ISIMIP3a.nc')
mask_025deg = xr.open_dataarray(r'Masks/fishMIP_regional_025degmask_ISIMIP3a.nc')

#We will now rename the coordinates in our mask so they match the ocean datasets
mask_1deg = mask_1deg.rename({'Longitude': 'lon', 'Latitude': 'lat'})
mask_025deg = mask_025deg.rename({'Longitude': 'lon', 'Latitude': 'lat'})


### Ocean data
#We will define the base directory containing the ocean data of interest
base_dir = '/work/bb0820/ISIMIP/ISIMIP3a/InputData/climate/ocean/obsclim/global/monthly/historical/GFDL-MOM6-COBALT2'

#Finding files for variables of interest
var_int = ['phydiaz-vint', 'phypico-vint']

#Looping through each variable of interest
for var in var_int:
    #Searching files matching variables of interest
    in_file = glob(os.path.join(base_dir, f'*{var}*.nc'))
    #Extracting base for output file name
    out_file = [re.split("COBALT2/", re.split("_global",f)[0])[1] for f in in_file]
    #Looping through each file of interest
    for i, f in enumerate(in_file):
        #Open dataset
        data = xr.open_mfdataset(f)
        data = data[var]
        #Extracting minimum and maximum years to add to output filename
        yr_max = str(data.time.dt.year.max().values)
        yr_min = str(data.time.dt.year.min().values)
        #Selecting correct regional mask depending on resolution of dataset
        if '60arc' in f:
            mask = mask_1deg
        elif '15arc' in f:
            mask = mask_025deg
        #Looping through each region within mask
        for reg in mask.RME_name:
            reg_m = mask.sel(RME_name = reg)
            #Applying as mask to extract data
            d_reg = data.where(reg_m == 1, drop = True).to_series()
            #Reorganising data, so there is one column per time step and all NA values are removed
            if len(d_reg) > 0:
                d_reg = d_reg[~pd.isna(d_reg)].unstack('time')
                #Creating full output file name
                fn_out = f'{out_file[i]}_{str(reg.values)}_monthly_{yr_min}_{yr_max}.csv'
                #Saving final output in a folder called Regional
                os.makedirs('Regional', exist_ok = True)
                d_reg.to_csv(os.path.join('Regional', fn_out))