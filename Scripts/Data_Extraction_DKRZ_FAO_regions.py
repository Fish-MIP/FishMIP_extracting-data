#!/usr/bin/python

#Libraries
import numpy as np
import xarray as xr
import pandas as pd
from glob import glob
import os
import re

#######################################################################################
#Variables between the hash lines can be edited
#Variable of interest - as it appears in the models
var_int = 'tcb'
#Keywords used to identified the files that will be processed
#These keywords must be present in all files across all models
file_key = '*nat_default_tcb_g*.nc'

#Base directory where outputs will be saved
base_out = 'FAO_data_extractions'

#Base directory where data is currently stored
base_dir = '/work/bb0820/ISIMIP/ISIMIP3b/OutputData/marine-fishery_global/'

#Indicate location of EEZ masks and area rasters
eez_mask_all = xr.open_dataset('Masks/EEZ-world-corrected_1degmask.nc')
area_all = xr.open_dataset('Masks/area_1deg.nc').area
#Masks for DBPM model
eez_mask_DBPM = xr.open_dataset('Masks/EEZ-world-corrected_1degmask_DBPM.nc')
area_DBPM = xr.open_dataset('Masks/area_1deg_DBPM.nc').area
#Masks for DBEM model
eez_mask_DBEM = xr.open_dataset('Masks/EEZ-world-corrected_05degmask.nc')
area_DBEM = xr.open_dataset('Masks/area_05deg.nc').area
#######################################################################################


#######################################################################################
#The section below will use the input above to find datasets of interest and calculate
#weighted means per year and sector.

#Ensuring base directory exists
os.makedirs(base_out, exist_ok = True)

#Go through each model/esm/activity and find netcdf files for variable of interest
file_key = f'*/*/*/{file_key}'
file_list = glob(os.path.join(base_dir, file_key))
#Removing any files for "picontrol" activity
file_list = [f for f in file_list if "picontrol" not in f]

#Saving names of experiments and ESMs to save results of work in the same directory structure
dir_str = []
for exp in os.listdir(base_dir):
    dir_list = [os.path.join(exp, e) for e in os.listdir(os.path.join(base_dir, exp))]
    for esm in dir_list:
        dir_list = [os.path.join(esm, a) for a in os.listdir(os.path.join(base_dir, esm))]
        for d in dir_list:
            dir_str.append(d)

#Loading future projections
def load_ds_noncf(fn, start, end):
    ds = xr.open_dataset(fn, decode_times = False)
    #Get start and end years for projections
    years = (end-start)+1
    if len(ds.time)/years == 1:
        freq = 'YS'
    elif len(ds.time)/years == 12:
        freq = 'MS'
    ds['time'] = pd.date_range(f'{start}-01-01', periods = len(ds.time), freq = freq)
    return ds

#Defining function to calculate weighted means
def weighted_means(ds, weight, mask):
    yr_means = []
    yr_anom = []
    for eez in mask.EEZ_regions:
        #Applying EEZ mask
        masked_ds = ds*eez
        #Calculate weights using grid cell area
        weights = weight*eez
        weights = weights/weights.sum()
        #Save results in list
        yr_mean_sec = (masked_ds*weights).groupby('time.year').sum(('lon', 'lat'))
        ref = yr_mean_sec.sel(time = slice('1990', '1999')).mean('time')
        yr_change = ((yr_mean_sec-ref)/ref)*100
        yr_means.append(yr_mean_sec)
        yr_anom.append(yr_change)
    yr_means = xr.concat(yr_means, dim = 'Country').to_dataset('Country').to_dataframe()
    yr_anom = xr.concat(yr_anom, dim = 'Country').to_dataset('Country').to_dataframe()
    return yr_means, yr_anom

#Getting list of historical and future projection experiments
file_hist = [f for f in file_list if "historical" in f]
file_non_hist = [f for f in file_list if "historical" not in f]
#Looping through list of files
for f in file_hist:
    #Find the correct folder to store files
    dir_out =  [d for d in dir_str if d in f]
    #Get the model and ESM to find the correct projection files
    exp, esm = re.split("/", dir_out[0])[:-1]
    print(exp, esm)
    future_paths = [d for d in file_non_hist if exp in d]
    future_paths = [d for d in future_paths if esm in d]
    dir_out_future = [d for d in dir_str if d in future_paths[0]]
    if len(dir_out) > 1:
        print("check the output directory folder")
    elif len(future_paths) > 2:
        print("wrong number of future projections")
    else:
        path_out = os.path.join(base_out, dir_out[0])
        path_out_future = os.path.join(base_out, dir_out_future[0])
        #Ensure folder exists
        os.makedirs(path_out, exist_ok = True)
        os.makedirs(path_out_future, exist_ok = True)
        #Extracting base file name to create output
        base_file = re.split("global_", re.split("/", f)[-1])[0]
        base_file_126 = re.split("global_", re.split("/", [d for d in future_paths if 'ssp126' in d][0])[-1])[0]
        base_file_585 = re.split("global_", re.split("/", [d for d in future_paths if 'ssp585' in d][0])[-1])[0]
        #Loading datasets
        #Historical
        #Get start and end years for data
        yr_min_hist, yr_max_hist = re.split("_", re.findall("\d{4}_\d{4}", re.split("/", f)[-1])[0])
        try:
            ds = xr.open_dataset(f).sel(time = slice('1950', '2015'))
        except:
            print('Time in historical data is not cf compliant. Fixing dates based on years in file name.')
            try:
                ds = load_ds_noncf(f, int(yr_min_hist), int(yr_max_hist)).sel(time = slice('1950', '2015'))
            except:
                print(f'{f} could not be opened.')
        #Get start and end years for projections
        yr_min_fut, yr_max_fut = re.split("_", re.findall("\d{4}_\d{4}", re.split("/", future_paths[0])[-1])[0])
        ds_126 = load_ds_noncf([d for d in future_paths if 'ssp126' in d][0], int(yr_min_fut), int(yr_max_fut))
        ds_585 = load_ds_noncf([d for d in future_paths if 'ssp585' in d][0], int(yr_min_fut), int(yr_max_fut))
        #Ensure flag values 1e20 are masked
        if (~np.isfinite(ds[var_int])).sum() == 0:
            ds = ds.where(ds < 1e20)
            ds_126 = ds_126.where(ds_126 < 1e20)
            ds_585 = ds_585.where(ds_585 < 1e20)
        #Load the correct grid area and mask rasters that match the model
        if (exp.lower() == "dbpm") or ('ipsl' in esm.lower() and exp.lower() == 'zoomss'):
            eez_mask = eez_mask_DBPM
            area = area_DBPM
        elif exp.lower() == "dbem":
            eez_mask = eez_mask_DBEM
            area = area_DBEM
        else:
            eez_mask = eez_mask_all
            area = area_all
        #Calculating mean weighted yearly values per EEZ
        ds_yr_mean, ds_anom_per = weighted_means(ds[var_int], area, eez_mask)
        ds_yr_mean_126, ds_anom_per_126 = weighted_means(ds_126[var_int], area, eez_mask)
        ds_yr_mean_585, ds_anom_per_585 = weighted_means(ds_585[var_int], area, eez_mask)
        #Getting output path
        mean_yr_out = os.path.join(path_out, (f'{base_file}_global_weighted_mean_abs_{yr_min_hist}_{yr_max_hist}.csv'))
        anom_out_per = os.path.join(path_out, (f'{base_file}_global_weighted_mean_per_{yr_min_hist}_{yr_max_hist}.csv'))
        mean_yr_out_126 = os.path.join(path_out_future, (f'{base_file_126}_global_weighted_mean_abs_{yr_min_fut}_{yr_max_fut}.csv'))
        anom_out_per_126 = os.path.join(path_out_future, (f'{base_file_126}_global_weighted_mean_per_{yr_min_fut}_{yr_max_fut}.csv'))
        mean_yr_out_585 = os.path.join(path_out_future, (f'{base_file_585}_global_weighted_mean_abs_{yr_min_fut}_{yr_max_fut}.csv'))
        anom_out_per_585 = os.path.join(path_out_future, (f'{base_file_585}_global_weighted_mean_per_{yr_min_fut}_{yr_max_fut}.csv'))
        #Saving outputs
        ds_yr_mean.to_csv(mean_yr_out, na_rep = np.nan)
        ds_anom_per.to_csv(anom_out_per, na_rep = np.nan)
        ds_yr_mean_126.to_csv(mean_yr_out_126, na_rep = np.nan)
        ds_anom_per_126.to_csv(anom_out_per_126, na_rep = np.nan)
        ds_yr_mean_585.to_csv(mean_yr_out_585, na_rep = np.nan)
        ds_anom_per_585.to_csv(anom_out_per_585, na_rep = np.nan)
        

