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
var_int = input('Write the name of the variable you want to process: ')
#Keywords used to identified the files that will be processed
#These keywords must be present in all files across all models
file_key = input('Write the common file pattern (e.g., *_default_tc_g*.nc): ')
#file_key = '*_default_tc_g*.nc'

#Base directory where outputs will be saved
base_out = 'FAO_data_extractions'

#Base directory where data is currently stored
base_dir = '/work/bb0820/ISIMIP/ISIMIP3b/OutputData/marine-fishery_global/'

#Indicate location of EEZ masks and area rasters and transform from km2 to m2
area_all = xr.open_dataarray('Masks/area_1deg_mask_FAO-EEZ.nc')*1e6
#Mask for DBPM model
area_DBPM = xr.open_dataarray('Masks/area_1deg_DBPM_mask_FAO-EEZ.nc')*1e6
#Mask for DBEM model
area_DBEM = xr.open_dataarray('Masks/area_05deg_mask_FAO-EEZ.nc')*1e6
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
#Getting list of experiments
for exp in os.listdir(base_dir):
    dir_list = [os.path.join(exp, e) for e in os.listdir(os.path.join(base_dir, exp))]
    #For each experiment get a list of ESMs
    for esm in dir_list:
        dir_list = [os.path.join(esm, a) for a in os.listdir(os.path.join(base_dir, esm))]
        #For each experiment and ESM combination, get a list of files
        for d in dir_list:
            dir_str.append(d)

#Loading data that is not CF compliant
def load_ds_noncf(fn, start, end):
    '''
    This function loads non-CF compliant datasets where dates cannot be read. It takes the following inputs:
    fn - ('string') refers to full filepath where the non-CF compliant dataset is located
    start - ('numeric') refers to the start year of the dataset
    end - ('numeric') refers to the end year of the dataset
    The start and end parameters are used to present dates correctly in the time dimension
    '''
    ds = xr.open_dataset(fn, decode_times = False)
    #Get start and end years for projections
    years = (end-start)+1
    if len(ds.time)/years == 1:
        freq = 'YS'
    elif len(ds.time)/years == 12:
        freq = 'MS'
    ds['time'] = pd.date_range(f'{start}-01-01', periods = len(ds.time), freq = freq)
    return ds

#Defining function to calculate weighted means and anomalies in relation to values in the 90s
def weighted_means(ds, weight, mask, file_out, **kwargs):
    '''
    This function calculates weighted means per year and yearly anomalies in relation to values 
    in the decade between 1990 and 1999. It takes the following inputs:
    ds - ('data array') refers to data array containing data upon which means will be calculated
    weight - ('data array') contains weights to be used in weighted mean calculation
    mask - ('data array') contains boundaries within which weighted means will be calculated
    file_out - ('string') contains the file path and base file name to be used to save results

    **Optional
    ref_ds - ('data array') containing the mean for the reference period. If non is provided, the
    function returns a data array containing this information.
    '''
    yr_means = []
    yr_anom = []
    if 'ref_ds' not in kwargs.keys():
        ref_ds = []
    else:
        ref_ds = kwargs.get('ref_ds')
    for eez in mask:
        #Masked grid cell area and calculate weights per EEZ
        weights = weight*eez
        weights = weights/weights.sum()
        #Save results in list
        yr_mean_sec = (ds*weights).groupby('time.year').sum(('lon', 'lat', 'time'))
        if 'ref_ds' in kwargs.keys():
            ref = ref_ds.sel(Country_EEZ = eez.Country_EEZ.values)
        else:
            ref = yr_mean_sec.sel(year = slice('1990', '1999')).mean('year')
        yr_change = ((yr_mean_sec-ref)/ref)*100
        yr_means.append(yr_mean_sec)
        yr_anom.append(yr_change)
        if 'ref_ds' not in kwargs.keys():
            ref_ds.append(ref)
    yr_means = xr.concat(yr_means, dim = 'Country_EEZ').to_dataset('Country_EEZ').to_dataframe()
    yr_anom = xr.concat(yr_anom, dim = 'Country_EEZ').to_dataset('Country_EEZ').to_dataframe()
    if 'ref_ds' not in kwargs.keys():
        ref_ds = xr.concat(ref_ds, dim = 'Country_EEZ')
    #Save results
    yr_min = str(yr_change.year.values.min())
    yr_max = str(yr_change.year.values.max())
    path_out_mean = f'{file_out}global_weighted_mean_abs_{yr_min}_{yr_max}.csv'
    path_out_anom = f'{file_out}global_weighted_mean_per_{yr_min}_{yr_max}.csv'
    yr_means.to_csv(path_out_mean, na_rep = np.nan)
    yr_anom.to_csv(path_out_anom, na_rep = np.nan)
    if 'ref_ds' in kwargs.keys():
        return yr_means, yr_anom
    else:
        return yr_means, yr_anom, ref_ds

#Defining function to calculate sum of values per FAO and EEZ area per month
def monthly_sum(ds, mask, file_out):
    '''
    This function calculates the sum of biomass per month per regions included in area data array.
    It takes the following inputs:
    ds - ('data array') refers to data array containing data upon which means will be calculated
    mask - ('data array') contains area per pixel and boundaries within which weighted means will 
    be calculated
    file_out - ('string') contains the file path and base file name to be used to save results
    '''
    #Creating empty array to store results
    month_sum = []
    #Multiplying dataset by area in m2
    ds_area = ds*mask
    #Calculating sums per year, per month and per region
    for yr, da in ds_area.groupby('time.year'):
        for mth, da_m in da.groupby('time.month'):
            month_sum.append(da_m.groupby('mask_FAO_EEZ').sum())
    #Create data array with results and transforming to tonnes
    month_sum = xr.concat(month_sum, dim = 'time')*1e-6
    #Saving results
    yr_min = str(ds.time.dt.year.values.min())
    yr_max = str(ds.time.dt.year.values.max())
    path_out = f'{file_out}global_tonnes_{yr_min}_{yr_max}.csv'
    month_sum.to_pandas().to_csv(path_out, na_rep = np.nan)

#Defining function to calculate range of values per EEZ area
def range_ds(ds, mask, file_out):
    '''
    This function calculates minimum and maximum values per year per EEZ area. It takes the following inputs:
    ds - ('data array') refers to data array containing data upon which means will be calculated
    mask - ('data array') contains boundaries within which weighted means will be calculated
    file_out - ('string') contains the file path and base file name to be used to save results
    '''
    min_vals = []
    max_vals = []
    for eez in mask:
        #Save results in list
        yr_min_sec = (ds*eez).groupby('time.year').min(('lon', 'lat', 'time'))
        yr_max_sec = (ds*eez).groupby('time.year').max(('lon', 'lat', 'time'))
        min_vals.append(yr_min_sec)
        max_vals.append(yr_max_sec)
    min_vals = xr.concat(min_vals, dim = 'Country_EEZ').to_dataset('Country_EEZ').to_dataframe()
    max_vals = xr.concat(max_vals, dim = 'Country_EEZ').to_dataset('Country_EEZ').to_dataframe()
    #Save results
    yr_min = str(yr_min_sec.year.values.min())
    yr_max = str(yr_min_sec.year.values.max())
    path_out_min = f'{file_out}global_min_{yr_min}_{yr_max}.csv'
    path_out_max = f'{file_out}global_max_{yr_min}_{yr_max}.csv'
    min_vals.to_csv(path_out_min, na_rep = np.nan)
    max_vals.to_csv(path_out_max, na_rep = np.nan)
    return min_vals, max_vals

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
        if exp.lower() == "feisty" and 'gfdl' in esm.lower():
            ds = ds.chunk({'time': 12})
            ds_126 = ds_126.chunk({'time': 12})
            ds_585 = ds_585.chunk({'time': 12})
        #Ensure flag values 1e20 are masked
        if (~np.isfinite(ds[var_int])).sum() == 0:
            ds = ds.where(ds < 1e20)
            ds_126 = ds_126.where(ds_126 < 1e20)
            ds_585 = ds_585.where(ds_585 < 1e20)
        #Load the correct grid area and mask rasters that match the model
        if (exp.lower() == "dbpm") or ('ipsl' in esm.lower() and exp.lower() == 'zoomss'):
            area = area_DBPM
        elif exp.lower() == "dbem":
            area = area_DBEM
        else:
            area = area_all
        #Calculating monthly sums per EEZ and FAO regions and saving to disk
        #Historical
        monthly_sum(ds[var_int], area, os.path.join(path_out, base_file))
        #SSP126
        monthly_sum(ds_126[var_int], area, os.path.join(path_out_future, base_file_126))
        #SSP585
        monthly_sum(ds_585[var_int], area, os.path.join(path_out_future, base_file_585))
