#Libraries
import matplotlib.pyplot as plt
import requests
import xarray as xr
import json
from pathlib import Path
import time
import zipfile

#Accessing FishMIP data
response = requests.get('https://data.isimip.org/api/v1/datasets', params={
    'simulation_round': 'ISIMIP3b',
    'product': 'OutputData',
    'climate_forcing': 'ipsl-cm6a-lr',
    'climate_scenario': 'historical',
    'sector': 'marine-fishery_global',
    'model': 'apecosm',
    'variable': 'tcb'
})

#Connecting to API
response.raise_for_status()
response_data = response.json()
dataset = response_data['results'][0]
#This dataset contains one file only
paths = [file['path'] for file in dataset['files']]


#Cut out area and create mask based 
operations = [
    #first, cut-out around Prydz Bay
    {
        'operation': 'cutout_bbox',
        'bbox': [
             60,  # west
             90,  # east
            -70,  # south
            -55   # north
        ]
    },
    # mask the file using the created mask
    {
        'operation': 'mask_mask',
        'mask': 'test_da.nc',
	      #even if 'var' is commented out, the request does not work
	      'var': 'region'
    }
    ]

#Connect to API again to create mask
response = requests.post('https://files.isimip.org/api/v2', files={
    'data': json.dumps({
        'paths': paths,
        'operations': operations
    }),
    'test_da.nc': Path('test_da.nc').read_bytes(),  # needs to be the same as in the create_mask operation
})
response.raise_for_status()
job = response.json() #Status is 'failed' when using a mask not produced by/modified from API in operations

#Send request to API
while True:
  job = requests.get(job['job_url']).json()
  print(job['status'], job['meta'])

  time.sleep(10)
  if job['status'] not in ['queued', 'started']:
      break
      
#Download request
zip_path = Path(job['file_name'])
with requests.get(job['file_url'], stream=True) as response:
    with zip_path.open('wb') as fp:
        for chunk in response.iter_content(chunk_size=8192):
            fp.write(chunk)

#Unzip data
out_path = zip_path.with_suffix('')
with zipfile.ZipFile(zip_path, 'r') as zf:
    zf.extractall(out_path)

da = xr.open_dataarray('isimip-download-a8c74c0eac5f7885d711ab6f35e1508dd0273af9/apecosm_ipsl-cm6a-lr_nobasd_historical_nat_default_tcb_lon60.0to90.0lat-70.0to-55.0+test_da_monthly_1850_2014.nc')

da.isel(time = 0).plot()
plt.show()
