# Accessing and plotting data from the ISIMIP Repository
These collection of notebooks include a step by step guide on how to access data directly from the [Inter-Sectoral Impact Model Intercomparison Project (ISIMIP) Repository](https://data.isimip.org/) using the [`isimip-client`](https://github.com/ISI-MIP/isimip-client) library developed in Python by ISIMIP. These notebooks use the `reticulate` package to use the Python scripts directly in the R notebook.  

Notebooks showing examples on how to create simple maps making use of the `raster` and `ggplot2` packages are included in this repository.

## Installing `Python` libraries
To install all `Python` libraries used in this repository, use the `environment.yml` file. You can do this from the command line, make sure you navigate to the folder containing the `environment.yml` file, and then type the following line of code: `conda env create -f environment.yml -n CMIP6_data`. The `-f` argument points to the file containing all libraries to be installed, while the `-n` contains the name to be given to the newly created environment. You can change the name of this environment if you prefer, but make sure this name is updated in all scripts when calling Python. Look for this line of code: `use_condaenv("CMIP6_data")` and update the name of the environment if you chose to use a different name.  
  
You can check if the conda environment has been successfully created by typing this line of code in your command line: `conda env list`. The output will include the name and location of your newly created environment and it should look something similar to the following:  
  
```
# conda environments:
#
base                     /file_path_to_environment/base
CMIP6_data            *  /file_path_to_environment/envs/CMIP6_data
```
  
Once we have successfully install the environment following the steps above, we will use the command line once again and type the following commands:  
- `conda activate CMIP6_data` (if you used a different name, remember to change `CMIP6_data` for whatever you named your environment). This will activate the environment you created in the previous steps.  
- `which python` (for Linux or Mac) or `where python` (for Windows) to locate where `Python` has been installed.  
  
Open the file named `.Rprofile`, replace the folder path for the `RETICULATE_PYTHON` variable in line 2 with the `Python `folder path that was printed in the command line.