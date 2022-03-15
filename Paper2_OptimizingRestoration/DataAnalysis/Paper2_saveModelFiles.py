# -*- coding: utf-8 -*-
"""
Load and save Paper2_OptimizingRestoration model data
from the postProcessing folders in each model folder and
save the data out to disk.

This script is a "Pythonic" version of a MATLAB script I wrote
last year.

BKN - USGS PCMSC 2022
"""
import os
import re
import zipfile
from pathlib import Path
import pandas as pd

version_no = 'V1'  # file version; appends to model save file

# Define file paths:
model_source = Path('h:/Models/Paper2_OptimizingRestoration/ModelRuns/Scenarios/Reprocessed/')
model_dest = Path('c:/Users/bknorris/Documents/Models/Paper2_OptimizingRestoration/ModelRuns/Scenarios/')
model_info = 'adjustModelSampling_TEST.csv'

# Define model default settings ### MOVE THIS TO UTILITIES SCRIPT!
x_bnd = (0, 2.6)  # lower/upper
z_bnd = (-5, -0.005)  # lower/upper
model_timestep = 0.5  # timestep in s

# Begin program
# 1. Load CSV
csv_file = model_dest / model_info
model_info = pd.read_csv(csv_file.open())

# 2. Get model names from the model directory
filenames = next(os.walk(str(model_source)), (None, None, []))[2]

# 3. Run data processing on files in model_source
for idx, file in enumerate(filenames):
    # 3a. Check if file is in model_info, only process if so
    scenario_number = re.search('\_(.*)\.', file).group(1)
    if int(scenario_number) in model_info['orgScenarioNumber']:
        
        print('Processing ' + file)
        
        # 3b. Unzip model file to model_source
        with zipfile.ZipFile(str(model_source / file), 'r') as zip_file:
            zip_file.extractall(str(model_dest / file.split('.')[0]))
        
        # 3c. Run processing functions on model files
        