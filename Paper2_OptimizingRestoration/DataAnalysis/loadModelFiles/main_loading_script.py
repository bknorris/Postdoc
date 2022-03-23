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
import time
import zipfile
import shutil
import numpy as np
import pandas as pd
import loadModelFiles
from pathlib import Path

version_no = 'V1'  # file version; appends to model save file

# Define file paths:
model_source = Path('h:/Models/Paper2_OptimizingRestoration/ModelRuns/Scenarios/Reprocessed/')
model_dest = Path('c:/Users/bknorris/Documents/Models/Paper2_OptimizingRestoration/ModelRuns/Scenarios/')
save_data_path = str(model_dest) + "\\DataAnalysis\\"
model_info = 'adjustModelSampling_RESTART2.csv'

# Begin program
# 1. Load CSV
t = time.localtime()
current = time.strftime("%H:%M:%S", t)
print('The current time is: ' + current)
csv_file = model_dest / model_info
model_info = pd.read_csv(csv_file.open())

# 2. Get model names from the model directory
filenames = next(os.walk(str(model_source)))[2]  # For later, [1] is dirs and [2] is files

# 3. Load files in model_source directory
t1 = time.time()
for file in filenames:
    t2 = time.time()
    
    # Only process ZIP files!
    if '.zip' in file:
        # 3a. Check if file is in model_info, only process if so
        scenario_number = re.search('\_(.*)\.', file).group(1)
        idx = np.where(int(scenario_number) == model_info['orgScenarioNumber'])
        if np.any(idx):
            idx = int(idx[0])
            print('Processing ' + file)
            
            # 3b. Unzip model file to model_source
            print('Unzipping file from source directory...')
            model_name = file.split('.')[0]  # model folder without extension
            with zipfile.ZipFile(str(model_source / file), 'r') as zip_file:
                zip_file.extractall(str(model_dest / model_name))
            
            # 3c. Run loading functions on model files
            model = loadModelFiles.loadModelFiles(str(model_dest), model_name, model_info['wavePeriod'][idx])
            timestep = model.createTimestep()
            freeSurf = model.loadFreeSurface(timestep)
            fields = model.loadFields(timestep)
            
            # 3d. Save out binary files to disk
            model.saveData(freeSurf, 'freeSurf', save_data_path, version_no)
            model.saveData(fields, 'fields', save_data_path, version_no)
            
            # 3e. Clear memory, remove model folder on disk to preserve space
            print('Removing model folder from dest directory...')
            del model, timestep, freeSurf, fields
            full_path = str(model_dest / model_name)
            shutil.rmtree(full_path)
            
            # Timer
            executionTime = (time.time() - t2) / 60
            print(f'File {file} processed in: {executionTime:.2f} minutes\n')

# End timer
executionTime = (time.time() - t1) / 60
print(f'Total runtime: {executionTime:.2f} minutes\n')
        
