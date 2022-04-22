# -*- coding: utf-8 -*-
"""
Load and process Paper2_OptimizingRestoration model data
from the postProcessing folders in each model folder and
save the data out to disk. This script combines functionality
from main_loading_script.py and main_processing_script.py
as some of the model RAW files are too large.


BKN - USGS PCMSC 2022
"""
import os
import re
import time
import zipfile
import shutil
import pickle
import pandas as pd
import loadModelFiles
import processModels
from analysisUtils import find_files
from pathlib import Path

version_no = 'V2'  # file version; appends to model save file
comment = 'Tw_5_10_20s'  # Output file comment

# Define file paths:
model_source = Path('e:/BKN-FIELD/Models/Paper2_OptimizingRestoration/ModelRuns/Scenarios/')
model_dest = Path('c:/Users/bknorris/Documents/Models/Paper2_OptimizingRestoration/ModelRuns/Scenarios/')
save_data_path = str(model_dest) + "\\postProcessed\\"
model_info = 'modelPostProcessing.csv'

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
for idx, scenario in enumerate(model_info['orgScenarioNumber']):
    t2 = time.time()
    print(f'Loading file {idx+1} of {len(model_info)}: Scenario_{scenario}')
    
    # Get the scenarios from the CSV file in the model directory
    regex = re.compile('Scenario_' + str(scenario) + '.zip')
    file = list(filter(regex.match, filenames))
    
    # Check to see if models have been processed already; unzip if no; load if yes
    model_files = find_files(save_data_path, f'Scenario_{scenario}_*_{version_no}*', 'files')
    if not model_files:  # list is empty
        # Unzip model file to model_source
        print('Unzipping file from source directory...')
        model_name = file[0].split('.')[0]  # model folder without extension
        with zipfile.ZipFile(str(model_source / file[0]), 'r') as zip_file:
            zip_file.extractall(str(model_dest / model_name))
    
        # Run loading functions on model files
        model_load = loadModelFiles.loadModelFiles(str(model_dest), model_name, model_info['wavePeriod'][idx])
        timestep = model_load.createTimestep()
        freeSurf = model_load.loadFreeSurface(timestep)
        fields = model_load.loadFields(timestep)
        
        # Model post processing steps
        model_pprocess = processModels.processModels(model_source, fields, freeSurf, model_info, idx)
        wef, Ab, fer, ubr = model_pprocess.computeWaveEnergy()
        avg_fields = model_pprocess.avgFields()
        eps_norm = model_pprocess.feddersenDissipation(avg_fields['eps'], wef['epsj'][:-1])
        Lambda_not = model_pprocess.hendersonDamping(ubr)
        calc_values = {'Ab': Ab, 'fer': fer, 'ubr': ubr, 'lambda_not': Lambda_not}
        
        # Save out binary files to disk
        model_load.saveData(freeSurf, 'freeSurf', save_data_path, version_no)
        model_load.saveData(wef, 'WEF', save_data_path, version_no)
        model_load.saveData(avg_fields, 'AVG-FIELDS', save_data_path, version_no)
        model_load.saveData(eps_norm, 'EPS-NORM', save_data_path, version_no)
        model_load.saveData(calc_values, 'CALC-VALUES', save_data_path, version_no)
        
        # Clear memory, remove model folder on disk to preserve space
        print('Removing model folder from dest directory...')
        del model_load, model_pprocess, timestep, freeSurf, fields
        full_path = str(model_dest / model_name)
        shutil.rmtree(full_path)
        
    elif model_files:  # load already processed binary files
        print('This model has been processed already...')
        file = open(save_data_path + model_files[0], 'rb')
        print(f'Reading {model_files[0]}')
        avg_fields = pickle.load(file)
        file.close()
        
        file = open(save_data_path + model_files[1], 'rb')
        print(f'Reading {model_files[1]}')
        calc_values = pickle.load(file)
        file.close()
        
        file = open(save_data_path + model_files[2], 'rb')
        print(f'Reading {model_files[2]}')
        eps_norm = pickle.load(file)
        file.close()
        
        file = open(save_data_path + model_files[3], 'rb')
        print(f'Reading {model_files[3]}')
        freeSurf = pickle.load(file)
        file.close()
        
        file = open(save_data_path + model_files[4], 'rb')
        print(f'Reading {model_files[4]}')
        wef = pickle.load(file)
        file.close()
        
    # Append data to dicts for each scenario num
    if idx == 0:
        wef_dict = dict()  # WEF values
        avf_dict = dict()  # avg_fields values
        epn_dict = dict()  # eps_norm values
        cal_dict = dict()  # other values
        wef_dict[str(scenario)] = wef
        avf_dict[str(scenario)] = avg_fields
        epn_dict[str(scenario)] = eps_norm
        cal_dict[str(scenario)] = calc_values
    else:
        wef_dict[str(scenario)] = wef
        avf_dict[str(scenario)] = avg_fields
        epn_dict[str(scenario)] = eps_norm
        cal_dict[str(scenario)] = calc_values
    
    # Print elapsed time
    executionTime = (time.time() - t2) / 60
    print(f'Completed in: {executionTime:.2f} minutes\n')
    
# Save data to file
print('Saving processed data...')
data = {'wef': wef_dict, 'avg_fields': avf_dict, 'eps_norm': epn_dict, 'calc_values': cal_dict}
save_file = False
output_file = 'modelPostProcessing_' + comment + '_' + version_no + ".dat"
if os.path.isfile(save_data_path + output_file):
    print('File already exists in destination directory!')
    result = input('Overwrite? (y/n) ')
    if result == 'y':
        save_file = True
        
else:
    save_file = True
    
if save_file:
    file_obj = open(save_data_path + output_file, mode='wb')
    pickle.dump(data, file_obj)
    file_obj.close()

# Print elapsed time
executionTime = (time.time() - t1) / 60
print(f'Data processed in: {executionTime:.2f} minutes')
