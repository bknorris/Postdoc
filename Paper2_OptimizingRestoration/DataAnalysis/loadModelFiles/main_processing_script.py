# -*- coding: utf-8 -*-
"""
Load and post process the RAW model data created by main_loading_script.py
and its dependencies. This script (and dependencies) creates reduced-size
data files to simplify figure making. This script computes key variables:
    
Total wave energy dissipation between the model inlet and patch;
Total mean TKE dissipation rate;
Gamma: wave breaking parameter;
Wavelength; and
Wave steepness.

TKE dissipation rate is binned (averaged) along the x-z axes into 100 bins to
reduce storage space.

This script is a "Pythonic" version of a MATLAB script I wrote
last year.

BKN - USGS PCMSC 2022
"""

import os
import re
import pandas as pd
import time
import pickle
import processModels
from pathlib import Path

version_no = 'V1'  # file version; appends to model save file
comment = '5s_10s'  # Output file comment

# Define file paths:
model_source = Path('h:/Models/Paper2_OptimizingRestoration/ModelRuns/Scenarios/rawResults')
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
    regex = re.compile('Scenario_' + str(scenario) + '_.*')
    bin_files = list(filter(regex.match, filenames))
    model_files = [bin_files, idx]
    
    # Model post processing steps
    model = processModels.processModels(model_source, model_files, model_info)
    wef, Ab, fer, ubr = model.computeWaveEnergy()
    avg_fields = model.avgFields()
    eps_norm = model.feddersenDissipation(avg_fields['eps'], wef['epsj'][:-1])
    Lambda_not = model.hendersonDamping(ubr)
    calc_values = {'Ab': Ab, 'fer': fer, 'ubr': ubr, 'lambda_not': Lambda_not}
    
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
