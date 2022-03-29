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
import processModels
from pathlib import Path


# Define file paths:
model_source = Path('h:/Models/Paper2_OptimizingRestoration/ModelRuns/Scenarios/rawResults')
model_dest = Path('c:/Users/bknorris/Documents/Models/Paper2_OptimizingRestoration/ModelRuns/Scenarios/')
save_data_path = str(model_dest) + "\\postProcessed\\"
model_info = 'adjustModelSampling_TEST.csv'

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
#t1 = time.time()
for idx, scenario in enumerate(model_info['orgScenarioNumber']):
    #t2 = time.time()
    print(f'Loading file {idx+1} of {len(model_info)}: Scenario_{scenario}\n')
    
    # Get the scenarios from the CSV file in the model directory
    regex = re.compile('Scenario_' + str(scenario) + '_.*')
    bin_files = list(filter(regex.match, filenames))
    model_files = [bin_files, idx]
    
    # Model post processing steps
    # NEED TO ADD CANOPY PARAMS INTO CSV FILE!
    
    