# -*- coding: utf-8 -*-
"""
Calculate Cd for Paper2_optimizingRestoration
using method in MacDonald et al., 2006

BKN - USGS 2022
"""

import pickle
import os
import sys
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import signal
from pathlib import Path
sys.path.insert(1, 'c:\\Users\\bknorris\\Documents\\Scripts\\Postdoc\\Paper2_OptimizingRestoration\\DataAnalysis\\loadModelFiles\\')
import analysisUtils 
from processModels import processModels

# Plot commands -- run this at the top of plotting routines
plt.close('all')
matplotlib.rcParams['mathtext.fontset'] = 'custom'
matplotlib.rcParams['mathtext.it'] = 'Arial:italic'
matplotlib.rcParams['mathtext.rm'] = 'Arial'
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42

# Define file paths:
csv_source = Path('c:/Users/bknorris/Documents/Models/Paper2_OptimizingRestoration/ModelRuns/Scenarios/')
model_source = Path('c:/Users/bknorris/Documents/Models/Paper2_OptimizingRestoration/ModelRuns/Scenarios/postProcessed')
save_fig_dir = Path('c:/Users/bknorris/Documents/Models/Paper2_OptimizingRestoration/Figures')
csv_file = 'modelPostProcessing_mod1.csv'
version_no = 'V3'

save_figures = False

# Load CSV
csv_file = csv_source / csv_file
model_info = pd.read_csv(csv_file.open())

f1 = plt.figure(figsize=(6, 5))
ax = f1.add_axes([0.1, 0.1, 0.7, 0.8])

# Load binary results file
filenames = next(os.walk(str(model_source)))[2]  # For later, [1] is dirs and [2] is files
for idx, scenario in enumerate(model_info['orgScenarioNumber']):
    # Check to see if models have been processed already; unzip if no; load if yes
    model_files = analysisUtils.find_files(model_source, f'Scenario_{scenario}_*_{version_no}*', 'files')
     
    # Calculate ub
    file = open(model_source / model_files[3], 'rb')
    freeSurf = pickle.load(file)
    file.close()
    ub = []
    for gauge, data in freeSurf.iteritems():
        if gauge != 'TimeStep':
            fs = 1 / 0.083  # sample frequency from models
            data = signal.detrend(data)
            seg = round(len(data) / 2, 0)
            noverlap = seg / 2
            f, Pxx_eta = signal.welch(data, fs,
                                      window='hanning',
                                      nperseg=seg,
                                      nfft=seg,
                                      detrend=False,
                                      scaling='spectrum',
                                      noverlap=noverlap)
            kh = analysisUtils.qkhf(f, 0.028)
            cutoff = np.argmax(f >= 2)
            
            # Estimate near-bottom orbital velocity with LWT
            # Lowe et al., 2007; Henderson et al., 2017
            kh = analysisUtils.qkhf(f, 0.028)
            cutoff = np.argmax(f >= 2)
            c = np.sqrt(9.81 * np.tanh(kh[1:cutoff])) / np.sqrt(kh[1:cutoff] / 0.028)
            n = (1 / 2) * (1 + ((2 * kh[1:cutoff]) / (np.sinh(2 * kh[1:cutoff]))))
            cg = c * n  # Wave celerity 
            F = 1025 * 9.81 * Pxx_eta[1:cutoff] * cg  # Wave energy spectrum
            urms = np.sqrt(2 * np.trapz(F) * np.mean(np.diff(f[1:cutoff]))) # Lowe eq 10
            ub.append(((8 / np.pi)**(1/2))*urms) 
    
    # Compute free surface slope eta/x
    if model_info.wavePeriod[idx] == 5:
        S = np.abs(freeSurf['WG1'] - freeSurf['WG10']) / 2.6
    elif model_info.wavePeriod[idx] == 10:
        S = np.abs(freeSurf['WG1'] - freeSurf['WG10']) / 5.2
    else:
        S = np.abs(freeSurf['WG1'] - freeSurf['WG10']) / 10.4
    
    H = (np.mean(freeSurf['WG10']) + 0.027)
    
    Cd = ((2 * 9.81 * H * np.mean(S)) / (np.mean(ub)**2))
    a = (model_info.deltaS[idx]  * 36) / 0.4
    ax.plot(a, Cd, '.')
    # print(f'Cd: {Cd}')
    print(f'{Cd}')
    # print(f'{np.mean(S)}')
    # print(f'Ub: {np.mean(ub)}')