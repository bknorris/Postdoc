# -*- coding: utf-8 -*-
"""
Plot different wave attenuation parameters to try and work out what's going on with the data...

BKN - USGS 2022
"""

import pickle
import os
import sys
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import signal
from pathlib import Path
sys.path.insert(1, 'c:\\Users\\bknorris\\Documents\\Scripts\\Postdoc\\Paper2_OptimizingRestoration\\DataAnalysis\\loadModelFiles\\')
import analysisUtils 
from processModels import processModels

plt.close('all')

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

# Load binary results file
filenames = next(os.walk(str(model_source)))[2]  # For later, [1] is dirs and [2] is files
for idx, scenario in enumerate(model_info['orgScenarioNumber']):
    # Check to see if models have been processed already; unzip if no; load if yes
    model_files = analysisUtils.find_files(model_source, f'Scenario_{scenario}_*_{version_no}*', 'files')
        
    file = open(model_source / model_files[3], 'rb')
    freeSurf = pickle.load(file)
    file.close()
    
    # Parameters
    Cd = 0.24
    D = 0.011
    h = 0.028 
    s = 0.0056 
    deltaS = model_info.deltaS[idx]
    a = (D / deltaS)
    N = (a / D)
    alpha = []
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
                                      scaling='density',
                                      noverlap=noverlap)

            cutoff = np.argmax(f >= 2)
            kh = analysisUtils.qkhf(f, h)
            m0 = Pxx_eta[1:cutoff] * np.mean(np.diff(f))
            Hrms = np.sqrt(8 * m0)
            k = kh[1:cutoff] / h
            kh = kh[1:cutoff]
            ks = k * s
            term1 = (1 / (3 * np.sqrt(np.pi))) * Cd * D * N * Hrms * k
            term2 = ((np.sinh(ks)**3) + 3 * np.sinh(ks)) / ((np.sinh(2 * kh) + 2 * kh) * np.sinh(kh))
            alpha.append(np.trapz(36 * term1 * term2))
    
    if idx == 0:
        a_dict = dict()
        a_dict[str(scenario)] = alpha
    else:
        a_dict[str(scenario)] = alpha
        
        
# Get data indices from CSV -- filter by waveHeight
five = np.where(model_info.waveHeight == 0.001388889)[0].tolist()
ten = np.where(model_info.waveHeight == 0.004166667)[0].tolist()
twenty = np.where(model_info.waveHeight == 0.008333333)[0].tolist()
waves = [five, ten, twenty]

# Begin figure
plt.style.use('_mpl-gallery')
f1 = plt.figure(figsize=(8, 10.5))
axes1 = f1.add_axes([0.1, 0.2, 0.25, 0.22])
axes2 = f1.add_axes([0.4, 0.2, 0.25, 0.22])
axes3 = f1.add_axes([0.7, 0.2, 0.25, 0.22])
axes4 = f1.add_axes([0.1, 0.47, 0.25, 0.22])
axes5 = f1.add_axes([0.4, 0.47, 0.25, 0.22])
axes6 = f1.add_axes([0.7, 0.47, 0.25, 0.22])
axes7 = f1.add_axes([0.1, 0.74, 0.25, 0.22])
axes8 = f1.add_axes([0.4, 0.74, 0.25, 0.22])
axes9 = f1.add_axes([0.7, 0.74, 0.25, 0.22])
ax = [[axes1, axes2, axes3],
      [axes4, axes5, axes6],
      [axes7, axes8, axes9]]

# Colormaps from ColorBrewer
BuGn = ['#ccece6', '#93ccb1', '#63aa83', '#3c885b', '#1b6539', '#00441b']
OrRd = ['#fdd49e', '#f5a677', '#e27b55', '#c65336', '#a42e1b', '#7f0000']
PuBu = ['#d0d1e6', '#99b3d2', '#6993b8', '#40749b', '#1d567a', '#023858']
cmap = [BuGn, OrRd, PuBu]

# Marker spec
markers = ['o', '^', 's', 'd', 'X', '>']

# Plot routine
xs = np.linspace(-0.025, 0.025, 5)
wave_gauges = np.linspace(0, 1, 10)

for i in range(0, len(waves)):
    scenario = model_info.orgScenarioNumber[waves[i]].tolist()
    scenario = [str(x) for x in scenario]
    wavePeriod = model_info.wavePeriod[waves[i]].tolist()
    wave_unq = np.unique(wavePeriod)
    spacing = model_info.deltaS[waves[i]]
    spce_unq = np.unique(spacing)
    
    for j in range(0, len(wave_unq)):
        for k in range(0, len(spce_unq)):
            idx1 = np.where(wavePeriod == wave_unq[j])[0].tolist()
            idx2 = np.where(spacing == spce_unq[k])[0].tolist()
            idx = np.intersect1d(idx1, idx2)[0]
        
            wave = a_dict[scenario[idx]]
        
            x = wave_gauges
            y = wave
            
            x_scaled = (spce_unq[k] * 36) / 0.4
            leg = f'${x_scaled:.0f}D$'
            ax[i][j].plot(x[1:] + xs[j], y[1:], color=cmap[j][k],
                          marker=markers[k],
                          markersize=5,
                          label=leg)

# Plot Adjustments:
# Axis scaling
for i in range(0, 3):
    [ax[i][j].set_xlim(0, 1.03) for j in range(0, 3)]
    [ax[i][j].set_ylim(10e-5, 10e-2) for j in range(0, 3)]
