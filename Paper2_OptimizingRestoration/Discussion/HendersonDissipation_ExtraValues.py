# -*- coding: utf-8 -*-
"""
Calculate a greater range of water depth/wave height/period combinations
using the Henderson dissipation model

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
csv_file = 'modelPostProcessing_HendersonExtra.csv'
version_no = 'V3'

save_figures = False

# Load CSV
csv_file = csv_source / csv_file
model_info = pd.read_csv(csv_file.open())

plt.style.use('_mpl-gallery')
f1 = plt.figure(figsize=(6, 6))
axes1 = f1.add_axes([0.15, 0.1, 0.8, 0.8])
cmap = ['#d0d1e6', '#99b3d2', '#6993b8', '#40749b', '#1d567a', '#023858']

# Marker spec
markers = ['o', '^', 's', 'd', 'X', '>']

# Load binary results file
filenames = next(os.walk(str(model_source)))[2]  # For later, [1] is dirs and [2] is files
for idx, scenario in enumerate(model_info['orgScenarioNumber']):
    # Check to see if models have been processed already; unzip if no; load if yes
    model_files = analysisUtils.find_files(model_source, f'Scenario_{scenario}_*_{version_no}*', 'files')
        
    file = open(model_source / model_files[3], 'rb')
    freeSurf = pickle.load(file)
    file.close()
    ub = []
    Hs = []
    hMod = 2
    h = 1 * hMod
    HsMod = 1
    TwMod = 1
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
            kh = analysisUtils.qkhf(f, (h / 36))
            cutoff = np.argmax(f >= 2)
            if TwMod != 1:
                f = f
                Pxx_eta = (Pxx_eta*HsMod*1.1)/TwMod
            else:
                Pxx_eta = Pxx_eta*HsMod*1.1
            Hs.append(36 * 4 * np.sqrt(np.sum(Pxx_eta[1:cutoff] * np.mean(np.diff(f)))))
            # Estimate near-bottom orbital velocity with LWT
            # Lowe et al., 2007; Henderson et al., 2017
            c = np.sqrt(9.81 * np.tanh(kh[1:cutoff])) / np.sqrt(kh[1:cutoff] / 0.028)
            n = (1 / 2) * (1 + ((2 * kh[1:cutoff]) / (np.sinh(2 * kh[1:cutoff]))))
            cg = c * n  # Wave celerity 
            F = 36 * 1025 * 9.81 * Pxx_eta[1:cutoff] * cg  # Wave energy spectrum
            urms = np.sqrt(2 * np.trapz(F) * np.mean(np.diff(f[1:cutoff]))) # Lowe eq 10
            ub.append(((8 / np.pi)**(1/2))*urms) 
    
    if idx == 0:
        # Plot Henderson model
        Lambda_not = np.linspace(1e-4, 1e2, 2000)
        Chi = ((((1 + 4 * Lambda_not**2)**0.5) - 1)**(3 / 2)) / ((2**(3 / 2)) * Lambda_not**2)
        axes1.loglog(Lambda_not, Chi, color='k',
                      lw=1)
        
    ubr = np.max(ub)   
    a = (0.0111 / (model_info.deltaS[idx]**2)) / 36
    kappa = 0.4  # Von Karman's constant
    Pi = 0.2  # Cole's wake strength
    z_not = 0.2 / h  # hc/h
    Cd = (kappa**2) * (np.log(h / z_not) + (Pi - 1))**-2
    Lambda_not = (Cd * a * ubr * (model_info.wavePeriod[idx] * 36)) / (4 * np.pi)
    Chi = ((((1 + 4 * Lambda_not**2)**0.5) - 1)**(3 / 2)) / ((2**(3 / 2)) * Lambda_not**2)   
    axes1.loglog(Lambda_not, Chi, color=cmap[idx], marker=markers[idx], markersize=7, mec='k', lw=0)   
    
    # Print out values
    print(f'Scenario: {model_info.orgScenarioNumber[idx]}')
    print(f'Water Depth: {h:.1f}')
    print(f'Wave Height: {model_info.waveHeight[idx]*HsMod*36:.2f}')
    print(f'Wave Period: {model_info.wavePeriod[idx]*TwMod*6}')
    print(f'Spacing: {model_info.spacing[idx]}')
    print(f'Lambda_not: {Lambda_not:.4f}')
    print(f'Chi: {Chi:.4f}\n')
    
axes1.set_ylabel(r'$\chi$')
axes1.set_xlabel(r'$\Lambda_0$')
axes1.set_xlim(1e-4, 10)
axes1.set_ylim(1e-4, 8e-1)
axes1.grid(visible=None)