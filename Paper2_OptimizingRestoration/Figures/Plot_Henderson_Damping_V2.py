# -*- coding: utf-8 -*-
"""
Plot Henderson Damping parameter for each model scenario
for Paper2_optimizingRestoration

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

h = 1  # water depth (m)
Hs = 1  # wave height multiplier
# Load binary results file
filenames = next(os.walk(str(model_source)))[2]  # For later, [1] is dirs and [2] is files
for idx, scenario in enumerate(model_info['orgScenarioNumber']):
    # Check to see if models have been processed already; unzip if no; load if yes
    model_files = analysisUtils.find_files(model_source, f'Scenario_{scenario}_*_{version_no}*', 'files')
        
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
            kh = analysisUtils.qkhf(f, (h / 36))
            cutoff = np.argmax(f >= 2)
            Pxx_eta = Pxx_eta * Hs
            
            # Estimate near-bottom orbital velocity with LWT
            # Lowe et al., 2007; Henderson et al., 2017
            c = np.sqrt(9.81 * np.tanh(kh[1:cutoff])) / np.sqrt(kh[1:cutoff] / (h / 36))
            n = (1 / 2) * (1 + ((2 * kh[1:cutoff]) / (np.sinh(2 * kh[1:cutoff]))))
            cg = c * n  # Wave celerity 
            F = 36 * 1025 * 9.81 * Pxx_eta[1:cutoff] * cg  # Wave energy spectrum
            urms = np.sqrt(2 * np.trapz(F) * np.mean(np.diff(f[1:cutoff]))) # Lowe eq 10
            ub.append(((8 / np.pi)**(1/2))*urms) 
    
    if idx == 0:
        cal_dict = dict()
        cal_dict[str(scenario)] = ub
        
    else:
        cal_dict[str(scenario)] = ub
        
        
# Get data indices from CSV -- filter by waveHeight
five = np.where(model_info.waveHeight == 0.001388889)[0].tolist()
ten = np.where(model_info.waveHeight == 0.004166667)[0].tolist()
twenty = np.where(model_info.waveHeight == 0.008333333)[0].tolist()
waves = [five, ten, twenty]

plt.style.use('_mpl-gallery')
f1 = plt.figure(figsize=(10, 6))
axes1 = f1.add_axes([0.08, 0.35, 0.28, 0.55])
axes2 = f1.add_axes([0.39, 0.35, 0.28, 0.55])
axes3 = f1.add_axes([0.7, 0.35, 0.28, 0.55])
ax = [axes1, axes2, axes3]

# Colormaps from ColorBrewer
BuGn = ['#ccece6', '#93ccb1', '#63aa83', '#3c885b', '#1b6539', '#00441b']
OrRd = ['#fdd49e', '#f5a677', '#e27b55', '#c65336', '#a42e1b', '#7f0000']
PuBu = ['#d0d1e6', '#99b3d2', '#6993b8', '#40749b', '#1d567a', '#023858']
cmap = [BuGn, OrRd, PuBu]

# Marker spec
markers = ['o', '^', 's', 'd', 'X', '>']


for i in range(0, len(waves)):
    scenario = model_info.orgScenarioNumber[waves[i]].tolist()
    scenario = [str(x) for x in scenario]
    wavePeriod = model_info.wavePeriod[waves[i]].tolist()
    wave_unq = np.unique(wavePeriod)
    spacing = model_info.deltaS[waves[i]]
    spce_unq = np.unique(spacing)
    
    
    # Plot Henderson model
    Lambda_not = np.linspace(1e-4, 1e2, 2000)
    Chi = ((((1 + 4 * Lambda_not**2)**0.5) - 1)**(3 / 2)) / ((2**(3 / 2)) * Lambda_not**2)
    ax[i].loglog(Lambda_not, Chi, color='k',
                 lw=1)
    
    for j in range(0, len(wave_unq)):
        for k in range(0, len(spce_unq)):
            idx1 = np.where(wavePeriod == wave_unq[j])[0].tolist()
            idx2 = np.where(spacing == spce_unq[k])[0].tolist()
            idx = np.intersect1d(idx1, idx2)[0]
            ubr = np.max(cal_dict[scenario[idx]])
            
            a = (0.0111 / (spce_unq[k]**2)) / 36
            # kappa = 0.4  # Von Karman's constant
            # Pi = 0.2  # Cole's wake strength
            # z_not = 0.2
            # Cd = (kappa**2) * (np.log(h / z_not) + (Pi - 1))**-2
            # print(f'Cd: {Cd}')
            Cd = model_info.Cd[idx]
            Lambda_not = (Cd * a * ubr * (wave_unq[j] * 36)) / (4 * np.pi)
            Chi = ((((1 + 4 * Lambda_not**2)**0.5) - 1)**(3 / 2)) / ((2**(3 / 2)) * Lambda_not**2)
            
            # KC number
            # KC = (ubr * wave_unq[j]) / 0.011
            # print(f'KC number: {KC}')
            
            # khc
            # wavelength = np.sqrt(9.81 * (1 / 36)) * wave_unq[j]
            # khc = 2 * np.pi * wavelength * 0.0056
            # print(f'k*hc: {khc}')
            
            x_scaled = (spce_unq[k] * 36) / 0.4
            leg = f'${x_scaled:.0f}D$'
            ax[i].loglog(Lambda_not, Chi, color=cmap[j][k],
                         lw=0,
                         marker=markers[k],
                         markersize=7,
                         mec='k',
                         label=leg)


# Plot Adjustments:
# Axis scaling
for i in range(0, 3):
    ax[i].set_xlim(1e-4, 10)
    ax[i].set_ylim(1e-4, 8e-1)
    ax[i].grid(visible=None)

# Labeling
ax[1].yaxis.set_ticklabels([])
ax[2].yaxis.set_ticklabels([])
ax[0].set_ylabel(r'$\chi$')
ax[0].set_xlabel(r'$\Lambda_0$')
ax[1].set_xlabel(r'$\Lambda_0$')
ax[2].set_xlabel(r'$\Lambda_0$')
ax[0].set_title(r'$H_s = \mathrm{0.05 \ m \ Models}$')
ax[1].set_title(r'$H_s = \mathrm{0.15 \ m \ Models}$')
ax[2].set_title(r'$H_s = \mathrm{0.30 \ m \ Models}$')

# Multiple legend titles
handles, labels = ax[2].get_legend_handles_labels()
leg1 = ax[0].legend(handles[0:6], labels[0:6], bbox_to_anchor=(1.35, -0.13),
                       frameon=False,
                       title=r'$T_w = \mathrm{30 \ s}$')
title = leg1.get_title()
title.set_size(12)
title.set_weight("bold")

handles, labels = ax[1].get_legend_handles_labels()
leg2 = ax[1].legend(handles[6:12], labels[6:12], bbox_to_anchor=(0.66, -0.13),
                       frameon=False,
                       title=r'$T_w = \mathrm{60 \ s}$')
title = leg2.get_title()
title.set_size(12)
title.set_weight("bold")

handles, labels = ax[2].get_legend_handles_labels()
leg3 = ax[2].legend(handles[12:18], labels[12:18], bbox_to_anchor=(0, -0.13),
                       frameon=False,
                       title=r'$T_w = \mathrm{120 \ s}$')
title = leg3.get_title()
title.set_size(12)
title.set_weight("bold")

# Global font adjustments
SMALL_SIZE = 8
MEDIUM_SIZE = 10
LARGE_SIZE = 12
            
# Save figure
if save_figures:
    fname = 'Henderson_Damping_V2.pdf'
    plt.savefig(save_fig_dir / fname, dpi=300, format='pdf', metadata=None,
                bbox_inches=None, pad_inches=0.1,
                facecolor='auto', edgecolor='auto',
                backend=None)
    