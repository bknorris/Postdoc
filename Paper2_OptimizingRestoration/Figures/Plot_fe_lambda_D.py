# -*- coding: utf-8 -*-
"""
Plot Henderson Damping parameter for each model scenario
for Paper2_optimizingRestoration

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


def find_nearest(array, value):
    array = np.asarray(array)
    idx = np.nanargmin(np.abs(array - value))
    return idx, array[idx]

def find_thresh(array, value, thresh):
    low = value - (value * thresh)
    high = value + (value * thresh)
    idx = np.where((array >= low) & (array <= high))
    idx[0].tolist()
    return idx

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
    
    file = open(model_source / model_files[0], 'rb')
    avgFields = pickle.load(file)
    file.close()
        
    file = open(model_source / model_files[3], 'rb')
    freeSurf = pickle.load(file)
    file.close()
    F = []
    ub = []
    ubj = []
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

            cutoff = np.argmax(f >= 2)
            kh = analysisUtils.qkhf(f, 0.028)
            c = np.sqrt(9.81 * np.tanh(kh[1:cutoff])) / np.sqrt(kh[1:cutoff] / 0.028)
            n = (1 / 2) * (1 + ((2 * kh[1:cutoff]) / (np.sinh(2 * kh[1:cutoff]))))
            cg = n * c  # Wave celerity
            F.append(1025 * 9.81 * Pxx_eta[1:cutoff] * cg)  # Wave energy
            
            # Estimate near-bottom orbital velocity with linear wave theory
            aj = np.sqrt(2 * Pxx_eta[1:cutoff])
            omegaj = 2 * np.pi * f[1:cutoff]
            ubj.append((aj * omegaj) / np.sinh(kh[1:cutoff]))
    
    xBin = avgFields['xBins'][0]
    eps = np.stack(avgFields['eps'])
    if idx == 0:
        x_dict = dict()
        x_dict[str(scenario)] = xBin
        eps_dict = dict()
        eps_dict[str(scenario)] = eps
        ubj_dict = dict()
        ubj_dict[str(scenario)] = ubj
    else:
        x_dict[str(scenario)] = xBin
        eps_dict[str(scenario)] = eps
        ubj_dict[str(scenario)] = ubj
        
        
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
ax = [[axes7, axes8, axes9],
      [axes4, axes5, axes6],
      [axes1, axes2, axes3]]

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
    spacing = model_info.deltaS[waves[i]] * 36
    spce_unq = np.unique(spacing)
    
    for j in range(0, len(wave_unq)):
        for k in range(0, len(spce_unq)):
            idx1 = np.where(wavePeriod == wave_unq[j])[0].tolist()
            idx2 = np.where(spacing == spce_unq[k])[0].tolist()
            idx = np.intersect1d(idx1, idx2)[0]
            
            # Reconstruct wave gauges
            if wave_unq[j] == 5:
                x_bnd = (0, 2.6)  # lower/upper
            elif wave_unq[j] == 10:
                x_bnd = (0, 5.2)  # lower/upper
            else:
                x_bnd = (0, 10.4)  # lower/upper
            wave_gauges = np.linspace(x_bnd[0], x_bnd[1], 10)  # 10 evenly spaced WGs
            xs = x_dict[scenario[idx]]
            
            eps_norm = []
            for m in range(0, len(wave_gauges)):
                bin_id = find_thresh(xs, wave_gauges[m], 0.25)
                eps = np.mean(eps_dict[scenario[idx]][:, bin_id], axis=1)
                ubj = ubj_dict[scenario[idx]][m]
                ubr = np.sqrt(sum(ubj**2))
                # eps_norm.append(6 * 4 * eps / (1025 * ubr * sum(ubj**2)))
                eps_norm.append(np.mean(eps))
            
            y = eps_norm
            x = np.linspace(0, 1, 10)
            
            x_scaled = (spce_unq[k] * 36) / 0.4
            leg = f'${x_scaled:.0f}D$'
            ax[i][j].semilogy(x[1:] + xs[j], y[1:], color=cmap[j][k],
                          marker=markers[k],
                          markersize=5,
                          label=leg)

# Plot Adjustments:
# Axis scaling
for i in range(0, 3):
    ax[i].set_xlim(1, 500)
    ax[i].set_ylim(1e-5, 5e-2)

# # Labeling
# ax[1].yaxis.set_ticklabels([])
# ax[2].yaxis.set_ticklabels([])

# ax[0].set_ylabel(r'$F$')
# ax[0].set_xlabel(r'$\lambda \left/ D \right.$')
# ax[1].set_xlabel(r'$\lambda \left/ D \right.$')
# ax[2].set_xlabel(r'$\lambda \left/ D \right.$')
# ax[0].set_title(r'$H_s = \mathrm{0.05 \ m \ Models}$')
# ax[1].set_title(r'$H_s = \mathrm{0.15 \ m \ Models}$')
# ax[2].set_title(r'$H_s = \mathrm{0.30 \ m \ Models}$')

# # Multiple legend titles
# handles, labels = ax[2].get_legend_handles_labels()
# leg1 = ax[0].legend(handles[0:6], labels[0:6], bbox_to_anchor=(1.35, -0.14),
#                        frameon=False,
#                        title=r'$T_w = \mathrm{30 \ s}$')
# title = leg1.get_title()
# title.set_size(12)
# title.set_weight("bold")

# handles, labels = ax[1].get_legend_handles_labels()
# leg2 = ax[1].legend(handles[6:12], labels[6:12], bbox_to_anchor=(0.66, -0.14),
#                        frameon=False,
#                        title=r'$T_w = \mathrm{60 \ s}$')
# title = leg2.get_title()
# title.set_size(12)
# title.set_weight("bold")

# handles, labels = ax[2].get_legend_handles_labels()
# leg3 = ax[2].legend(handles[12:18], labels[12:18], bbox_to_anchor=(0, -0.14),
#                        frameon=False,
#                        title=r'$T_w = \mathrm{120 \ s}$')
# title = leg3.get_title()
# title.set_size(12)
# title.set_weight("bold")

# # Global font adjustments
# SMALL_SIZE = 8
# MEDIUM_SIZE = 10
# LARGE_SIZE = 12
            
# # Save figure
# if save_figures:
#     fname = 'fe_vs_lambdaD_V1.png'
#     plt.savefig(save_fig_dir / fname, dpi=300, format=None, metadata=None,
#                 bbox_inches=None, pad_inches=0.1,
#                 facecolor='auto', edgecolor='auto',
#                 backend=None)
    