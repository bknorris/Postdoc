# -*- coding: utf-8 -*-
"""
Plot Hs/Hs0 vs x (colored by Hs/deltaS) for Paper2_optimizingRestoration

BKN - USGS 2022
"""

import pickle
import os
import sys
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from pathlib import Path
sys.path.insert(1, 'c:\\Users\\bknorris\\Documents\\Scripts\\Postdoc\\Paper2_OptimizingRestoration\\DataAnalysis\\loadModelFiles\\')
from analysisUtils import find_files 
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
data = dict()  # Initialize Dictionary
for idx, scenario in enumerate(model_info['orgScenarioNumber']):
    # Check to see if models have been processed already; unzip if no; load if yes
    model_files = find_files(model_source, f'Scenario_{scenario}_*_{version_no}*', 'files')
    
    file = open(model_source / model_files[0], 'rb')
    print(f'Reading {model_files[0]}')
    avg_fields = pickle.load(file)
    file.close()
        
    file = open(model_source / model_files[3], 'rb')
    print(f'Reading {model_files[3]}')
    freeSurf = pickle.load(file)
    file.close()
    
    model_pprocess = processModels(model_source, [], freeSurf, model_info, idx)
    wef, Ab, fer, ubr = model_pprocess.computeWaveEnergy()
    
    data[str(scenario)] = wef
    # output_file = f'Scenario_{scenario}_WEF_' + version_no + '.dat'
    # file_obj = open(model_source / output_file, mode='wb')
    # pickle.dump(wef, file_obj)
    # file_obj.close()
    
    # eps_norm = model_pprocess.feddersenDissipation(avg_fields['eps'], wef['epsj'][:-1])
    # output_file = f'Scenario_{scenario}_EPS-NORM_' + version_no + '.dat'
    # file_obj = open(model_source / output_file, mode='wb')
    # pickle.dump(eps_norm, file_obj)
    # file_obj.close()

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
        
            wave = data[scenario[idx]].loc['Hs'][2:].tolist()
            HsNorm = ([x / 0.028 for x in wave])
        
            x = wave_gauges
            y = HsNorm
            
            x_scaled = (spce_unq[k] * 36) / 0.4
            leg = f'${x_scaled:.0f}D$'
            ax[i][j].semilogy(x[1:] + xs[j], y[1:], color=cmap[j][k],
                          marker=markers[k],
                          markersize=5,
                          label=leg)


# Plot Adjustments:
# Axis scaling
for i in range(0, 3):
    [ax[i][j].set_xlim(0, 1.03) for j in range(0, 3)]
    [ax[i][j].set_ylim(0.001, 0.08) for j in range(0, 3)]

# Labeling
# ax[0][0].xaxis.set_ticklabels([])
# ax[0][1].xaxis.set_ticklabels([])
# ax[0][1].yaxis.set_ticklabels([])
# ax[0][2].xaxis.set_ticklabels([])
# ax[0][2].yaxis.set_ticklabels([])
# ax[1][0].xaxis.set_ticklabels([])
# ax[1][1].xaxis.set_ticklabels([])
# ax[1][1].yaxis.set_ticklabels([])
# ax[1][2].xaxis.set_ticklabels([])
# ax[1][2].yaxis.set_ticklabels([])
# ax[2][1].yaxis.set_ticklabels([])
# ax[2][2].yaxis.set_ticklabels([])

ax[2][0].set_xlabel(r'$x/x_0$')
ax[2][1].set_xlabel(r'$x/x_0$')
ax[2][2].set_xlabel(r'$x/x_0$')
ax[0][0].set_ylabel(r'$H_s/H_{s_0}$')
ax[1][0].set_ylabel(r'$H_s/H_{s_0}$')
ax[2][0].set_ylabel(r'$H_s/H_{s_0}$')
ax[0][1].set_title(r'$H_s = \mathrm{0.05 \ m \ Models}$')
ax[1][1].set_title(r'$H_s = \mathrm{0.15 \ m \ Models}$')
ax[2][1].set_title(r'$H_s = \mathrm{0.30 \ m \ Models}$')

# Multiple legend titles
handles, labels = ax[2][0].get_legend_handles_labels()
leg1 = ax[2][0].legend(handles, labels, bbox_to_anchor=(1.45, -0.2),
                       frameon=False,
                       title=r'$T_w = \mathrm{30 \ s}$')
title = leg1.get_title()
title.set_size(12)
title.set_weight("bold")

handles, labels = ax[2][1].get_legend_handles_labels()
leg2 = ax[2][1].legend(handles, labels, bbox_to_anchor=(0.75, -0.2),
                       frameon=False,
                       title=r'$T_w = \mathrm{60 \ s}$')
title = leg2.get_title()
title.set_size(12)
title.set_weight("bold")

handles, labels = ax[2][2].get_legend_handles_labels()
leg3 = ax[2][2].legend(handles, labels, bbox_to_anchor=(0.08, -0.2),
                       frameon=False,
                       title=r'$T_w = \mathrm{120 \ s}$')
title = leg3.get_title()
title.set_size(12)
title.set_weight("bold")

# Global font adjustments
SMALL_SIZE = 8
MEDIUM_SIZE = 10
LARGE_SIZE = 12

plt.rc('axes', titlesize=LARGE_SIZE)     # fontsize of the axes title
plt.rc('axes', labelsize=MEDIUM_SIZE)    # fontsize of the x and y labels
plt.rc('xtick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('ytick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('legend', fontsize=SMALL_SIZE)    # legend fontsize

# Save figure
if save_figures:
    fname = 'HsNorm_vs_x_Tw_D_finalV1.png'
    plt.savefig(save_fig_dir / fname, dpi=300, format=None, metadata=None,
                bbox_inches=None, pad_inches=0.1,
                facecolor='auto', edgecolor='auto',
                backend=None)
