# -*- coding: utf-8 -*-
"""
Eps_norm(z) vs Z (colored by Hs & deltaS) for Paper2_optimizingRestoration

BKN - USGS 2022
"""

import pickle
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from pathlib import Path

plt.close('all')

# Define file paths:
csv_source = Path('c:/Users/bknorris/Documents/Models/Paper2_OptimizingRestoration/ModelRuns/Scenarios/')
model_source = Path('c:/Users/bknorris/Documents/Models/Paper2_OptimizingRestoration/ModelRuns/Scenarios/postProcessed')
save_fig_dir = Path('c:/Users/bknorris/Documents/Models/Paper2_OptimizingRestoration/Figures')
csv_file = 'modelPostProcessing.csv'
data_file = 'modelPostProcessing_5s_10s_V2.dat'
save_figures = True

# Load binary results file
file = open(model_source / data_file, 'rb')
data = pickle.load(file)
file.close()

# Load CSV
csv_file = csv_source / csv_file
model_info = pd.read_csv(csv_file.open())

# Get data indices from CSV -- filter by wavePeriod
five = np.where(model_info.wavePeriod == 5)[0].tolist()
ten = np.where(model_info.wavePeriod == 10)[0].tolist()
#twenty = np.where(model_info.wavePeriod == 20)[0].tolist()
waves = [five, ten]

# Begin figure
plt.style.use('_mpl-gallery')
f1 = plt.figure(figsize=(8, 12))
axes1 = f1.add_axes([0.1, 0.08, 0.25, 0.28])
axes2 = f1.add_axes([0.4, 0.08, 0.25, 0.28])
axes3 = f1.add_axes([0.7, 0.08, 0.25, 0.28])
axes4 = f1.add_axes([0.1, 0.39, 0.25, 0.28])
axes5 = f1.add_axes([0.4, 0.39, 0.25, 0.28])
axes6 = f1.add_axes([0.7, 0.39, 0.25, 0.28])
axes7 = f1.add_axes([0.1, 0.7, 0.25, 0.28])
axes8 = f1.add_axes([0.4, 0.7, 0.25, 0.28])
axes9 = f1.add_axes([0.7, 0.7, 0.25, 0.28])
ax = [[axes7, axes8, axes9],
      [axes4, axes5, axes6],
      [axes1, axes2, axes3]]

# Colormaps from ColorBrewer
BuGn = ['#ccece6', '#99d8c9', '#66c2a4', '#41ae76', '#238b45']
OrRd = ['#fdd49e', '#fdbb84', '#fc8d59', '#ef6548', '#d7301f']
PuBu = ['#d0d1e6', '#a6bddb', '#74a9cf', '#3690c0', '#0570b0']
cmap = [BuGn, OrRd, PuBu]

# Marker spec
markers = ['o', '^', 's', 'd', 'P']
# Plot routine
for i in range(0, len(waves)):
    scenario = model_info.orgScenarioNumber[waves[i]].tolist()
    scenario = [str(x) for x in scenario]
    waveHeight = model_info.waveHeight[waves[i]]
    wave_unq = np.unique(waveHeight)
    spacing = model_info.deltaS[waves[i]]
    spce_unq = np.unique(spacing)

    for j in range(0, len(wave_unq)):
        for k in range(0, len(spce_unq)):
            idx1 = np.where(waveHeight == wave_unq[j])[0].tolist()
            idx2 = np.where(spacing == spce_unq[k])[0].tolist()
            idx = np.intersect1d(idx1, idx2)[0]
            
            x = data['eps_norm'][scenario[idx]]
            y = np.linspace(-0.028, -0.001, len(x)) * 36
            
            x_scaled = spce_unq[j] * 36
            leg = f'$\Delta S = {x_scaled:.1f}$'
            ax[i][j].plot(x, y, color=cmap[j][k],
                        marker=markers[j],
                        markersize=6,
                        label=leg)

        
# Plot Adjustments
# ax[0].set_xlabel(r'$x/x_0 \mathrm{(m)}$')
# [ax[i].set_ylabel(r'$H_s/H_{s_0}$') for i in range(0, 3)]
# ax[0].set_title(r'$T_w = \mathrm{120 \ s \ Models}$')
# ax[1].set_title(r'$T_w = \mathrm{60 \ s \ Models}$')
# ax[2].set_title(r'$T_w = \mathrm{30 \ s \ Models}$')
# ax[1].xaxis.set_ticklabels([])
# ax[2].xaxis.set_ticklabels([])
for i in range(0, 3):
    [ax[i][j].set_xlim(-1e-4, 5e-3) for j in range(0, 3)]
    [ax[i][j].set_ylim(-1.1, -0.6) for j in range(0, 3)]
