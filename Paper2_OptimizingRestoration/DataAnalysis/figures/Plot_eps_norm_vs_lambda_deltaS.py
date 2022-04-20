# -*- coding: utf-8 -*-
"""
Eps_norm vs waveLength/DeltaS (colored by Hs) for Paper2_optimizingRestoration

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
waves = [five, ten]

# Begin figure
plt.style.use('_mpl-gallery')
f1 = plt.figure(figsize=(12, 5))
axes1 = f1.add_axes([0.1, 0.1, 0.27, 0.8])
axes2 = f1.add_axes([0.4, 0.1, 0.27, 0.8])
axes3 = f1.add_axes([0.7, 0.1, 0.27, 0.8])
ax = [axes1, axes2, axes3]

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
    # wavePeriod = model_info.wavePeriod[waves[i]].tolist()
    waveHeight = model_info.waveHeight[waves[i]]
    wave_unq = np.unique(waveHeight)
    spacing = model_info.deltaS[waves[i]] * 36
    spce_unq = np.unique(spacing)

    for j in range(0, len(wave_unq)):
        for k in range(0, len(spce_unq)):
            idx1 = np.where(waveHeight == wave_unq[j])[0].tolist()
            idx2 = np.where(spacing == spce_unq[k])[0].tolist()
            idx = np.intersect1d(idx1, idx2)[0]
            # Tp = 1 / wavePeriod[idx]
            gauges = data['wef'][scenario[idx]].keys()[2:]
            eta_mean = np.mean(data['wef'][scenario[idx]][gauges][:-1], axis=1)
            freq = data['wef'][scenario[idx]]['frequency'][:-1]
            Tp = 1 / freq[np.argmax(eta_mean)]
            waveLength = (9.81 * (Tp**2)) / (2 * np.pi)
            
            print(scenario[idx])
            print(Tp)
            print(waveLength)
            
            y = np.max(data['eps_norm'][scenario[idx]])
            x = waveLength / spce_unq[k]
            
            leg = f'$\Delta S = {spce_unq[j]:.1f}$'
            ax[i].loglog(x, y, color=cmap[j][k],
                       marker=markers[k],
                       markersize=6,
                       label=leg)
            
            
# x-axis is wavelength / deltaS [-]
# y-axis is epsilon/(g^3h)^0.5 [-]
# try plotting all wave heights in a single plot, and plot different Tw on different plots. 