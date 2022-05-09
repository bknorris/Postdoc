# -*- coding: utf-8 -*-
"""
Eps_norm vs waveLength/D (colored by Hs) for Paper2_optimizingRestoration

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
data_file = 'modelPostProcessing_Tw_5_10_20s_V2.dat'
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
twenty = np.where(model_info.wavePeriod == 20)[0].tolist()
waves = [five, ten, twenty]

# Begin figure
plt.style.use('_mpl-gallery')
f1 = plt.figure(figsize=(12, 5))
axes1 = f1.add_axes([0.1, 0.1, 0.28, 0.8])
axes2 = f1.add_axes([0.4, 0.1, 0.28, 0.8])
axes3 = f1.add_axes([0.7, 0.1, 0.28, 0.8])
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
    wavePeriod = model_info.wavePeriod[waves[i]].tolist()
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
            # eta_mean = np.mean(data['wef'][scenario[idx]][gauges][:-1], axis=1)
            # freq = data['wef'][scenario[idx]]['frequency'][:-1]
            # Tp = 1 / freq[np.argmax(eta_mean)]
            Tp = wavePeriod[idx]
            waveLength = (9.81 * (Tp**2)) / (2 * np.pi)
            
            # print(scenario[idx])
            # print(Tp)
            # print(waveLength)
            x_scaled = spce_unq[k] / 0.4
            x = waveLength / x_scaled
            y = np.max(data['eps_norm'][scenario[idx]])
            
            leg = f'${x_scaled:.1f}D$'
            ax[i].loglog(x, y, color=cmap[j][k],
                         marker=markers[k],
                         markersize=8,
                         label=leg)
            
# Plot Adjustments:
# Axis scaling
ax[0].set_xlim(10, 800)
ax[0].set_ylim(2e-4, 1e-2)
ax[1].set_xlim(10, 800)
ax[1].set_ylim(2e-4, 1e-2)
ax[2].set_xlim(10, 800)
ax[2].set_ylim(2e-4, 1e-2)

# Labeling
ax[1].yaxis.set_ticklabels([])
ax[2].yaxis.set_ticklabels([])

ax[0].set_ylabel(r'$\epsilon \left/ (g^3h)^{1/2} \right.$')
ax[0].set_xlabel(r'$\lambda \left/ \Delta S \right.$')
ax[1].set_xlabel(r'$\lambda \left/ \Delta S \right.$')
ax[2].set_xlabel(r'$\lambda \left/ \Delta S \right.$')
ax[0].set_title(r'$T_w = \mathrm{30 \ s \ Models}$')
ax[1].set_title(r'$T_w = \mathrm{60 \ s \ Models}$')
ax[2].set_title(r'$T_w = \mathrm{120 \ s \ Models}$')

# Legends
handles, labels = ax[2].get_legend_handles_labels()
leg1 = ax[0].legend(handles[:5], labels[:5], loc=(0.005, 0),
                    frameon=False,
                    title=r'$H_s = \mathrm{0.05 \ m}$')
title = leg1.get_title()
title.set_size(12)
title.set_weight("bold")

leg2 = ax[0].legend(handles[5:10], labels[5:10], loc=(0.335, 0),
                    frameon=False,
                    title=r'$H_s = \mathrm{0.15 \ m}$')
title = leg2.get_title()
title.set_size(12)
title.set_weight("bold")

leg3 = ax[0].legend(handles[10:15], labels[10:15], loc=(0.67, 0),
                    frameon=False,
                    title=r'$H_s = \mathrm{0.30 \ m}$')
title = leg3.get_title()
title.set_size(12)
title.set_weight("bold")

ax[0].add_artist(leg1)
ax[0].add_artist(leg2)

# Add legend bounding box
ax[0].add_patch(matplotlib.patches.Rectangle((10.5, 2.02e-4), 780, 5.5e-4,
                edgecolor='black',
                facecolor='white',
                lw=1))
            
# Save figure
if save_figures:
    fname = 'Eps_norm_vs_deltaS.png'
    plt.savefig(save_fig_dir / fname, dpi=300, format=None, metadata=None,
                bbox_inches=None, pad_inches=0.1,
                facecolor='auto', edgecolor='auto',
                backend=None)

# try plotting all wave heights in a single plot.