# -*- coding: utf-8 -*-
"""
Plot Hs/Hs0 vs x (colored by Hs/deltaS) for Paper2_optimizingRestoration

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
twenty = np.where(model_info.wavePeriod == 20)[0].tolist()

# Begin figure
plt.style.use('_mpl-gallery')
f1 = plt.figure(figsize=(8, 12))
axes1 = f1.add_axes([0.15, 0.08, 0.8, 0.27])
axes2 = f1.add_axes([0.15, 0.39, 0.8, 0.27])
axes3 = f1.add_axes([0.15, 0.7, 0.8, 0.27])
ax = [axes1, axes2, axes3]

# Colormaps from ColorBrewer
BuGn = ['#ccece6', '#99d8c9', '#66c2a4', '#41ae76', '#238b45']
OrRd = ['#fdd49e', '#fdbb84', '#fc8d59', '#ef6548', '#d7301f']
PuBu = ['#d0d1e6', '#a6bddb', '#74a9cf', '#3690c0', '#0570b0']
cmap = [BuGn, OrRd, PuBu]

# Marker spec
markers = ['o', '^', 's', 'd', 'P']

# Plot 30 s waves on axes1
wave_gauges = np.linspace(0, 2.6, 10) / 2.6
scenario = [str(x) for x in model_info.orgScenarioNumber[five].tolist()]
waveHeight = model_info.waveHeight[five]
wave_unq = np.unique(waveHeight)
spacing = model_info.deltaS[five]
spce_unq = np.unique(spacing)

size = np.arange(6, 11, 1)
xs = np.linspace(-0.025, 0.025, 5)

for i in range(0, len(wave_unq)):
    for j in range(0, len(spce_unq)):
        idx1 = np.where(waveHeight == wave_unq[i])[0].tolist()
        idx2 = np.where(spacing == spce_unq[j])[0].tolist()
        idx = np.intersect1d(idx1, idx2)[0]
        
        wave = data['wef'][scenario[idx]].loc['Hs'][2:].tolist()
        HsNorm = ([x / wave[0] for x in wave])
    
        x = wave_gauges
        y = HsNorm
        
        x_scaled = spce_unq[j] * 36
        leg = f'$\Delta S = {x_scaled:.1f}$'
        ax[2].plot(x[1:] + xs[j], y[1:], color=cmap[i][j],
                   marker=markers[j],
                   markersize=size[j],
                   label=leg)
    
# Plot 60 s waves on axes2
wave_gauges = np.linspace(0, 5.2, 10) / 5.2
scenario = [str(x) for x in model_info.orgScenarioNumber[ten].tolist()]
waveHeight = model_info.waveHeight[ten]
wave_unq = np.unique(waveHeight)
spacing = model_info.deltaS[ten]
spce_unq = np.unique(spacing)

size = np.arange(6, 11, 1)
xs = np.linspace(-0.025, 0.025, 5)

for i in range(0, len(wave_unq)):
    for j in range(0, len(spce_unq)):
        idx1 = np.where(waveHeight == wave_unq[i])[0].tolist()
        idx2 = np.where(spacing == spce_unq[j])[0].tolist()
        idx = np.intersect1d(idx1, idx2)[0]
        
        wave = data['wef'][scenario[idx]].loc['Hs'][2:].tolist()
        HsNorm = ([x / wave[0] for x in wave])
    
        x = wave_gauges
        y = HsNorm
        
        x_scaled = spce_unq[j] * 36
        leg = f'$\Delta S = {x_scaled:.1f}$'
        ax[1].plot(x[1:] + xs[j], y[1:], color=cmap[i][j],
                   marker=markers[j],
                   markersize=size[j],
                   label=leg)

# Plot 120 s waves on axes3
# wave_gauges = np.linspace(0, 5.2, 10) / 5.2
# scenario = [str(x) for x in model_info.orgScenarioNumber[twenty].tolist()]
# waveHeight = model_info.waveHeight[twenty]
# wave_unq = np.unique(waveHeight)
# spacing = model_info.deltaS[twenty]
# spce_unq = np.unique(spacing)

# size = np.arange(6, 11, 1)
# xs = np.linspace(-0.025, 0.025, 5)

# for i in range(0, len(wave_unq)):
#     for j in range(0, len(spce_unq)):
#         idx1 = np.where(waveHeight == wave_unq[i])[0].tolist()
#         idx2 = np.where(spacing == spce_unq[j])[0].tolist()
#         idx = np.intersect1d(idx1, idx2)[0]
        
#         wave = data['wef'][scenario[idx]].loc['Hs'][2:].tolist()
#         HsNorm = ([x / wave[0] for x in wave])
    
#         x = wave_gauges
#         y = HsNorm
        
#         x_scaled = spce_unq[j] * 36
#         leg = f'$\Delta S = {x_scaled:.1f}$'
#         ax[0].plot(x[1:] + xs[j], y[1:], color=cmap[i][j],
#                    marker=markers[j],
#                    markersize=size[j],
#                    label=leg)

# Multiple legend titles
handles, labels = ax[2].get_legend_handles_labels()
leg1 = ax[0].legend(handles[:5], labels[:5], loc=(0.01, 0),
                    frameon=False,
                    title=r'$H_s = \mathrm{0.05 \ m}$')
title = leg1.get_title()
title.set_size(12)
title.set_weight("bold")

leg2 = ax[0].legend(handles[5:10], labels[5:10], loc=(0.18, 0),
                    frameon=False,
                    title=r'$H_s = \mathrm{0.15 \ m}$')
title = leg2.get_title()
title.set_size(12)
title.set_weight("bold")

leg3 = ax[0].legend(handles[10:15], labels[10:15], loc=(0.35, 0),
                    frameon=False,
                    title=r'$H_s = \mathrm{0.30 \ m}$')
title = leg3.get_title()
title.set_size(12)
title.set_weight("bold")

ax[0].add_artist(leg1)
ax[0].add_artist(leg2)

# Add legend bounding box
ax[0].add_patch(matplotlib.patches.Rectangle((0.01, 0.455), 0.535, 0.24,
                edgecolor='black',
                facecolor='white',
                lw=1))

# Plot Adjustments
ax[0].set_xlabel(r'$x/x_0 \mathrm{(m)}$')
[ax[i].set_ylabel(r'$H_s/H_{s_0}$') for i in range(0, 3)]
ax[0].set_title(r'$T_w = \mathrm{120 \ s \ Models}$')
ax[1].set_title(r'$T_w = \mathrm{60 \ s \ Models}$')
ax[2].set_title(r'$T_w = \mathrm{30 \ s \ Models}$')
ax[1].xaxis.set_ticklabels([])
ax[2].xaxis.set_ticklabels([])
[ax[i].set_xlim(0, 1.05) for i in range(0, 3)]
[ax[i].set_ylim(0.45, 1) for i in range(0, 3)]

# Save figure
if save_figures:
    fname = 'HsNorm_vs_x_Tw_deltaS.png'
    plt.savefig(save_fig_dir / fname, dpi=300, format=None, metadata=None,
                bbox_inches=None, pad_inches=0.1,
                facecolor='auto', edgecolor='auto',
                backend=None)
