# -*- coding: utf-8 -*-
"""
Created on Mon Apr  4 12:47:43 2022

@author: bknorris
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
data_file = 'modelPostProcessing_5s_10s_V1.dat'
save_figures = False

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

# Figures:
# Hs/Hs0 vs x (colored by Hs/deltaS)
plt.style.use('_mpl-gallery')
f1 = plt.figure(figsize=(10, 6.5))
axes = f1.add_axes([0.12, 0.1, 0.8, 0.8])

# Plot 5 s waves first
wave_gauges = np.linspace(0, 2.6, 10)[::-1] / 2.6
scenario = model_info.orgScenarioNumber[five].tolist()
scenario = [str(x) for x in scenario]
Hs = []
for index in scenario:
    wave = data['wef'][index].loc['Hs'][2:].tolist()
    Hs.append([x / wave[0] for x in wave])

x = wave_gauges
y = np.array(Hs)
waveHeight = model_info.waveHeight[five].tolist()
wave_unq = np.unique(waveHeight)
spacing = model_info.deltaS[five].tolist()
spce_unq = np.unique(spacing)

# Colormaps from ColorBrewer
BuGn = ['#ccece6', '#99d8c9', '#66c2a4', '#41ae76', '#238b45']
OrRd = ['#fdd49e', '#fdbb84', '#fc8d59', '#ef6548', '#d7301f']
PuBu = ['#d0d1e6', '#a6bddb', '#74a9cf', '#3690c0', '#0570b0']
cmap = [BuGn, OrRd, PuBu]

# Marker spec
markers = ['o', '^', 's', 'd', 'P']
size = np.arange(6, 11, 1)
xs = np.linspace(-0.025, 0.025, 5)

for i in range(0, len(wave_unq)):
    for j in range(0, len(spce_unq)):
        ys = np.intersect1d(y[waveHeight == wave_unq[i]], y[spacing == spce_unq[j]])
        x_scaled = spce_unq[j] * 36
        leg = f'$\Delta S = {x_scaled:.1f}$'
        plt.plot(x[:-1] + xs[j], ys[:-1], color=cmap[i][j],
                 marker=markers[j],
                 markersize=size[j],
                 label=leg)

# Multiple legend titles
handles, labels = axes.get_legend_handles_labels()
leg1 = axes.legend(handles[:5], labels[:5], loc=(0.01, 0),
                   frameon=False,
                   title=r'$H_s = \mathrm{0.05 \ m}$')
title = leg1.get_title()
title.set_size(12)
title.set_weight("bold")

leg2 = axes.legend(handles[5:10], labels[5:10], loc=(0.18, 0),
                   frameon=False,
                   title=r'$H_s = \mathrm{0.15 \ m}$')
title = leg2.get_title()
title.set_size(12)
title.set_weight("bold")

leg3 = axes.legend(handles[10:15], labels[10:15], loc=(0.35, 0),
                   frameon=False,
                   title=r'$H_s = \mathrm{0.30 \ m}$')
title = leg3.get_title()
title.set_size(12)
title.set_weight("bold")

axes.add_artist(leg1)
axes.add_artist(leg2)

# Add legend bounding box
axes.add_patch(matplotlib.patches.Rectangle((0.01, 0.5525), 0.505, 0.119,
               edgecolor='black',
               facecolor='white',
               lw=1))

# Plot Adjustments
axes.set_xlabel(r'$x/x_0 \mathrm{(m)}$')
axes.set_ylabel(r'$H_s/H_{s_0}$')
axes.set_title(r'$T_w = \mathrm{30 \ s \ Models}$')
plt.xlim(0, 1.05)
plt.ylim(0.55, 1)

# Save figure
if save_figures:
    fname = 'HsNorm_vs_x_Tw_5s.png'
    plt.savefig(save_fig_dir / fname, dpi=300, format=None, metadata=None,
                bbox_inches=None, pad_inches=0.1,
                facecolor='auto', edgecolor='auto',
                backend=None)
  
    
# eps_norm vs deltaS (colored by Hs/deltaS)
plt.style.use('_mpl-gallery')
f2 = plt.figure(figsize=(5, 4))
axes = f2.add_axes([0.175, 0.12, 0.8, 0.8])

# Plot 5 s waves first
scenario = model_info.orgScenarioNumber[five].tolist()
scenario = [str(x) for x in scenario]
WEF = []
for index in scenario:
    wave = data['eps_norm'][index][:-1]
    WEF.append(np.max(wave))
    
x = model_info.deltaS[five] * 36
y = np.array(WEF)
waveHeight = model_info.waveHeight[five].tolist()
wave_unq = np.unique(waveHeight)
spacing = model_info.deltaS[five].tolist()
spce_unq = np.unique(spacing)

# Colormaps from ColorBrewer
BuGn = ['#ccece6', '#99d8c9', '#66c2a4', '#41ae76', '#238b45']
OrRd = ['#fdd49e', '#fdbb84', '#fc8d59', '#ef6548', '#d7301f']
PuBu = ['#d0d1e6', '#a6bddb', '#74a9cf', '#3690c0', '#0570b0']
cmap = [BuGn, OrRd, PuBu]

# Marker spec
markers = ['o', '^', 's', 'd', 'P']

# Plot
for i in range(0, len(wave_unq)):
    for j in range(0, len(spce_unq)):
        ys = np.intersect1d(y[waveHeight == wave_unq[i]], y[spacing == spce_unq[j]])
        xs = np.intersect1d(x[waveHeight == wave_unq[i]], x[spacing == spce_unq[j]])
        leg = f'$\Delta S = {x_scaled:.1f}$'
        plt.semilogy(xs, ys, color=cmap[i][j],
                     marker=markers[j],
                     markersize=8,
                     label=leg)

# Plot Adjustments
axes.set_ylabel(r'$\epsilon \left/ (g^3h)^{1/2} \right.$')
axes.set_xlabel(r'$\Delta S \ \mathrm{(m)}$')
axes.set_title(r'$T_w = \mathrm{30 \ s \ Models}$')
plt.ylim(1e-3, 5e-3)
plt.xlim(0.55, 1.45)
plt.yticks(np.linspace(1e-3, 5e-3, 5))

# Save figure
if save_figures:
    fname = 'Eps_norm_vs_deltaS_Tw_5s.png'
    plt.savefig(save_fig_dir / fname, dpi=300, format=None, metadata=None,
                bbox_inches=None, pad_inches=0.1,
                facecolor='auto', edgecolor='auto',
                backend=None)


# Eps_norm(z) vs Z (colored by deltaS)
plt.style.use('_mpl-gallery')
f3 = plt.figure(figsize=(10, 6.5))
axes1 = f3.add_axes([0.1, 0.1, 0.25, 0.8])
axes2 = f3.add_axes([0.4, 0.1, 0.25, 0.8])
axes3 = f3.add_axes([0.7, 0.1, 0.25, 0.8])
ax = [axes1, axes2, axes3]

# Colormaps from ColorBrewer
BuGn = ['#ccece6', '#99d8c9', '#66c2a4', '#41ae76', '#238b45']
OrRd = ['#fdd49e', '#fdbb84', '#fc8d59', '#ef6548', '#d7301f']
PuBu = ['#d0d1e6', '#a6bddb', '#74a9cf', '#3690c0', '#0570b0']
cmap = [BuGn, OrRd, PuBu]

# Marker spec
markers = ['o', '^', 's', 'd', 'P']

# Plot 5 s waves first
scenario = model_info.orgScenarioNumber[five].tolist()
scenario = [str(x) for x in scenario]
waveHeight = model_info.waveHeight[five]
wave_unq = np.unique(waveHeight)
spacing = model_info.deltaS[five]
spce_unq = np.unique(spacing)

count = 1
for i in range(0, len(wave_unq)):
    for j in range(0, len(spce_unq)):
        idx1 = np.where(waveHeight == wave_unq[i])[0].tolist()
        idx2 = np.where(spacing == spce_unq[j])[0].tolist()
        idx = np.intersect1d(idx1, idx2)[0]
        
        x = data['eps_norm'][scenario[idx]]
        y = np.linspace(-0.028, -0.001, len(x)) * 36
        
        x_scaled = spce_unq[j] * 36
        leg = f'$\Delta S = {x_scaled:.1f}$'
        ax[i].plot(x, y, color=cmap[i][j],
                 marker=markers[j],
                 markersize=6,
                 label=leg)
    
    count += 1
    
# Plot Adjustments
axes1.set_ylim(-1.1, 0)
axes1.set_xlim(-1e-4, 5e-3)
axes1.set_xticks(np.linspace(1e-3, 5e-3, 3))
axes1.legend(loc='upper right', frameon=True, framealpha=1)
axes1.set_title(r'$T_w = \ \mathrm{0.05 m}$')
axes1.set_ylabel('Depth (m)')
axes1.set_xlabel(r'$\epsilon \left/ (g^3h)^{1/2} \right.$')

axes2.set_ylim(-1.1, 0)
axes2.set_xlim(-1e-4, 5e-3)
axes2.set_xticks(np.linspace(1e-3, 5e-3, 3))
axes2.legend(loc='upper right')
axes2.set_title(r'$T_w = \ \mathrm{0.15 m}$')
axes2.yaxis.set_ticklabels([])
axes2.set_xlabel(r'$\epsilon \left/ (g^3h)^{1/2} \right.$')

axes3.set_ylim(-1.1, 0)
axes3.set_xlim(-1e-4, 5e-3)
axes3.set_xticks(np.linspace(1e-3, 5e-3, 3))
axes3.legend(loc='upper right')
axes3.set_title(r'$T_w = \ \mathrm{0.30 m}$')
axes3.yaxis.set_ticklabels([])
axes3.set_xlabel(r'$\epsilon \left/ (g^3h)^{1/2} \right.$')

# Save figure
if save_figures:
    fname = 'Eps_norm_vs_z_5s.png'
    plt.savefig(save_fig_dir / fname, dpi=300, format=None, metadata=None,
                bbox_inches=None, pad_inches=0.1,
                facecolor='auto', edgecolor='auto',
                backend=None)


# max(Eps_norm) vs Tp
plt.style.use('_mpl-gallery')
f4 = plt.figure(figsize=(5, 4))
axes = f4.add_axes([0.175, 0.12, 0.8, 0.8])

# Colormaps from ColorBrewer
BuGn = ['#ccece6', '#99d8c9', '#66c2a4', '#41ae76', '#238b45']
OrRd = ['#fdd49e', '#fdbb84', '#fc8d59', '#ef6548', '#d7301f']
PuBu = ['#d0d1e6', '#a6bddb', '#74a9cf', '#3690c0', '#0570b0']
cmap = [BuGn, OrRd, PuBu]

# Marker spec
markers = ['o', '^', 's', 'd', 'P']

# Plot 5 s waves first
scenario = model_info.orgScenarioNumber[five].tolist()
scenario = [str(x) for x in scenario]
waveHeight = model_info.waveHeight[five]
wave_unq = np.unique(waveHeight)
wavePeriod = model_info.wavePeriod[five]
spacing = model_info.deltaS[five]
spce_unq = np.unique(spacing)

count = 1
for i in range(0, len(wave_unq)):
    for j in range(0, len(spce_unq)):
        idx1 = np.where(waveHeight == wave_unq[i])[0].tolist()
        idx2 = np.where(spacing == spce_unq[j])[0].tolist()
        idx = np.intersect1d(idx1, idx2)[0]
        
        y = np.max(data['eps_norm'][scenario[idx]])
        x = wavePeriod.tolist()[idx]
        
        x_scaled = spce_unq[j] * 36
        leg = f'$\Delta S = {x_scaled:.1f}$'
        plt.plot(x, y, color=cmap[i][j],
                 marker=markers[j],
                 markersize=6,
                 label=leg)
    
    count += 1


# Runup at shore vs Tp
