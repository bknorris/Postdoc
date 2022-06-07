# -*- coding: utf-8 -*-
"""
Plot Henderson Damping parameter for each model scenario
for Paper2_optimizingRestoration

BKN - USGS 2022
"""

import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from pathlib import Path

# Plot commands -- run this at the top of plotting routines
plt.close('all')
matplotlib.rcParams['mathtext.fontset'] = 'custom'
matplotlib.rcParams['mathtext.it'] = 'Arial:italic'
matplotlib.rcParams['mathtext.rm'] = 'Arial'
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42

# Define file paths:
csv_source = Path('c:/Users/bknorris/Documents/Models/Paper2_OptimizingRestoration/ModelRuns/Scenarios/')
save_fig_dir = Path('c:/Users/bknorris/Documents/Models/Paper2_OptimizingRestoration/Figures')
csv_file = 'Henderson_ExtraValues.csv'

save_figures = True

# Load CSV
csv_file = csv_source / csv_file
models = pd.read_csv(csv_file.open())

# Get data indices from CSV -- filter by waveHeight
one = np.where(models.waveHeight == 0.05)[0].tolist()
two = np.where(models.waveHeight == 0.15)[0].tolist()
three = np.where(models.waveHeight == 0.30)[0].tolist()
four = np.where(models.waveHeight == 0.6)[0].tolist()
five = np.where(models.waveHeight == 0.9)[0].tolist()
waves = [one, two, three, four, five]

plt.style.use('_mpl-gallery')
f1 = plt.figure(figsize=(6, 10))
axes1 = f1.add_axes([0.1, 0.76, 0.35, 0.2])
axes2 = f1.add_axes([0.1, 0.53, 0.35, 0.2])
axes3 = f1.add_axes([0.1, 0.30, 0.35, 0.2])
axes4 = f1.add_axes([0.58, 0.65, 0.35, 0.2])
axes5 = f1.add_axes([0.58, 0.42, 0.35, 0.2])
ax = [axes3, axes2, axes1, axes5, axes4]

# Colormaps from ColorBrewer
BuGn = ['#ccece6', '#93ccb1', '#63aa83', '#3c885b', '#1b6539', '#00441b']
OrRd = ['#fdd49e', '#f5a677', '#e27b55', '#c65336', '#a42e1b', '#7f0000']
PuBu = ['#d0d1e6', '#99b3d2', '#6993b8', '#40749b', '#1d567a', '#023858']
Reds = ['#fee0d2', '#fcbba1', '#fc9272', '#fb6a4a', '#ef3b2c', '#cb181d']
Purples = ['#dadaeb', '#bcbddc', '#9e9ac8', '#807dba', '#6a51a3', '#54278f']
cmap1 = [BuGn, OrRd, PuBu]
cmap2 = [PuBu, Reds, Purples]

# Marker spec
markers = ['o', '^', 's', 'd', 'X', '>']


for i in range(0, len(waves)):
    wavePeriod = models.wavePeriod[waves[i]].tolist()
    wave_unq = np.unique(wavePeriod)
    spacing = models.spacing[waves[i]]
    spce_unq = np.unique(spacing)
    if i < 3:
        cmap = cmap1
    else:
        cmap = cmap2
        
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
            Lambda_not = models.lambdaNot[waves[i][idx]]
            Chi = models.chi[waves[i][idx]]
            
            x_scaled = spce_unq[k]
            leg = f'${x_scaled:.0f}D$'
            ax[i].loglog(Lambda_not, Chi, color=cmap[j][k],
                         lw=0,
                         marker=markers[k],
                         markersize=7,
                         alpha=1,
                         mec='k',
                         label=leg)

# Plot Adjustments:
# Axis scaling
for i in range(0, 5):
    ax[i].set_xlim(1e-4, 10)
    ax[i].set_ylim(1e-4, 8e-1)
    ax[i].grid(visible=None)

# Labeling
ax[2].xaxis.set_ticklabels([])
ax[1].xaxis.set_ticklabels([])
ax[4].xaxis.set_ticklabels([])

ax[0].set_ylabel(r'$\chi$')
ax[1].set_ylabel(r'$\chi$')
ax[2].set_ylabel(r'$\chi$')
ax[3].set_ylabel(r'$\chi$')
ax[4].set_ylabel(r'$\chi$')

ax[0].set_xlabel(r'$\Lambda_0$')
ax[3].set_xlabel(r'$\Lambda_0$')
ax[2].set_title(r'$H_{s0} = \mathrm{0.05 \ m \ Models}$')
ax[1].set_title(r'$H_{s0} = \mathrm{0.15 \ m \ Models}$')
ax[0].set_title(r'$H_{s0} = \mathrm{0.30 \ m \ Models}$')
ax[4].set_title(r'$H_{s0} = \mathrm{0.60 \ m \ Calculated}$')
ax[3].set_title(r'$H_{s0} = \mathrm{0.90 \ m \ Calculated}$')


# Multiple legends
# THIS IS DONE IN ILLUSTRATOR!

# Global font adjustments
SMALL_SIZE = 8
MEDIUM_SIZE = 10
LARGE_SIZE = 12
            
# Save figure
if save_figures:
    fname = 'Henderson_DampingExtraValues_V2.pdf'
    plt.savefig(save_fig_dir / fname, dpi=300, format='pdf', metadata=None,
                bbox_inches=None, pad_inches=0.1,
                facecolor='auto', edgecolor='auto',
                backend=None)
    