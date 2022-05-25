# -*- coding: utf-8 -*-
"""
Eps vs waveLength/D (colored by Tw) for Paper2_optimizingRestoration

BKN - USGS 2022
"""

import pickle
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from pathlib import Path
from scipy import stats

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
data_file = 'modelPostProcessing_D2-32_V3.dat'
save_figures = True

# Load binary results file
file = open(model_source / data_file, 'rb')
data = pickle.load(file)
file.close()

# Load CSV
csv_file = csv_source / csv_file
model_info = pd.read_csv(csv_file.open())

# Get data indices from CSV -- filter by waveHeight
five = np.where(model_info.waveHeight == 0.001388889)[0].tolist()
ten = np.where(model_info.waveHeight == 0.004166667)[0].tolist()
twenty = np.where(model_info.waveHeight == 0.008333333)[0].tolist()
waves = [five, ten, twenty]

# Begin figure
plt.style.use('_mpl-gallery')
f1 = plt.figure(figsize=(10, 6))
axes1 = f1.add_axes([0.1, 0.35, 0.28, 0.55])
axes2 = f1.add_axes([0.4, 0.35, 0.28, 0.55])
axes3 = f1.add_axes([0.7, 0.35, 0.28, 0.55])
ax = [axes1, axes2, axes3]

# Colormaps from ColorBrewer
BuGn = ['#ccece6', '#93ccb1', '#63aa83', '#3c885b', '#1b6539', '#00441b']
OrRd = ['#fdd49e', '#f5a677', '#e27b55', '#c65336', '#a42e1b', '#7f0000']
PuBu = ['#d0d1e6', '#99b3d2', '#6993b8', '#40749b', '#1d567a', '#023858']
cmap = [BuGn, OrRd, PuBu]

# Marker spec
markers = ['o', '^', 's', 'd', 'X', '>']

# Plot routine
for i in range(0, len(waves)):
    scenario = model_info.orgScenarioNumber[waves[i]].tolist()
    scenario = [str(x) for x in scenario]
    wavePeriod = model_info.wavePeriod[waves[i]].tolist()
    waveHeight = model_info.waveHeight[waves[i]]
    wave_unq = np.unique(wavePeriod)
    spacing = model_info.deltaS[waves[i]] * 36
    spce_unq = np.unique(spacing)

    for j in range(0, len(wave_unq)):
        for k in range(0, len(spce_unq)):
            idx1 = np.where(wavePeriod == wave_unq[j])[0].tolist()
            idx2 = np.where(spacing == spce_unq[k])[0].tolist()
            idx = np.intersect1d(idx1, idx2)[0]
            Tp = wavePeriod[idx]
            waveLength = (9.81 * (Tp**2)) / (2 * np.pi)

            x_scaled = spce_unq[k] / 0.4
            x = waveLength / x_scaled
            
            eps = np.stack(data['avg_fields'][scenario[idx]]['eps']) * 36
            eps_med = [np.median(eps[x]) for x in range(len(eps))]
            eps_med[0] = np.min(eps[0])

            y = eps_med
            z = np.linspace(0.027, 0.001, len(y)) * 36
            yBins = stats.binned_statistic(z, y, 'mean', bins=22)
            y_mean = np.mean(yBins.statistic[-4:])
            
            leg = f'${x_scaled:.1f}D$'
            ax[i].loglog(x, y_mean, color=cmap[j][k],
                         lw=0,
                         marker=markers[k],
                         markersize=7,
                         mec='k',
                         label=leg)
            
# Plot Adjustments:
# Axis scaling
ax[0].set_xlim(1, 500)
ax[0].set_ylim(2e-6, 5e-3)
ax[1].set_xlim(1, 500)
ax[1].set_ylim(2e-6, 5e-3)
ax[2].set_xlim(1, 500)
ax[2].set_ylim(2e-6, 5e-3)

# Labeling
ax[1].yaxis.set_ticklabels([])
ax[2].yaxis.set_ticklabels([])

ax[0].set_ylabel(r'$\langle \langle \overline{\epsilon} \rangle \rangle \quad \mathrm{(W\left/m^2\right.)}$')
ax[0].set_xlabel(r'$\lambda \left/ D \right.$')
ax[1].set_xlabel(r'$\lambda \left/ D \right.$')
ax[2].set_xlabel(r'$\lambda \left/ D \right.$')
ax[0].set_title(r'$H_s = \mathrm{0.05 \ m \ Models}$')
ax[1].set_title(r'$H_s = \mathrm{0.15 \ m \ Models}$')
ax[2].set_title(r'$H_s = \mathrm{0.30 \ m \ Models}$')

# Multiple legend titles
handles, labels = ax[2].get_legend_handles_labels()
leg1 = ax[0].legend(handles[0:6], labels[0:6], bbox_to_anchor=(1.35, -0.14),
                        frameon=False,
                        title=r'$T_w = \mathrm{30 \ s}$')
title = leg1.get_title()
title.set_size(12)
title.set_weight("bold")

handles, labels = ax[1].get_legend_handles_labels()
leg2 = ax[1].legend(handles[6:12], labels[6:12], bbox_to_anchor=(0.66, -0.14),
                        frameon=False,
                        title=r'$T_w = \mathrm{60 \ s}$')
title = leg2.get_title()
title.set_size(12)
title.set_weight("bold")

handles, labels = ax[2].get_legend_handles_labels()
leg3 = ax[2].legend(handles[12:18], labels[12:18], bbox_to_anchor=(0, -0.14),
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
    fname = 'Eps_vs_lambdaD_V1.pdf'
    plt.savefig(save_fig_dir / fname, dpi=300, format='pdf', metadata=None,
                bbox_inches=None, pad_inches=0.1,
                facecolor='auto', edgecolor='auto',
                backend=None)