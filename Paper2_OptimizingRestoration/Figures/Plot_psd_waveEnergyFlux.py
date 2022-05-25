# -*- coding: utf-8 -*-
"""
Power spectra of F at WG2 and WG9 (10% and 90% of domain)
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
from scipy.stats import chi2
from pathlib import Path
sys.path.insert(1, 'c:\\Users\\bknorris\\Documents\\Scripts\\Postdoc\\Paper2_OptimizingRestoration\\DataAnalysis\\loadModelFiles\\')
import analysisUtils 

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
csv_file = 'modelPostProcessing_PowerSpectraF.csv'
version_no = 'V3'

save_figures = True

# Load CSV
csv_file = csv_source / csv_file
model_info = pd.read_csv(csv_file.open())
spacing = model_info.spacing.tolist()

# Begin figure
plt.style.use('_mpl-gallery')
f1 = plt.figure(figsize=(8, 8))
axes1 = f1.add_axes([0.1, 0.58, 0.4, 0.3])
axes2 = f1.add_axes([0.1, 0.22, 0.4, 0.3])
axes3 = f1.add_axes([0.58, 0.58, 0.4, 0.3])
axes4 = f1.add_axes([0.58, 0.22, 0.4, 0.3])
ax = [axes1, axes2, axes3, axes4]

# Colormaps from ColorBrewer
BuGn = ['#93ccb1', '#1b6539']
PuBu = ['#99b3d2', '#1d567a']
cmap = [BuGn, BuGn, PuBu, PuBu]

# Load binary results file
filenames = next(os.walk(str(model_source)))[2]  # For later, [1] is dirs and [2] is files
data = dict()  # Initialize Dictionary
for idx, scenario in enumerate(model_info['orgScenarioNumber']):
    # Check to see if models have been processed already; unzip if no; load if yes
    model_files = analysisUtils.find_files(model_source, f'Scenario_{scenario}_*_{version_no}*', 'files')
        
    file = open(model_source / model_files[3], 'rb')
    print(f'Reading {model_files[3]}')
    freeSurf = pickle.load(file)
    file.close()
    
    F = []
    F_lo = []  # 95% CI lower bound
    F_hi = []  # 95% CI upper bound
    for gauge, data in freeSurf.iteritems():
        if gauge != 'TimeStep':
            fs = 1 / 0.083  # sample frequency from models
            data = signal.detrend(data)
            seg = round(len(data) / 6, 0)
            noverlap = seg / 2
            f, Pxx_eta = signal.welch(data, fs,
                                      window='hanning',
                                      nperseg=seg,
                                      nfft=seg * 6,
                                      detrend=False,
                                      scaling='density',
                                      noverlap=noverlap)
            kh = analysisUtils.qkhf(f, 0.028)
            c = np.sqrt(9.81 * np.tanh(kh)) / np.sqrt(kh / 0.028)
            n = (1 / 2) * (1 + ((2 * kh) / (np.sinh(2 * kh))))
            cg = c * n  # Wave celerity
            cutoff = np.argmax(f >= 2)
            DOF = 2 * np.round(len(data) / seg, 0)
            alfa = 1 - 0.95
            v = 2 * DOF
            chi = chi2.ppf([1 - alfa / 2, alfa / 2], v)
            chi = v / chi
            
            F_lo.append(36 * 1025 * 9.81 * Pxx_eta[1:cutoff] * chi[0] * cg[1:cutoff])
            F.append(36 * 1025 * 9.81 * Pxx_eta[1:cutoff] * cg[1:cutoff])  # Wave energy
            F_hi.append(36 * 1025 * 9.81 * Pxx_eta[1:cutoff] * chi[1] * cg[1:cutoff])
     
    leg = f'${spacing[idx]:.1f}D$'
    ax[idx].plot(f[1:cutoff] / 6, F_lo[1], color=cmap[idx][0], ls=':', lw=1, label=leg)
    ax[idx].plot(f[1:cutoff] / 6, F[1], color=cmap[idx][0], ls='-', lw=2, label=leg)
    ax[idx].plot(f[1:cutoff] / 6, F_hi[1], color=cmap[idx][0], ls=':', lw=1, label=leg)
    
    ax[idx].plot(f[1:cutoff] / 6, F_lo[8], color=cmap[idx][1], ls=':', lw=1, label=leg)
    ax[idx].plot(f[1:cutoff] / 6, F[8], color=cmap[idx][1], ls='-', lw=2, label=leg)
    ax[idx].plot(f[1:cutoff] / 6, F_hi[8], color=cmap[idx][1], ls=':', lw=1, label=leg)
    
    ax[idx].plot(f[1:cutoff] / 6, F[1] - F[8], color='k', ls='-', lw=1.5, label=leg)
           
# Plot Adjustments:
# Axis scaling
[ax[i].set_xlim(0, 0.33) for i in range(0, 2)]
[ax[i].set_ylim(0, 0.22) for i in range(0, 2)]
[ax[i].set_xlim(0, 0.16) for i in range(2, 4)]
[ax[i].set_ylim(0, 0.16) for i in range(2, 4)]
[ax[i].grid(False) for i in range(0, 4)]
#[ax[i].invert_xaxis() for i in range(0, 4)]


# Labeling
ax[0].xaxis.set_ticklabels([])
ax[2].xaxis.set_ticklabels([])

ax[0].set_ylabel(r'$F \quad \mathrm{(kg \cdot m \cdot s^{-3} \left/ Hz \right.)}$')
ax[1].set_ylabel(r'$F \quad \mathrm{(kg \cdot m \cdot s^{-3} \left/ Hz \right.)}$')
ax[1].set_xlabel(r'$T_w \quad \mathrm{(s)}$')
ax[3].set_xlabel(r'$T_w \quad \mathrm{(s)}$')

ax[0].set_title(r'$H_s = \mathrm{0.15 \ m}, T_w = \mathrm{30 \ s}, 2D$')
ax[1].set_title(r'$H_s = \mathrm{0.15 \ m}, T_w = \mathrm{30 \ s}, 32D$')
ax[2].set_title(r'$H_s = \mathrm{0.15 \ m}, T_w = \mathrm{120 \ s}, 2D$')
ax[3].set_title(r'$H_s = \mathrm{0.15 \ m}, T_w = \mathrm{120 \ s}, 32D$')


# Multiple legend titles
handles, s = ax[0].get_legend_handles_labels()
labels = [r'$0.1 \cdot W_f$', '95% CI', r'$0.9 \cdot W_f$', '95% CI', r'$\Delta F$']
h = [handles[1], handles[0], handles[4], handles[3], handles[6]]
leg1 = ax[0].legend(h, labels, frameon=True)

handles, s = ax[1].get_legend_handles_labels()
h = [handles[1], handles[0], handles[4], handles[3], handles[6]]
leg1 = ax[1].legend(h, labels, frameon=True)

handles, s = ax[2].get_legend_handles_labels()
h = [handles[1], handles[0], handles[4], handles[3], handles[6]]
leg1 = ax[2].legend(h, labels, frameon=True)

handles, s = ax[3].get_legend_handles_labels()
h = [handles[1], handles[0], handles[4], handles[3], handles[6]]
leg1 = ax[3].legend(h, labels, frameon=True)

# Global font adjustments
SMALL_SIZE = 8
MEDIUM_SIZE = 10
LARGE_SIZE = 12
            
# Save figure
if save_figures:
    fname = 'F_PSD_Tw_30_120_Hs_0_15_V2.pdf'
    plt.savefig(save_fig_dir / fname, dpi=300, format='pdf', metadata=None,
                bbox_inches=None, pad_inches=0.1,
                facecolor='auto', edgecolor='auto',
                backend=None)