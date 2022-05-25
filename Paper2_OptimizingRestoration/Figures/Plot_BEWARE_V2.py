# -*- coding: utf-8 -*-
"""
Plot Hs/Tp/Reef Flat Width from BEWARE database
for Paper2_optimizingRestoration

BKN - USGS 2022
"""

import netCDF4 as nc
import numpy as np
import pandas as pd
from pathlib import Path
import matplotlib.pyplot as plt
import matplotlib

# Plot commands -- run this at the top of plotting routines
plt.close('all')
matplotlib.rcParams['mathtext.fontset'] = 'custom'
matplotlib.rcParams['mathtext.it'] = 'Arial:italic'
matplotlib.rcParams['mathtext.rm'] = 'Arial'
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42

# Load BEWARE dataset
file_dir = Path('c:/Users/bknorris/Documents/Data/')
file_name = 'BEWARE_Database.nc'
beware = nc.Dataset(file_dir / file_name)

# Load datasets
csv_source = Path('c:/Users/bknorris/Documents/Models/Paper2_OptimizingRestoration/ModelRuns/Scenarios/')
save_fig_dir = Path('c:/Users/bknorris/Documents/Models/Paper2_OptimizingRestoration/Figures')
csv_file = 'BEWARE_OpenFOAM_Model_Data.csv'
save_figures = True

# Load CSV
csv_file = csv_source / csv_file
models = pd.read_csv(csv_file.open())


# Find nearest function
def find_thresh(array, value, thresh):
    low = value - (value * thresh)
    high = value + (value * thresh)
    idx = np.where((array >= low) & (array <= high))
    idx[0].tolist()
    return idx


def find_nearest(array, value):
    array = np.asarray(array)
    idx = np.nanargmin(np.abs(array - value))
    return idx, array[idx]


# BEWARE datasets of interest
bw_Hs_IG = beware['Hm0_IG'][:].tolist()
bw_H0L0 = beware['H0L0'][:].tolist()
bw_Tm = beware['Tm1_0'][:].tolist()
bw_wf = beware['W_reef'][:].tolist()
bw_eta0 = beware['eta0'][:].tolist()
bw_H0 = beware['H0'][:].tolist()

# Estimate deep-water wave period from wavelength
L0 = [bw_H0[x] / bw_H0L0[x] for x in range(0, len(bw_H0))]
bw_T0 = []
for i in range(0, len(bw_H0)):
    if bw_eta0[i] <= 0:
        continue
    if bw_eta0[i] < (L0[i] / 20):
        c = (np.sqrt(9.81 * bw_eta0[i]))
    else:
        c = (np.sqrt((9.81 * L0[i]) / (2 * np.pi)))
    bw_T0.append(L0[i] / c)
    
# Loop through OpenFOAM models and extract values
Wf = []
H0 = []
T0 = []
Hsig = []
Tw = []
for i in range(0, len(models.waterDepth)):
    # Extract Tp and Hs from WG10 in model domain
    Tp = models.wavePeriod[i]
    Hs = models.waveHeight[i]

    # Threshold Hs and Tp values, find nearest Hs value
    idx1 = find_thresh(bw_Hs_IG, Hs, 0.5)
    idx2 = find_thresh(bw_Tm, Tp, 0.5)
    idx3 = np.intersect1d(idx1, idx2)
    idx, Hs0 = find_nearest(np.array(bw_Hs_IG)[idx3.astype(int)], Hs)
    
    Wf.append(bw_wf[idx])
    T0.append(bw_T0[idx])
    H0.append(bw_H0[idx])
    Hsig.append(Hs0)
    Tw.append(Tp)
    
plt.style.use('_mpl-gallery')
f1 = plt.figure(figsize=(6, 5))
ax = f1.add_axes([0.1, 0.1, 0.7, 0.8])
cm = plt.cm.get_cmap('BuPu')
cmaplist = [cm(i) for i in range(cm.N)]
cmap = matplotlib.colors.LinearSegmentedColormap.from_list(
    'Custom cmap', cmaplist, cm.N)
bounds = [150, 200, 250, 300, 350, 400, 450, 500]
norm = matplotlib.colors.BoundaryNorm(bounds, cmap.N)
scatt = ax.scatter(T0, H0, c=Wf, s=60,
                   cmap=cmap,
                   norm=norm,
                   alpha=0.65,
                   linewidths=0.5,
                   edgecolors='k')

ax.set_xlabel(r'$T_0 \ \mathrm{(s)}$')
ax.set_ylabel(r'$H_0 \ \mathrm{(m)}$')
ax.set_title('BEWARE Model Results')

# Colorbar axes
ax2 = f1.add_axes([0.85, 0.1, 0.03, 0.8])
cb = matplotlib.colorbar.ColorbarBase(ax2, cmap=cmap, norm=norm,
                                      spacing='proportional', ticks=bounds, boundaries=bounds, format='%1i')

cb.set_label(r'$W_f \ \mathrm{(m)}$')

# Global font adjustments
SMALL_SIZE = 8
MEDIUM_SIZE = 10
LARGE_SIZE = 12

# Save figure
if save_figures:
    fname = 'BEWARE_values_V2.pdf'
    plt.savefig(save_fig_dir / fname, dpi=300, format='pdf', metadata=None,
                bbox_inches=None, pad_inches=0.1,
                facecolor='auto', edgecolor='auto',
                backend=None)