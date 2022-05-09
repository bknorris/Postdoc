# -*- coding: utf-8 -*-
"""
Plot Hs/Tp/Reef Flat Width from BEWARE database
for Paper2_optimizingRestoration

BKN - USGS 2022
"""

import pickle
import netCDF4 as nc
import numpy as np
from pathlib import Path
import matplotlib.pyplot as plt
import matplotlib as mpl

plt.close('all')

# Load BEWARE dataset
file_dir = Path('c:/Users/bknorris/Documents/Data/')
file_name = 'BEWARE_Database.nc'
beware = nc.Dataset(file_dir / file_name)

# Load OpenFOAM model dataset
model_source = Path('c:/Users/bknorris/Documents/Models/Paper2_OptimizingRestoration/ModelRuns/Scenarios/postProcessed')
save_fig_dir = Path('c:/Users/bknorris/Documents/Models/Paper2_OptimizingRestoration/Figures')
data_file = 'modelPostProcessing_Final_V1.dat'
save_figures = False

file = open(model_source / data_file, 'rb')
openfoam = pickle.load(file)
file.close()


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
bw_Tm = beware['Tm1_0'][:].tolist()
bw_wf = beware['W_reef'][:].tolist()
bw_eta0 = beware['eta0'][:].tolist()
bw_H0 = beware['H0'][:].tolist()

# Loop through OpenFOAM models and extract values
Wf = []
Tm = []
Hs_IG = []
H0 = []
scenarios = openfoam['wef'].keys()
for num in scenarios:
    # Extract Tp and Hs from WG10 in model domain
    f = openfoam['wef'][num]['frequency'][:-1]
    Tp = 1 / f[np.argmax(openfoam['wef'][num]['WG10'][:-1])] * 6
    Hs = openfoam['wef'][num].loc['Hs'][-1] * 36

    # Threshold Hs and Tp values, find nearest Hs value
    idx1 = find_thresh(bw_Hs_IG, Hs, 0.1)
    idx2 = find_thresh(bw_Tm, Tp, 0.2)
    idx3 = np.intersect1d(idx1, idx2)
    
    idx, Hs0 = find_nearest(np.array(bw_Hs_IG)[idx3.astype(int)], Hs)
    # idx, Hs0 = find_nearest(bw_Tm, Tp)
    
    print('Scenario ' + num)
    # print(Hs)
    # print(Hs0)
    # print(bw_H0[idx])
    print(bw_eta0[idx])
    # print(bw_wf[idx])
    # print('\n')
    Wf.append(bw_wf[idx])
    Tm.append(bw_Tm[idx])
    Hs_IG.append(Hs0)
    H0.append(bw_H0[idx])
    
    
plt.style.use('_mpl-gallery')
f1 = plt.figure(figsize=(6, 5))
ax = f1.add_axes([0.1, 0.1, 0.7, 0.8])
cm = plt.cm.get_cmap('BuPu')
cmaplist = [cm(i) for i in range(cm.N)]
cmap = mpl.colors.LinearSegmentedColormap.from_list(
    'Custom cmap', cmaplist, cm.N)
bounds = [0, 25, 50, 75, 100, 125, 150, 175, 200, 225, 250, 275, 300, 325, 350, 375, 400]
norm = mpl.colors.BoundaryNorm(bounds, cmap.N)
Hs_scaled = ([x * 60 for x in H0])
scatt = ax.scatter(Tm, Hs_IG, c=Wf, s=Hs_scaled,
                   cmap=cmap,
                   norm=norm,
                   alpha=0.65,
                   linewidths=0.5,
                   edgecolors='k')

ax.set_xlabel(r'$T_m \ \mathrm{(s)}$')
ax.set_ylabel(r'$H_{m0,IG} \ \mathrm{(m)}$')

# Colorbar axes
ax2 = f1.add_axes([0.85, 0.1, 0.03, 0.8])
cb = mpl.colorbar.ColorbarBase(ax2, cmap=cmap, norm=norm,
                               spacing='proportional', ticks=bounds, boundaries=bounds, format='%1i')
cb.set_label(r'$W_f \ \mathrm{(m)}$')

# Size legend
inset = f1.add_axes([0.65, 0.15, 0.13, 0.16])
s_min = np.argmin(Hs_scaled)
s_max = np.argmax(Hs_scaled)
s_med = np.argsort(Hs_scaled)[len(Hs_scaled)//2]
xs = [1, 1, 1]
ys = [0.1, 0.235, 0.4]
size = np.unique(Hs_scaled)
inset.scatter(xs, ys, c='w', s=size, linewidths=0.5, edgecolors='k')
labeltext = ['1.0', '2.0', '3.0']
for i in range (0,3):
    inset.text(xs[i]+0.02, ys[i]-0.035, labeltext[i])

inset.text(0.99, 0.53, r'$H_0 \ \mathrm{(m)}$', fontsize=11)
inset.set_xlim(0.98, 1.05)
inset.set_ylim(0.02, 0.7)
inset.xaxis.set_ticklabels([])
inset.yaxis.set_ticklabels([])
inset.grid(False)

# Global font adjustments
SMALL_SIZE = 8
MEDIUM_SIZE = 10
LARGE_SIZE = 12

plt.rc('axes', titlesize=LARGE_SIZE)     # fontsize of the axes title
plt.rc('axes', labelsize=LARGE_SIZE)    # fontsize of the x and y labels
plt.rc('xtick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('ytick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('legend', fontsize=SMALL_SIZE)    # legend fontsize

# Save figure
if save_figures:
    fname = 'BEWARE_Tm_Hm0_H0_Wf.png'
    plt.savefig(save_fig_dir / fname, dpi=300, format=None, metadata=None,
                bbox_inches=None, pad_inches=0.1,
                facecolor='auto', edgecolor='auto',
                backend=None)