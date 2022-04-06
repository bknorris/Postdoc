# -*- coding: utf-8 -*-
"""
Created on Mon Apr  4 12:47:43 2022

@author: bknorris
"""

import pickle
import matplotlib.pyplot as plt
from mycolorpy import colorlist as mcp
import numpy as np
import pandas as pd
from pathlib import Path

# Define file paths:
csv_source = Path('c:/Users/bknorris/Documents/Models/Paper2_OptimizingRestoration/ModelRuns/Scenarios/')
model_source = Path('c:/Users/bknorris/Documents/Models/Paper2_OptimizingRestoration/ModelRuns/Scenarios/postProcessed')
csv_file = 'modelPostProcessing.csv'
data_file = 'modelPostProcessing_5s_10s_V1.dat'

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
# Hs/Hs0 vs x (colored by Tp)
plt.style.use('_mpl-gallery')
f1 = plt.figure(figsize=(6.5, 6.5))
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

cmap = mcp.gen_color(cmap="Set1", n=3)
markers = ['o', 's', '^', 'd', 'p']
for i in range(0, len(wave_unq)):
    for j in range(0, len(spce_unq)):
        ys = np.intersect1d(y[waveHeight == wave_unq[i]], y[spacing == spce_unq[j]])
        scaled = spce_unq[j]*36
        leg = f'$\Delta S: {scaled:.3f}$'
        plt.plot(x, ys, color=cmap[i], marker=markers[j], label=leg)

# axes.invert_xaxis()
axes.set_xlabel(r'$x/x_0 \mathrm{(m)}$')
axes.set_ylabel(r'$H_s/H_{s_0}$')    
plt.legend()

  
# delF/delx vs x (colored by Tp)
# Eps_norm(z) vs Z (colored by TP)
# max(Eps_norm) vs Tp
# Runup at shore vs Tp