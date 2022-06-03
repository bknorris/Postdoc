# -*- coding: utf-8 -*-
"""
Calculate z0 for Cd calculations for Paper2_optimizingRestoration

BKN - USGS 2022
"""

import pickle
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats
from sklearn import linear_model
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
model_source = Path('c:/Users/bknorris/Documents/Models/Paper2_OptimizingRestoration/ModelRuns/Scenarios/postProcessed')
csv_file = 'modelPostProcessing_mod1.csv'
data_file = 'modelPostProcessing_D2-32_V3.dat'

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

# Plot routine
for i in range(0, len(waves)):
    scenario = model_info.orgScenarioNumber[waves[i]].tolist()
    scenario = [str(x) for x in scenario]
    wavePeriod = model_info.wavePeriod[waves[i]].tolist()
    wave_unq = np.unique(wavePeriod)
    spacing = model_info.deltaS[waves[i]]
    spce_unq = np.unique(spacing)

    for j in range(0, len(wave_unq)):
        for k in range(0, len(spce_unq)):
            idx1 = np.where(wavePeriod == wave_unq[j])[0].tolist()
            idx2 = np.where(spacing == spce_unq[k])[0].tolist()
            idx = np.intersect1d(idx1, idx2)[0]
            x_scaled = (spce_unq[k] * 36) / 0.4
            Umag = np.mean(np.stack(data['avg_fields'][scenario[idx]]['Umag']), axis=1)
            zBins = data['avg_fields'][scenario[idx]]['zBins'][0]
            
            y = Umag
            x = np.log(-1 * zBins)
            
            y = y[-19:-8]
            x = x[-19:-8]
            regr = linear_model.LinearRegression()
            regr.fit(x.reshape(-1, 1), y.reshape(-1, 1))
            znot = np.exp(-1 * regr.intercept_[0] / regr.coef_[0][0])
            # print(f'Scenario: {scenario[idx]}')
            # print(f'Spacing: {x_scaled:.2f}')
            print(f'{znot * 36}')
            # print(f'D/z_not: {0.028 / znot}\n')

 
