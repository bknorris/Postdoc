# -*- coding: utf-8 -*-
"""
Plot Model vs Field Scale Umag and EPS for Paper2_optimizingRestoration

BKN - USGS 2022
"""

import pickle
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from sklearn import linear_model
from sklearn.metrics import r2_score
from pathlib import Path

plt.close('all')

# Define file paths:
csv_source = Path('c:/Users/bknorris/Documents/Models/Paper2_OptimizingRestoration/ModelRuns/Scenarios/')
model_source = Path('c:/Users/bknorris/Documents/Models/Paper2_OptimizingRestoration/ModelRuns/Scenarios/postProcessed')
save_fig_dir = Path('c:/Users/bknorris/Documents/Models/Paper2_OptimizingRestoration/Figures')
csv_file = 'modelPostProcessing_fieldScale.csv'
model_file = 'modelPostProcessing_D2-32_V3.dat'
field_file = 'modelPostProcessing_FS_V1.dat'
save_figures = True

# Load binary results file
file = open(model_source / model_file, 'rb')
model = pickle.load(file)
file.close()

file = open(model_source / 'FieldScale' / field_file, 'rb')
field = pickle.load(file)
file.close()

# Load CSV
csv_file = csv_source / csv_file
model_info = pd.read_csv(csv_file.open())

plt.style.use('_mpl-gallery')
f1 = plt.figure(figsize=(10, 6))
axes1 = f1.add_axes([0.1, 0.1, 0.38, 0.8])
axes2 = f1.add_axes([0.58, 0.1, 0.38, 0.8])

# Colormaps from ColorBrewer
BuGn = ['#ccece6', '#93ccb1', '#63aa83', '#3c885b', '#1b6539', '#00441b',
        '#ccece6', '#93ccb1', '#63aa83', '#3c885b', '#1b6539', '#00441b']
PuBu = ['#d0d1e6', '#99b3d2', '#6993b8', '#40749b', '#1d567a', '#023858',
        '#d0d1e6', '#99b3d2', '#6993b8', '#40749b', '#1d567a', '#023858']

# Marker spec
markers = ['o', '^', 's', 'd', 'X', '>',
           'o', '^', 's', 'd', 'X', '>']

U_R2 = []
E_R2 = []
for i in range(0, len(model_info.orgScenarioNumber)):
    scenario = str(model_info.orgScenarioNumber[i])
    waveHeight = model_info.waveHeight[i]
    spacing = model_info.spacing[i]
    U_model = np.mean(np.stack(model['avg_fields'][scenario]['Umag']), axis=1)
    U_field = np.mean(np.stack(field['avg_fields'][scenario]['Umag']), axis=1)
    E_model = np.mean(np.stack(model['avg_fields'][scenario]['eps']), axis=1)
    E_field = np.mean(np.stack(field['avg_fields'][scenario]['eps']), axis=1)
    
    print(f'Scenario: {scenario}')
    print(f'Spacing: {spacing}')
    regr = linear_model.LinearRegression()
    regr.fit(U_model.reshape(-1, 1), U_field.reshape(-1, 1))
    print(f'Umag regression coef: {regr.coef_[0][0]:.2f}')
    
    leg = f'${spacing:.0f}D$'
    axes1.plot(U_model, U_field, color=BuGn[i],
               lw=0,
               marker=markers[i],
               markersize=5,
               label=leg)
    U_R2.append(regr.coef_[0][0])
    
    regr = linear_model.LinearRegression()
    regr.fit(E_model.reshape(-1, 1), E_field.reshape(-1, 1))
    print(f'Eps regression coef: {regr.coef_[0][0]:.2f}\n')
    axes2.loglog(E_model, E_field, color=PuBu[i],
               lw=0,
               marker=markers[i],
               markersize=5,
               label=leg)
    E_R2.append(regr.coef_[0][0])
    