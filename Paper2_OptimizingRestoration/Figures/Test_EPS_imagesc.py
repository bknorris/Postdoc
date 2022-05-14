# -*- coding: utf-8 -*-
"""
TEST PLOT EPS
"""

from pathlib import Path
import pickle
import os
import sys
import pandas as pd
from plottingUtils import imagesc
import numpy as np
import matplotlib.pyplot as plt
sys.path.insert(1, 'c:\\Users\\bknorris\\Documents\\Scripts\\Postdoc\\Paper2_OptimizingRestoration\\DataAnalysis\\loadModelFiles\\')
from analysisUtils import find_files 

plt.close('all')
# Define file paths:
csv_source = Path('c:/Users/bknorris/Documents/Models/Paper2_OptimizingRestoration/ModelRuns/Scenarios/')
model_source = Path('c:/Users/bknorris/Documents/Models/Paper2_OptimizingRestoration/ModelRuns/Scenarios/postProcessed')
save_fig_dir = Path('c:/Users/bknorris/Documents/Models/Paper2_OptimizingRestoration/Figures')
csv_file = 'modelPostProcessing_TEST3.csv'
version_no = 'V2'
save_figures = False


def extents(f):
    delta = f[1] - f[0]
    return [f[0] - delta / 2, f[-1] + delta / 2]


# Load CSV
csv_file = csv_source / csv_file
model_info = pd.read_csv(csv_file.open())

# Create figure
plt.style.use('_mpl-gallery')
f1 = plt.figure(figsize=(8, 10.5))
axes1 = f1.add_axes([0.1, 0.03, 0.8, 0.15])
axes2 = f1.add_axes([0.1, 0.19, 0.8, 0.15])
axes3 = f1.add_axes([0.1, 0.35, 0.8, 0.15])
axes4 = f1.add_axes([0.1, 0.51, 0.8, 0.15])
axes5 = f1.add_axes([0.1, 0.67, 0.8, 0.15])
axes6 = f1.add_axes([0.1, 0.83, 0.8, 0.15])

ax = [axes1, axes2, axes3,
      axes4, axes5, axes6]
caxis = [1e-7, 5e-5]
c_map = 'viridis'

# Load binary results file
filenames = next(os.walk(str(model_source)))[2]  # For later, [1] is dirs and [2] is files
data = dict()  # Initialize Dictionary
for idx, scenario in enumerate(model_info['orgScenarioNumber']):
    # Check to see if models have been processed already; unzip if no; load if yes
    model_files = find_files(model_source, f'Scenario_{scenario}_*_{version_no}*', 'files')
    
    file = open(model_source / model_files[0], 'rb')
    print(f'Reading {model_files[0]}')
    avg_fields = pickle.load(file)
    file.close()
    
    eps = np.stack(avg_fields['eps'])
    xBins = avg_fields['xBins'][0]
    zBins = avg_fields['zBins'][0]
    
    ax[idx].imshow(eps, aspect='auto', extent=extents(xBins) + extents(zBins),
                   vmin=caxis[0], vmax=caxis[1], origin='lower', cmap=c_map)

[ax[x].xaxis.set_ticklabels([]) for x in range(1,6)]

fname = 'Eps_Tw_30_Hs_0_15.png'
plt.savefig(save_fig_dir / fname, dpi=300, format=None, metadata=None,
            bbox_inches=None, pad_inches=0.1,
            facecolor='auto', edgecolor='auto',
            backend=None)