# -*- coding: utf-8 -*-
"""
Test for free surface data processing!
"""

import pickle
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy import signal
from pathlib import Path


plt.close('all')
model_source = Path('c:/Users/bknorris/Documents/Models/Paper2_OptimizingRestoration/ModelRuns/Scenarios/postProcessed')
model_file = 'Scenario_27_freeSurf_V1.dat'
file = open(model_source / model_file, 'rb')
fs = pickle.load(file)
file.close()
wave_gauges = np.linspace(0, 5.2, 10)


fig1 = plt.figure(figsize=(8,6))
axes1 = fig1.add_axes([0.1, 0.1, 0.8, 0.8])
axes1.plot(fs['WG4'])
axes1.plot(fs['WG5'])
axes1.plot(fs['WG6'])
# axes1.plot(fs['TimeStep'],fs['WG4'])
# axes1.plot(fs['TimeStep'],fs['WG5'])
# axes1.plot(fs['TimeStep'],fs['WG6'])

axes1.set_xlabel(r'$Timestep \mathrm{(s)}$')
axes1.set_ylabel(r'$\eta \mathrm{(m)}$')

fig2 = plt.figure(figsize=(6,6))
axes2 = fig2.add_axes([0.1, 0.1, 0.8, 0.8])
# f0, Pxx_eta0 = signal.welch(fs['WG1'], 1 / 0.083, nperseg=round(len(fs['WG1']) / 4, 0))
# f1, Pxx_eta1 = signal.welch(fs['WG4'], 1 / 0.083, nperseg=round(len(fs['WG4']) / 4, 0))
# f2, Pxx_eta2 = signal.welch(fs['WG5'], 1 / 0.083, nperseg=round(len(fs['WG5']) / 4, 0))
# f3, Pxx_eta3 = signal.welch(fs['WG6'], 1 / 0.083, nperseg=round(len(fs['WG6']) / 4, 0))
f0, Pxx_eta0 = signal.welch(fs['WG1'][1200:], 1 / 0.083, nperseg=round(len(fs['WG1'][1200:]) / 8, 0))
f1, Pxx_eta1 = signal.welch(fs['WG3'][1200:], 1 / 0.083, nperseg=round(len(fs['WG3'][1200:]) / 8, 0))
f2, Pxx_eta2 = signal.welch(fs['WG4'][1200:], 1 / 0.083, nperseg=round(len(fs['WG4'][1200:]) / 8, 0))
f3, Pxx_eta3 = signal.welch(fs['WG5'][1200:], 1 / 0.083, nperseg=round(len(fs['WG5'][1200:]) / 8, 0))
cutoff = np.argmax(f1 >= 2)

Hs0 = (4 * np.sqrt(np.sum(Pxx_eta0[1:cutoff] * np.mean(np.diff(f0)))))
Hs1 = (4 * np.sqrt(np.sum(Pxx_eta1[1:cutoff] * np.mean(np.diff(f1)))))
Hs2 = (4 * np.sqrt(np.sum(Pxx_eta2[1:cutoff] * np.mean(np.diff(f2)))))
Hs3 = (4 * np.sqrt(np.sum(Pxx_eta3[1:cutoff] * np.mean(np.diff(f3)))))

print(f'WG4 = {Hs1/Hs0}')
print(f'WG5 = {Hs2/Hs0}')
print(f'WG6 = {Hs3/Hs0}')

axes2.plot(f1, Pxx_eta1)
axes2.plot(f2, Pxx_eta2)
axes2.plot(f3, Pxx_eta3)

axes2.set_xlabel(r'$Frequency \mathrm{(Hz)}$')
axes2.set_ylabel(r'$S_{\eta} \mathrm{(m)}$') 