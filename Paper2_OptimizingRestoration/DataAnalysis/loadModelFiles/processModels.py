# -*- coding: utf-8 -*-
"""
Modules for post processing OpenFOAM freeSurface and field data
for Paper2_OptimizingRestoration models.

BKN - USGS PCMSC 2022
"""

import pandas as pd
import numpy as np
import pickle
import time
import analysisUtils
from scipy import signal


class processModels: 
    def __init__(self, source_path, model_files, model_info):
        '''
        CHANGE THIS HERE!
        
        Processing utilities for OpenFOAM model output files. Note this script
        is designed to work with the Paper2_OptimizingRestoration model files only!
        
        Inputs:
            source_path: path of input file. Must be the base model folder (typically ".../ModelRuns/Scenarios/")
            file_name: model folder name (typically "Scenario_XX")
            wave_period: integer value for the wave period of each model (used to
                         determine which processing step to run)
            
        Outputs:
            pandas array containing model freeSurface and field data saved in BIN format
        '''
        
        # Load files
        print('Loading model files...')
        fid = open(source_path / model_files[0][0], mode='rb')
        fields = pickle.load(fid)
        fid = open(source_path / model_files[0][1], mode='rb')
        free_surf = pickle.load(fid)
        
        # Basic wave characteristics
        Tp = model_info['wavePeriod'][model_files[0]]
        Hs = model_info['waveHeight'][model_files[0]]
        h = model_info['waterDepth'][model_files[0]]
        gamma = Hs / h
        wavelength = (9.81 * (Tp**2)) / (2 * np.pi)
        steepness = Hs / wavelength
        
        # Define wave gauges
        if Tp == 5:
            x_bnd = (0, 2.6)  # lower/upper
        elif Tp == 10:
            x_bnd = (0, 5.2)  # lower/upper
        elif Tp == 20:
            x_bnd = (0, 10.4)  # lower/upper
        else:
            x_bnd = (0, 21.6)  # lower/upper
        wave_gauges = np.linspace(x_bnd[0], x_bnd[1], 10)  # 10 evenly spaced WGs
        
        self.fields = fields
        self.free_surf = free_surf
        self.Tp = Tp
        self.Hs = Hs
        self.h = h
        self.gamma = gamma
        self.wavelength = wavelength
        self.steepness = steepness
        self.wave_gauges = wave_gauges
    
    def computeWaveEnergy(self):
        print('Performing wave energy calculations...')
        F = []
        omegaj = []
        ubj = []
        for gauge, data in free_surf.iteritems():
            if gauge != 'TimeStep':
                fs = 1 / 0.083  # sample frequency from models
                f, Pxx_eta = signal.welch(data, fs, nperseg=round(len(data) / 4, 0))
                kh = analysisUtils.qkhf(f, h)
                c = np.sqrt(9.81 * np.tanh(kh[1:]) / np.sqrt(kh[1:] / h))
                n = (1 / 2) * (1 + ((2 * kh[1:]) / (np.sinh(2 * kh[1:]))))
                cg = c * n  # Wave celerity
                F.append(1025 * 9.81 * Pxx_eta[1:] * cg)  # Wave energy
                
                # Estimate near-bottom orbital velocity with LWT
                aj = np.sqrt(2 * Pxx_eta[1:])
                omegaj.append(2 * np.pi * f[1:])
                ubj.append((2 * Pxx_eta[1:] * 2 * np.pi * f[1:]) / np.sinh(kh[1:]))
                
        # Mean wave energy dissipation
        del_x = wave_gauges[-1] - wave_gauges[0]
        ubr = np.sqrt(np.sum(ubj[-1])**2)
        omegar = np.sum(omegaj[-1] * (ubj[-1]**2)) / np.sum(ubj[-1]**2)
        Ab = ubr / omegar
        epsj = -1 * (F[-1] - F[0]) / del_x
        fej = (4 * epsj) / (1025 * ubr * ubj[-1]**2)
        fer = np.sum(fej * ubj[-1]**2) / np.sum(ubj[-1]**2) ## I don't think this is right...
        
        # SEE: Paper1_plotWaveDissipationFactorVsTp_V2_alt2
        # Average TKE Dissipation in time
        