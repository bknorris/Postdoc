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
        Processing utilities for OpenFOAM model output files. Note this script
        is designed to work with the Paper2_OptimizingRestoration model files only!
        
        Inputs:
            source_path: path of input file. Must be the base model folder (typically ".../ModelRuns/Scenarios/")
            model_files: names of BIN files created by loadModelFiles.py in one list, and
                        the index position of these files in the model_info CSV.
            model_info: Data Frame containing CSV values from a data processing script.
            
        Outputs:
            wef, Ab, fer, ubr = computeWaveEnergy()
                wef: Data Frame containing wave energy flux values for 10 wave gauges;
                Ab: Wave orbital excursion at the bed [field scale, m]
                fer: Wave friction factor (fe) [Dimensionless]
                ubr: Near-bottom orbital velocity [field scale, m/s]
        '''
        
        # Load files
        print('Loading model files...')
        fid = open(source_path / model_files[0][0], mode='rb')
        fields = pickle.load(fid)
        fid = open(source_path / model_files[0][1], mode='rb')
        free_surf = pickle.load(fid)
        
        # Basic wave characteristics
        Tp = model_info['wavePeriod'][model_files[1]]
        Hs = model_info['waveHeight'][model_files[1]]
        h = model_info['waterDepth'][model_files[1]]
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
        Hs_spec = []
        np.seterr(divide='ignore', invalid='ignore')
        for gauge, data in self.free_surf.iteritems():
            if gauge != 'TimeStep':
                fs = 1 / 0.083  # sample frequency from models
                f, Pxx_eta = signal.welch(data, fs, nperseg=round(len(data) / 4, 0))
                kh = analysisUtils.qkhf(f, self.h)
                c = np.sqrt(9.81 * np.tanh(kh) / np.sqrt(kh / self.h))
                n = (1 / 2) * (1 + ((2 * kh) / (np.sinh(2 * kh))))
                cg = c * n  # Wave celerity
                cutoff = np.argmax(f >= 2)
                Hs_spec.append(4 * np.sqrt(np.sum(Pxx_eta[1:cutoff] * np.mean(np.diff(f)))))
                F.append(1025 * 9.81 * Pxx_eta[1:cutoff] * cg[1:cutoff])  # Wave energy
                
                # Estimate near-bottom orbital velocity with LWT
                omegaj.append(2 * np.pi * f[1:cutoff])
                ubj.append((36 * 10 * 2 * Pxx_eta[1:cutoff] * 2 * np.pi * f[1:cutoff]) / np.sinh(kh[1:cutoff]))
                # Note: 36 * 10 is the physical scaling * fudge factor for viscosity diff between
                # the model and field scales.
                
        # Mean wave energy dissipation
        del_x = self.wave_gauges[-1] - self.wave_gauges[0]
        ubr = np.sqrt(np.sum(ubj[-1])**2)
        omegar = np.sum(omegaj[-1] * (ubj[-1]**2)) / np.sum(ubj[-1]**2)
        Ab = ubr / omegar
        epsj = -1 * (F[-1] - F[0]) / del_x
        fej = (4 * epsj) / (1025 * ubr * ubj[-1]**2)
        fer = np.sum(fej * ubj[-1]**2) / np.sum(ubj[-1]**2)
        
        # Create data frame and return results
        wgs = [f'WG{x}' for x in list(range(1, len(self.wave_gauges) + 1))]
        wef = pd.DataFrame(dict(zip(wgs, F)))
        wef.insert(loc=0, column='epsj', value=epsj)
        wef.insert(loc=0, column='frequency', value=f[1:cutoff])
        wef.loc['Hs'] = [np.nan, np.nan] + Hs_spec
        
        return wef, Ab, fer, ubr
    
    def avgTKE(self):
        print('Averaging turbulence results...')
        
        