# -*- coding: utf-8 -*-
"""
Modules for post processing OpenFOAM freeSurface and field data
for Paper2_OptimizingRestoration models.

BKN - USGS PCMSC 2022
"""

import pandas as pd
import numpy as np
import pickle
import analysisUtils
from scipy import signal
from scipy import stats


class processModels:
    def __init__(self, source_path, fields, free_surf, model_info, idx):
        '''
        Processing utilities for OpenFOAM model output files. Note this script
        is designed to work with the Paper2_OptimizingRestoration model files only!
        
        Inputs:
            source_path: path of input file. Must be the base model folder (typically ".../ModelRuns/Scenarios/")
            fields: field data from the RAW OpenFOAM output; "Scenario_#/Model/postProcessing/bathySample.."
            free_surf: free surface data from the RAW OpenFOAM output; "Scenario_#/Model/postProcessing/freeSurface.."
            idx: index corresponding to the model scenario in model_info
            model_info: Data Frame containing CSV values from a data processing script.
            
        Outputs:
            wef, Ab, fer, ubr = computeWaveEnergy():
                wef: Data Frame containing wave energy flux values for 10 wave gauges;
                Ab: Wave orbital excursion at the bed [field scale, m]
                fer: Wave friction factor (fe) [Dimensionless]
                ubr: Near-bottom orbital velocity [field scale, m/s]
                
            avg_fields = avgFields():
                avg_fields: Contains the following x-z binned values:
                    eps: TKE Dissipation Rate [m^2/s^3]
                    Umag: Magnitude of U [field scale, m/s]
                    k: TKE [m^2/s^2]
                    xBins: constant bin size along the model x axis
                    zBins: 26 bins along the model z axis
                    
            eps_norm = feddersenDissipation():
                eps_norm: depth profile of surfzone normalized dissipation
                          (see Eq. 9 of Fedderson, 2012)
            
            Lambda_not = hendersonDamping():
                Lambda_not: Dimensionless damping parameter from
                            Henderson et al. 2017
        '''
        
        # Basic wave characteristics
        Tp = model_info['wavePeriod'][idx]
        Hs = model_info['waveHeight'][idx]
        h = model_info['waterDepth'][idx]
        gamma = Hs / h
        wavelength = (9.81 * (Tp**2)) / (2 * np.pi)
        steepness = Hs / wavelength
        
        # Canopy characteristics
        d = model_info['d'][idx]
        hc = model_info['hc'][idx]
        deltaS = model_info['deltaS'][idx]
        
        # Define wave gauges
        if Tp == 5:
            x_bnd = (0, 2.6)  # lower/upper
        elif Tp == 10:
            x_bnd = (0, 5.2)  # lower/upper
        elif Tp == 20:
            x_bnd = (0, 10.4)  # lower/upper
        else:
            x_bnd = (0, 21.6)  # lower/upper
        # x_bnd = (0, 94)
        wave_gauges = np.linspace(x_bnd[0], x_bnd[1], 10)  # 10 evenly spaced WGs
        
        self.fields = fields
        self.free_surf = free_surf
        self.Tp = Tp
        self.Hs = Hs
        self.h = h
        self.d = d
        self.hc = hc
        self.deltaS = deltaS
        self.gamma = gamma
        self.wavelength = wavelength
        self.steepness = steepness
        self.wave_gauges = wave_gauges
        self.fs = 0.083  # model sample frequency
        
    def computeWaveEnergy(self):
        print('Performing wave energy calculations...')
        F = []
        omegaj = []
        ubj = []
        Hs_spec = []
        np.seterr(divide='ignore', invalid='ignore')
        for gauge, data in self.free_surf.iteritems():
            if gauge != 'TimeStep':
                fs = 1 / self.fs  # sample frequency from models
                data = signal.detrend(data)
                seg = round(len(data) / 2, 0)
                noverlap = seg / 2
                f, Pxx_eta = signal.welch(data, fs,
                                          window='hanning',
                                          nperseg=seg,
                                          nfft=seg,
                                          detrend=False,
                                          scaling='spectrum',
                                          noverlap=noverlap)
                kh = analysisUtils.qkhf(f, 0.028)
                cutoff = np.argmax(f >= 2)
                c = np.sqrt(9.81 * np.tanh(kh[1:cutoff])) / np.sqrt(kh[1:cutoff] / 0.028)
                n = (1 / 2) * (1 + ((2 * kh[1:cutoff]) / (np.sinh(2 * kh[1:cutoff]))))
                cg = c * n  # Wave celerity
                Hs_spec.append(4 * np.sqrt(np.sum(Pxx_eta[1:cutoff] * np.mean(np.diff(f)))))
                F.append(1025 * 9.81 * Pxx_eta[1:cutoff] * cg)  # Wave energy
                
                # Estimate near-bottom orbital velocity with LWT
                omegaj.append(2 * np.pi * f[1:cutoff])
                ubj.append((2 * Pxx_eta[1:cutoff] * 2 * np.pi * f[1:cutoff]) / np.sinh(kh[1:cutoff]))
                
        # Mean wave energy dissipation
        del_x = (self.wave_gauges[-1] - self.wave_gauges[0])
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
    
    def avgFields(self):
        print('Averaging field data...')
        x = self.fields['x'][0]
        z = self.fields['z'][0]
        maxBin = round(self.wave_gauges[-1] / 0.094, 0)
        
        # Bin epsilon along x and z
        eps_avg = np.mean(self.fields['epsilon'])  # Avg. TKE dissipation in time
        bins = stats.binned_statistic_2d(z, x, eps_avg, 'mean', bins=[26, int(maxBin)])
        epsg = bins.statistic  # Binned result
        
        # Bin Umag along x and z
        Ux_avg = np.mean(self.fields['Ux'])
        Uy_avg = np.mean(self.fields['Uy'])
        Uz_avg = np.mean(self.fields['Uz'])
        Umag_avg = np.sqrt(Ux_avg**2 + Uy_avg**2 + Uz_avg**2)
        bins = stats.binned_statistic_2d(z, x, Umag_avg, 'mean', bins=[26, int(maxBin)])
        Umagg = bins.statistic  # Binned result
        
        # Bin k along x and z
        k_avg = np.mean(self.fields['k'])  # Avg. TKE dissipation in time
        bins = stats.binned_statistic_2d(z, x, k_avg, 'mean', bins=[26, int(maxBin)])
        kg = bins.statistic  # Binned result
        
        # Create data frame and return results
        avg_fields = pd.DataFrame()
        avg_fields['eps'] = pd.Series(list(epsg))
        avg_fields['Umag'] = pd.Series(list(Umagg))
        avg_fields['k'] = pd.Series(list(kg))
        
        # The x and z bins need some prep... this code reshapes them into (20, 80) arrays.
        xBins = np.tile(bins.y_edge[:-1].reshape(-1, int(maxBin)), (26, 1))
        zBins = np.tile(bins.x_edge[:-1].reshape(-1, 26), (int(maxBin), 1))
        avg_fields['xBins'] = pd.Series(list(xBins))
        avg_fields['zBins'] = pd.Series(list(zBins))
        
        return avg_fields
    
        # NOTE: to 'reconstruct' avg_fields values, a series of lists -> array,
        # use np.stack(), e.g. np.stack(avg_fields['eps'])
    
    def feddersenDissipation(self, epsilon, epsj):
        '''
        Calculate surf-zone averaged dissipation (Feddersen, 2012)
        Inputs:
            params: contained in 'self':
                h: water depth
            epsilon: averaged values from avgFields (as a series of lists)
            epsj: wave energy flux across the model (remember to use [:-1])
        
        Outputs:
            eps_norm: profile of normalized dissipation; e.g., eps_norm(z)
        '''
        dFdx = np.trapz(epsj)
        eps_norm = [np.median(epsilon[x]) / (float(self.h)**-1 * dFdx) for x in range(len(epsilon))]
        
        return eps_norm
    
    def hendersonDamping(self, ubr):
        '''
        Calculate Henderson dimensionless damping parameter Lambda
        (Henderson et al. 2017)
        
        Inputs:
            params: contained in 'self':
                d: mean diameter of objects in flow
                deltaS: mean object spacing
                Tw: wave period
                h: water depth
                hc: canopy height
            ubr: near-bottom wave orbital velocity
        
        Outputs:
            Lambda_not: Henderson damping parameter
        '''
        # Estimate the drag coefficient (Lentz et al. 2017)
        a = self.d / (self.deltaS**2)
        kappa = 0.4  # Von Karman's constant
        Pi = 0.2  # Cole's wake strength
        z_not = self.hc / self.h
        Cd = kappa**2 * (np.log(self.h / z_not) + (Pi - 1))**-2
        Lambda_not = (Cd * a * ubr * self.Tp) / (4 * np.pi)
        
        return Lambda_not
        