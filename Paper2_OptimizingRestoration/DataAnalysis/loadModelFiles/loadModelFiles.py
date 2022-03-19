# -*- coding: utf-8 -*-
"""
Class for loading and saving OpenFOAM model freeSurface and field data
for Paper2_OptimizingRestoration

BKN - USGS PCMSC 2022
"""
import os
from time import time
import pandas as pd
import numpy as np
import analysisUtils

# init should contain model paths "/Model/postProcessing/..."
# init should take in model settings from saveModelFiles.py
# Create timestep
# Load freeSurf
# Load U
# Load p_rgh
# Load k
# Load eps


class loadModelFiles:
    '''
    Processing utilities for OpenFOAM model output files. Note this script
    is designed to work with the Paper2_OptimizingRestoration model files only!
    
    Inputs:
        source_path: path of input file. Must be the base model folder (typically ".../ModelRuns/Scenarios/")
        dest_path: path where processed data file will be saved
        file_name: model folder name (typically "Scenario_XX")
        wave_period: integer value for the wave period of each model (used to
                     determine which processing step to run)
        
    Outputs:
        pandas array containing model freeSurface and field data saved in BIN format
    '''
    
    def __init__(self, source_path, dest_path, model_name, wave_period):
        if wave_period == 5:
            x_bnd = (0, 2.6)  # lower/upper
            z_bnd = (-0.028, -0.001)  # lower/upper
        elif wave_period == 10:
            x_bnd = (0, 5.2)  # lower/upper
            z_bnd = (-0.028, -0.001)  # lower/upper
        elif wave_period == 20:
            x_bnd = (0, 10.4)  # lower/upper
            z_bnd = (-0.028, -0.001)  # lower/upper
        else:
            x_bnd = (0, 21.6)  # lower/upper
            z_bnd = (-0.028, -0.001)  # lower/upper
        wave_gauges = np.linspace(x_bnd[0], x_bnd[1], 10)  # 10 evenly spaced WGs
        
        # Get the model folders from postProcessing
        full_path = dest_path + "\\" + model_name + "\\Model\\postProcessing\\freeSurface\\"
        folders = next(os.walk(full_path), (None, None, []))[1]
        
        self.source_path = source_path
        self.dest_path = dest_path
        self.model_name = model_name
        self.model_folders = folders[1:]  # ignore first directory (is 0 time)
        self.x_bnd = x_bnd
        self.z_bnd = z_bnd
        self.wave_gauges = wave_gauges
        
    def createTimestep(self):
        timesteps = [float(x) for x in self.model_folders]
        timesteps.sort()  # just in case
        return timesteps
        
    def loadFreeSurface(self, timesteps):
        # Set timer
        startTime = time()
        
        # Redefine full_path for freeSurface postProcessing folder
        surf_path = self.dest_path + "\\" + self.model_name + "\\Model\\postProcessing\\freeSurface\\"
        output = []
        for folder in self.model_folders:
            full_path = surf_path + folder + "\\"
            file = next(os.walk(full_path))[2]
            file_text = analysisUtils.read_model_file(full_path + file[0])
            
            # Find unique values of x with a tolerance of 1e-4
            x_vals = np.array(file_text['x'])
            z_vals = np.array(file_text['z'])
            x_sort, ix, iy = analysisUtils.uniquetol(x_vals, 1e-4)
            
            # Sort z values based on uniquetol of x
            z_sort = z_vals[ix]
            z_mean = analysisUtils.runningmean(z_sort, 5)
            
            # Output data into dicts for each wave_gauge
            wgs = [f'WG{x}' for x in list(range(1, len(self.wave_gauges) + 1))]
            d = []
            for i, gauge in enumerate(wgs):
                inst_x_id = np.argmin((np.abs(x_sort - self.wave_gauges[i])))
                d.append(z_mean[inst_x_id])
                
            output.append(dict(zip(wgs, d)))
            del file_text, x_vals, z_vals, x_sort, z_mean  # clear memory
        
        # Create data frame from list of dicts (output)
        freeSurf = pd.DataFrame(output)
        freeSurf.insert(loc=0, column='TimeStep', value=timesteps)
        
        # Print elapsed time and return freeSurf
        executionTime = (time() - startTime) / 60
        print(f'Data loaded in {executionTime:.2f} minutes')
        return freeSurf
    
    def loadFields(self, timesteps):
        # Set timer
        startTime = time()

        # Redefine full_path for bathySample postProcessing folder
        bathy_path = self.dest_path + "\\" + self.model_name + "\\Model\\postProcessing\\bathySample\\surface\\"
        
        print('Loading U field data')
        output = []
        for folder in self.model_folders:
            full_path = bathy_path + folder + "\\"
            file = next(os.walk(full_path))[2]
            fid = file.index('U_sampledSurface_bathymetrySample.raw')
            file_text = analysisUtils.read_model_file(full_path + file[fid])
            x_vals = np.array(file_text['x'])
            y_vals = np.array(file_text['y'])
            z_vals = np.array(file_text['z'])
            Ux_vals = np.array(file_text['U_x'])
            Uy_vals = np.array(file_text['U_y'])
            Uz_vals = np.array(file_text['U_z'])
            
            # Crop the data file by the area defined by the user
            crop_x = np.where(np.logical_and(x_vals >= self.x_bnd[0], x_vals <= self.x_bnd[1]))
            crop_z = np.where(np.logical_and(z_vals >= self.z_bnd[0], z_vals <= self.z_bnd[1]))
            crop_xz = np.intersect1d(crop_x, crop_z)
            
            # Output data into dicts
            keys = ['x', 'y', 'z', 'Ux', 'Uy', 'Uz']
            values = [x_vals[crop_xz], y_vals[crop_xz], z_vals[crop_xz], Ux_vals[crop_xz], Uy_vals[crop_xz], Uz_vals[crop_xz]]
            
            output.append(dict(zip(keys, values)))
            del file_text, x_vals, y_vals, z_vals, Ux_vals, Uy_vals, Uz_vals  # clear memory
            
        # Create data frame from list of dicts (output)
        U = pd.DataFrame(output)
        U.insert(loc=0, column='TimeStep', value=timesteps)
        
        # Print elapsed time and return U
        executionTime = (time() - startTime) / 60
        print(f'Data loaded in {executionTime:.2f} minutes')
        return U
    
        ### NOTE: CHECK OUTPUT OF U in SEABORN 3D SCATTER PLOT!