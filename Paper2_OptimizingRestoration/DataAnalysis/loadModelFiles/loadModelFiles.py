# -*- coding: utf-8 -*-
"""
Modules for loading and saving OpenFOAM freeSurface and field data
for Paper2_OptimizingRestoration models.

BKN - USGS PCMSC 2022
"""
import os
import pandas as pd
import numpy as np
import pickle
import analysisUtils
from time import time


class loadModelFiles:
    '''
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
    
    def __init__(self, source_path, model_name, wave_period):
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
        '''
        Update 4/21/22: Sometimes the folders were zipped recursively,
        i.e., Scenario_X/Scenario_X/...
        
        Some folders were renamed,
        e.g., Scenario_X/Scenario_Y/...
        
        Catch these exceptions.
        '''
        surf_path = source_path + "\\" + model_name + "\\Model\\postProcessing\\freeSurface\\"
        bathy_path = source_path + "\\" + model_name + "\\Model\\postProcessing\\bathySample\\surface\\"
        folders = next(os.walk(surf_path), (None, None, []))[1]
        if folders is None:
            subfolder = analysisUtils.find_files(source_path + "\\" + model_name + '\\', 'Scenario_*', 'folders')
            surf_path = source_path + "\\" + model_name + "\\" + subfolder[0] + "\\Model\\postProcessing\\freeSurface\\"
            bathy_path = source_path + "\\" + model_name + "\\" + subfolder[0] + "\\Model\\postProcessing\\bathySample\\surface\\"
            folders = next(os.walk(surf_path), (None, None, []))[1]
            
        self.surf_path = surf_path
        self.bathy_path = bathy_path
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
        
        print('Loading model free surface data...')
        output = []
        for folder in self.model_folders:
            full_path = self.surf_path + folder + "\\"
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
        print(f'Data loaded in: {executionTime:.2f} minutes')
        return freeSurf
    
    def loadFields(self, timesteps):
        # Set timer
        startTime = time()

        # Load the field data
        print('Loading model field data...')
        output = []
        for folder in self.model_folders:
            full_path = self.bathy_path + folder + "\\"
            file = next(os.walk(full_path))[2]
            
            # Read in text files
            fid = file.index('U_sampledSurface_bathymetrySample.raw')
            U_text = analysisUtils.read_model_file(full_path + file[fid])
            fid = file.index('p_rgh_sampledSurface_bathymetrySample.raw')
            p_rgh_text = analysisUtils.read_model_file(full_path + file[fid])
            fid = file.index('k_sampledSurface_bathymetrySample.raw')
            k_text = analysisUtils.read_model_file(full_path + file[fid])
            fid = file.index('epsilon_sampledSurface_bathymetrySample.raw')
            epsilon_text = analysisUtils.read_model_file(full_path + file[fid])
            
            # Create arrays from text data
            x_vals = np.array(U_text['x'])
            y_vals = np.array(U_text['y'])
            z_vals = np.array(U_text['z'])
            Ux_vals = np.array(U_text['U_x'])
            Uy_vals = np.array(U_text['U_y'])
            Uz_vals = np.array(U_text['U_z'])
            p_rgh = np.array(p_rgh_text['p_rgh'])
            k = np.array(k_text['k'])
            epsilon = np.array(epsilon_text['epsilon'])
            
            # Crop the data file by the area defined by the user
            crop_x = np.where(np.logical_and(x_vals >= self.x_bnd[0], x_vals <= self.x_bnd[1]))
            crop_z = np.where(np.logical_and(z_vals >= self.z_bnd[0], z_vals <= self.z_bnd[1]))
            crop_xz = np.intersect1d(crop_x, crop_z)
            
            # Output data into dicts
            keys = ['x', 'y', 'z', 'Ux', 'Uy', 'Uz', 'p_rgh', 'k', 'epsilon']
            values = [x_vals[crop_xz], y_vals[crop_xz], z_vals[crop_xz],
                      Ux_vals[crop_xz], Uy_vals[crop_xz], Uz_vals[crop_xz],
                      p_rgh[crop_xz], k[crop_xz], epsilon[crop_xz]]
            
            output.append(dict(zip(keys, values)))
            del U_text, x_vals, y_vals, z_vals,
            Ux_vals, Uy_vals, Uz_vals,
            p_rgh, k, epsilon  # clear memory
            
        # Create data frame from list of dicts (output)
        fields = pd.DataFrame(output)
        fields.insert(loc=0, column='timestep', value=timesteps)
        
        # Print elapsed time and return U
        executionTime = (time() - startTime) / 60
        print(f'Data loaded in: {executionTime:.2f} minutes')
        return fields
    
    def saveData(self, data, kind, file_path, version):
        '''
        Inputs:
            data: data to be saved as a binary file
            kind: string type of data file to be saved [e.g., freeSurf, fields]
            file_path: path to save out binary file
            version: string value for version no. [e.g., V1, V2, etc.]
        '''
        
        print('Saving ' + kind + ' file...')
        save_file = False
        output_file = self.model_name + "_" + kind + "_" + version + ".dat"
        if os.path.isfile(file_path + output_file):
            print('File already exists in destination directory!')
            result = input('Overwrite? (y/n) ')
            if result == 'y':
                save_file = True
                
        else:
            save_file = True
            
        if save_file:
            file_obj = open(file_path + output_file, mode='wb')
            pickle.dump(data, file_obj)
            file_obj.close()
