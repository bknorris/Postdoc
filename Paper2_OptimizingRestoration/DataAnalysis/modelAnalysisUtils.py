# -*- coding: utf-8 -*-
"""
Class for processing and saving OpenFOAM model freeSurface and field data
for Paper2_OptimizingRestoration

BKN - USGS PCMSC 2022
"""
import os
import pandas as pd
import numpy as np
import readModelUtils

# init should contain model paths "/Model/postProcessing/..."
# init should take in model settings from saveModelFiles.py
# Create timestep
# Load freeSurf
# Load U
# Load p_rgh
# Load k
# Load eps


class modelAnalysisUtils:
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
    postProcessing = '\\Model\\postProcessing\\'
    
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
        wave_gauges = np.linspace(x_bnd[0],x_bnd[1],10) # 10 evenly spaced WGs
        
        # Get the model folders from postProcessing
        full_path = dest_path + "\\" + model_name + modelAnalysisUtils.postProcessing + "freeSurface\\"
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
        
    def modelSurfAnalysis(self,timesteps):
        # Redefine full_path for freeSurface postProcessing folder
        surf_path = self.dest_path + "\\" + self.model_name + modelAnalysisUtils.postProcessing + "freeSurface\\"
        for folder in self.model_folders:
            full_path = surf_path + folder + "\\"
            file = next(os.walk(full_path))[2]
            file_text = readModelUtils.read_model_file(full_path + file[0])
            
            # Find unique values of x with a tolerance of 1e-4
            x_sort = np.unique(file_text['x'].round(decimals=4),return_index=True)
            y_ind = x_sort[1]
            z_sort = file_text['z'][x_sort[1]]
            
            # Accumulate array of z based on y_ind
            
            ### ACCUMARRAY DOESNT WORK! NEED EDITS!
            output = [x_sort[0], readModelUtils.accumarray(y_ind, z_sort)]
            data = pd.DataFrame()
            count = 1
            for i in range(0,len(self.wave_gauges)):
                inst_x_id = np.argmin(np.abs(output[0])-self.wave_gauges[i])
                data = pd.DataFrame(output[1][inst_x_id], columns=[f'WG{count}'])
                