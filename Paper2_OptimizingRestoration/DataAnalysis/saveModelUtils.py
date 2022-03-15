# -*- coding: utf-8 -*-
"""
Class for processing and saving OpenFOAM model freeSurface and field data
for Paper2_OptimizingRestoration

BKN - USGS PCMSC 2022
"""
import pandas

# init should contain model paths "/Model/postProcessing/..."
# init should take in model settings from saveModelFiles.py
# Create timestep
# Load freeSurf
# Load U
# Load p_rgh
# Load k
# Load eps

class modelUtils():
    ''' 
    Processing utilities for OpenFOAM model output files.
    
    Inputs:
        source_path: path of input file. Must be the base model folder (typically ".../ModelRuns/Scenarios/")
        dest_path: path where processed data file will be saved
        file_name: model folder name (typically "Scenario_XX")
        wave_period: integer value for the wave period of each model (used to
                     determine which processing step to run)
        
    Outputs:
        pandas array containing model freeSurface and field data saved in BIN format
    '''
    postProcessing = 'Model/postProcessing/'
    
    def __init__(self,source_path,dest_path,file_name,wave_period):
        if wave_period == 5:
            x_bnd = (0, 2.6)  # lower/upper
            z_bnd = (-5, -0.005)  # lower/upper
            model_timestep = 0.5  # timestep in s
            wave_gauges = []  # wave gauge positions for freeSurface
        elif wave_period == 10: ### CHANGE THESE VALUES!!! ###
            x_bnd = (0, 2.6)  # lower/upper
            z_bnd = (-5, -0.005)  # lower/upper
            model_timestep = 0.5  # timestep in s
            wave_gauges = []  # wave gauge positions for freeSurface
        elif wave_period == 20:
            x_bnd = (0, 2.6)  # lower/upper
            z_bnd = (-5, -0.005)  # lower/upper
            model_timestep = 0.5  # timestep in s
            wave_gauges = []  # wave gauge positions for freeSurface
        else:
            x_bnd = (0, 2.6)  # lower/upper
            z_bnd = (-5, -0.005)  # lower/upper
            model_timestep = 0.5  # timestep in s
            wave_gauges = []  # wave gauge positions for freeSurface
        
        self.source_path = source_path
        self.dest_path = dest_path
        self.file_name = file_name
        
    def createTimestep():
        
        
    