# -*- coding: utf-8 -*-
"""
Modules to support analyzing OpenFOAM postProcessing output files

BKN - USGS PCMSC 2022
"""
import pandas as pd
import numpy as np
from itertools import product

def read_model_file(file_path):
    # Read model data file and skip headerlines commented with "#"
    with open(file_path) as f:
        line = f.readline()
        count = 0
        while line.startswith('#'):
            prev_line = line
            line = f.readline()
            count += 1

    header = prev_line.strip().lstrip('# ').split()
    file_text = pd.read_csv(file_path, delimiter='\s+',
                            names=header,
                            skiprows=count)
    return file_text


def uniquetol(data, tol=1e-6):
    '''
    Performs unique comparison using a tolerance.
    An approximation of MATLAB's "uniquetol"
    
    Usage
    -----
    result, IA, IC = uniquetol(data, tol)

    Parameters
    ----------
    data : input data
    tol : tolerance (default 1e-6)

    Returns
    -------
    result : unique data based on tolerance
    IA : Index vector such that result = data(IA)
    IC : Index vector such that data ~ result(IC)

    '''
    idx = np.argsort(data.flat)
    d_flat = data.flat[idx]
    d = np.append(True, np.diff(d_flat))
    IA = idx[d > tol]
    IC = np.full(len(d_flat), np.nan)
    for ids, val in enumerate(d_flat[d > tol]):
        IC[np.where(d_flat == val)] = idx[ids]
        
    result = data.flat[IA]
    return result, IA, IC


def runningmean(a, n):
    return np.convolve(a, np.ones(n, dtype='float'),
                       'same') / np.convolve(np.ones(len(a)), np.ones(n), 'same')

def accumarray(subs, val):
    df = pd.DataFrame({"y": val, "x": subs})
    return pd.pivot_table(df, values='y', index='x', aggfunc='mean')