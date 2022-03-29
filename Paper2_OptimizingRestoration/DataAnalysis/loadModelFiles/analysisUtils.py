# -*- coding: utf-8 -*-
"""
Modules to support analyzing OpenFOAM postProcessing output files

BKN - USGS PCMSC 2022
"""
import pandas as pd
import numpy as np
import seaborn as sns
from matplotlib import pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.colors import ListedColormap
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


def sub2ind(array_shape, rows, cols):
    return rows * array_shape[1] + cols


def qkhf(w, h):
    '''
    QKHF  Quick explicit calculation of kh in dispersion relationship.
    
    kh = qkhf( w, h )
    
    Hard-wired for MKS units.
    Dean and Dalrymple (1991), p 72.
    
    Input:
        w Angular wave frequency = 2*pi/T where T = wave period [1/s]
        h Water depth [m]
        
    Returns:
        kh = wavenumber * depth [ ]
    
        Either w or h can be a vector, but not both.
    
    Chris Sherwood, USGS (March 17, 1999)
    '''
    D1 = 0.6666666666
    D2 = 0.3555555555
    D3 = 0.1608465608
    D4 = 0.0632098765
    D5 = 0.0217540484
    D6 = 0.0065407983
    G = 9.80665
    
    y = (w**2) * (h / G)
    
    # Calculate polynomial on bottom line
    kh2 = 1.0 + y * (D1 + y * (D2 + y * (D3 + y * (D4 + y * (D5 + y * D6)))))
    
    # Calculate final term
    kh2 = y * y + y / kh2
    
    kh = np.sqrt(kh2)
    return kh
    
    
def quick_plot(x, y, z, value, labels=['x', 'y', 'z']):
    # axes instance
    fig = plt.figure(figsize=(6, 6))
    ax = Axes3D(fig, auto_add_to_figure=False)
    fig.add_axes(ax)
    
    # get colormap from seaborn
    cmap = ListedColormap(sns.color_palette("husl", 256).as_hex())
    
    # plot
    sc = ax.scatter(x, y, z, s=40, c=value, marker='o', cmap=cmap, alpha=1)
    ax.set_xlabel(labels[0])
    ax.set_ylabel(labels[1])
    ax.set_zlabel(labels[2])
    
    # legend
    plt.legend(*sc.legend_elements(), bbox_to_anchor=(1.05, 1), loc=2)

    
def imagesc(x, y, image, fig_size=(6, 4), caxis=None, c_map='viridis'):
    if(caxis is None):
        caxis = [np.min(image), np.max(image)]
        
    def extents(f):
        delta = f[1] - f[0]
        return [f[0] - delta / 2, f[-1] + delta / 2]
    fig = plt.figure(figsize=fig_size)
    ax = fig.add_axes([0.1, 0.1, 0.8, 0.8])
    axes = ax.imshow(image, aspect='auto', extent=extents(x) + extents(y),
                     vmin=caxis[0], vmax=caxis[1], origin='lower', cmap=c_map)
    plt.colorbar(axes)

    