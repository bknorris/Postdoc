# -*- coding: utf-8 -*-
"""
Modules to support plotting OpenFOAM postProcessing output files

BKN - USGS PCMSC 2022
"""
import numpy as np
import seaborn as sns
import matplotlib
from matplotlib import pyplot as plt
from mpl_toolkits.mplot3d import Axes3D


def quick_plot(x, y, z, value, labels=['x', 'y', 'z']):
    # axes instance
    fig = plt.figure(figsize=(6, 6))
    ax = Axes3D(fig, auto_add_to_figure=False)
    fig.add_axes(ax)
    
    # get colormap from seaborn
    cmap = matplotlib.colors.ListedColormap(sns.color_palette("husl", 256).as_hex())
    
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


def cmap_discretize(cmap, N):
    """Return a discrete colormap from the continuous colormap cmap.
    
        cmap: colormap instance, eg. cm.jet.
        N: number of colors.
    """
    if type(cmap) == str:
        cmap = matplotlib.cm.get_cmap(cmap)
    colors_i = np.concatenate((np.linspace(0, 1., N), (0., 0., 0., 0.)))
    colors_rgba = cmap(colors_i)
    indices = np.linspace(0, 1., N + 1)
    cdict = {}
    for ki, key in enumerate(('red', 'green', 'blue')):
        cdict[key] = [(indices[i], colors_rgba[i - 1, ki], colors_rgba[i, ki]) for i in range(N + 1)]
    # Return colormap object.
    return matplotlib.colors.LinearSegmentedColormap(cmap.name + "_%d" % N, cdict, 1024)
