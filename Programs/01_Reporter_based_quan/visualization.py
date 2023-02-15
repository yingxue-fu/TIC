#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Aug 14 20:54:19 2022

@author: yfu
"""

import math
import numpy as np, pandas as pd
import matplotlib.pyplot as plt
import matplotlib
matplotlib.use('Agg')
# %matplotlib inline

plt.rcParams.update({'font.sans-serif':'Arial'})

# sample ratio histogram 
def quan_mtx_histogram(df, species, theoretical_sample_ratios, params, level, saveFile):
    
    mtx = df.loc[:, df.columns.str.contains(pat = 'sig1')]
    
    # read parameters
    used_channels = params["tmt_channels_used"].split(";")
        
    # normalize sample ratios
    mtx_norm = mtx.div(mtx.apply(np.sum, axis='columns'), axis='index')
    r1_index = theoretical_sample_ratios == 1
    normalize_factor = mtx_norm.loc[:, r1_index].mean().mean()
    ratios_norm = mtx_norm.div(normalize_factor) 
    
    # define No. of rows and columns for plot
    n_channels = len(used_channels)
    n_row = math.ceil(math.sqrt(n_channels))
    n_col = math.ceil(n_channels / n_row)
    
    # generate figure
    figure_size = (4*n_col, 3*n_row)
    fig, axs = plt.subplots(n_row, n_col, sharey=False, figsize = figure_size)
    axs = axs.flatten()
    for x in range(n_channels):
        axs[x].hist(ratios_norm.iloc[:,x], bins=30, density=False)
        y_hist, _ = np.histogram(ratios_norm.iloc[:,x], bins=30, density=False)
        axs[x].vlines(x=theoretical_sample_ratios[x], ymin=0, ymax=y_hist.max(), linewidth=2, color='r')
        axs[x].set_xlabel('Ratio', fontsize=16)
        axs[x].set_ylabel('# of '+level+'s', fontsize=16)
        axs[x].set_title(used_channels[x], fontsize=16)
        axs[x].set_ylim(0, y_hist.max()+30)
        for tick in axs[x].yaxis.get_major_ticks():
            tick.label.set_fontsize(14) 
        for tick in axs[x].xaxis.get_major_ticks():
            tick.label.set_fontsize(14) 
    fig.suptitle('Reporter-based quantification ('+ species + ', ' + level + '-level)', fontsize = 18)
    
    # save figure 
    fig.tight_layout()
    if saveFile != 0:
        fig.savefig(saveFile)
        

# sample ratio barplot
def quan_mtx_barplot(df, species, theoretical_sample_ratios, params, level, saveFile):
    
    mtx = df.loc[:, df.columns.str.contains(pat = 'sig1')]
    
    # read parameters
    used_channels = params["tmt_channels_used"].split(";")
    
    # normalize sample ratios
    mtx_norm = mtx.div(mtx.apply(np.sum, axis='columns'), axis='index')
    r1_index = theoretical_sample_ratios == 1
    normalize_factor = mtx_norm.loc[:, r1_index].mean().mean()
    ratios_norm = mtx_norm.div(normalize_factor) 
    
    # calculate mean and sd
    mean_channels = np.mean(ratios_norm, axis=0)
    sd_channels = np.std(ratios_norm, axis=0)
    
    # generate figure
    fig, ax = plt.subplots(figsize = (10, 4))
    n_channels = len(used_channels)
    bars = ax.bar(np.arange(n_channels), mean_channels, yerr=[[0]*n_channels,sd_channels], edgecolor='black', align='center', alpha=0.5, ecolor='gray', capsize=5)
    ax.bar_label(bars, labels=np.around(mean_channels,2), fontsize = 14)
    ax.set_ylabel('Sample Ratio', fontsize = 14)
    ax.set_xticks(np.arange(n_channels))
    x_labels = pd.Series(used_channels) + '_' + pd.Series(theoretical_sample_ratios).astype(str)
    ax.set_xticklabels(x_labels, rotation='vertical', fontsize = 14)
    ax.tick_params(axis='both', which='major', labelsize=14)
    ax.set_title('Reporter-based quantification ('+ species + ', ' + level + '-level)', fontsize = 18)
    ax.yaxis.grid(False)
    
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)

    # save figure 
    fig.tight_layout()
    if saveFile != 0:
        fig.savefig(saveFile)
        
    
