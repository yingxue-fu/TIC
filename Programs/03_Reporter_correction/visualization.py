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
import seaborn as sns
from scipy.stats import norm
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
        

def logfc_hist(df, Xrange, figsize, figtitle, saveFile):
    # draw histogram
    axes = df.hist(bins=20, density=True, figsize=figsize, sharex=True, sharey=True, range=Xrange, grid=False)
    axes = axes.flatten()
    for i,x in enumerate(axes):
        # Despine
        x.spines['right'].set_visible(False)
        x.spines['top'].set_visible(False)
        # Set x-axis label
        x.set_xlabel("log2FC", size=16)
        # Set y-axis label
        x.set_ylabel("Density", size=16)
        x.set_title(x.get_title(), fontsize=16)
        x.tick_params(axis='both', which='major', labelsize=14)
        
    # set title 
    plt.suptitle(figtitle, ha='center', fontsize=18)
    
    plt.tight_layout()
    if saveFile != 0:
        plt.savefig(saveFile)



def logfc_density(df, IonType, ax=None):
    
    channels = df.columns.values.tolist()
    
    if ax is None:
        ax = plt.gca()
    
    # Iterate through channels
    for channel in channels:
        # Subset to the channel
        df_sub = normTrimed(df[channel], 2)
        # Draw the density plot
        sns.kdeplot(df_sub, label = channel, ax=ax)
        
    # Plot formatting
    ax.legend(prop={'size': 10}, title = 'Channel', loc='upper left')
    ax.set_title('Density plot of log2FC ('+IonType+')', fontsize=16)
    ax.set_xlabel('log2FC', fontsize=16)
    ax.set_ylabel('Density', fontsize=16)
    
    return ax
    

# 
def normTrimed(logfc, SDfold):
    # calculate mean and sd
    mu, std = norm.fit(logfc)
    # get a subset
    logfc_sub = logfc[(logfc > (mu-SDfold*std)) & (logfc < (mu+SDfold*std))]
    
    return logfc_sub 



def averSDcompare(reptSD, tmtcSD, ax=None):
    
    if ax is None:
        ax = plt.gca()
    
    x = np.linspace(-4*reptSD,4*reptSD,100)
    ax.plot(x, norm.pdf(x, 0, reptSD), c='r')
    x = np.linspace(-4*tmtcSD,4*tmtcSD,100)
    ax.plot(x, norm.pdf(x, 0, tmtcSD), c='b')
    
    # Plot formatting
    ax.legend(['Reporter (SD='+str(np.around(reptSD,3))+')', 
                '   TMTc    (SD='+str(np.around(tmtcSD,3))+')'], 
               prop={'size': 10},loc='upper left', title='Ion type')
    ax.set_title('Density plot of log2FC (using average SD)', fontsize=16)
    ax.set_xlabel('log2FC', fontsize=16)
    ax.set_ylabel('Density', fontsize=16)
    
    return ax
    


def SDcomparePlots(rept_df, tmtc_df, reptSD, tmtcSD, saveFile):
    
    fig = plt.figure(figsize = (21,6))
    ax1 = fig.add_subplot(1, 3, 1)
    ax2 = fig.add_subplot(1, 3, 2, sharex = ax1, sharey = ax1)
    ax3 = fig.add_subplot(1, 3, 3)
    
    # reporter log2FC density plot
    logfc_density(rept_df, 'Reporter', ax=ax1)
    # TMTc log2FC density plot
    logfc_density(tmtc_df, 'TMTc', ax=ax2)
    # average SD comparison
    averSDcompare(reptSD, tmtcSD, ax=ax3)

    if saveFile != 0:
        plt.savefig(saveFile)