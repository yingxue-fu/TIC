#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Aug 14 20:54:19 2022

@author: yfu
"""

import math, os
import numpy as np
import matplotlib.pyplot as plt
from get_tmtc_groups import get_combined_ratios
import matplotlib
from matplotlib.backends.backend_pdf import PdfPages
matplotlib.use('Agg')
# %matplotlib inline
#import seaborn as sns

plt.rcParams.update({'font.sans-serif':'Arial'})

# Count # of extracted TMTc peaks 
def n_TMTc_plot(dfId_TMTc, TMTc_peak, figTitle, ax=None):
    # TMTc_peak = "mono" or "all"
    
    if ax is None:
        ax = plt.gca()
    
    if TMTc_peak == 'mono':
        n_TMTc_peaks_stat = dfId_TMTc['n_mono_TMTc_peaks'].value_counts().sort_index(ascending=True)
    elif TMTc_peak == 'all':
        n_TMTc_peaks_stat = dfId_TMTc['n_TMTc_peaks'].value_counts().sort_index(ascending=True)
    
    # for mono TMTc peaks
    x = n_TMTc_peaks_stat.index
    bar = ax.bar(x, n_TMTc_peaks_stat, 0.8, edgecolor='black', align='center', alpha=0.5, capsize=10)
    ax.bar_label(bar, fontsize = 14)
    ax.set_xticks(x)
    ax.set_xlabel("# of "+TMTc_peak+" TMTc peaks", fontsize = 16)
    ax.set_ylabel('# of PSMs', fontsize = 16)
    ax.set_title(figTitle, fontsize=16)
    ax.tick_params(axis='both', which='major', labelsize=14)
    
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    
    return ax

        
def n_TMTc_fractions(dfId_TMTc, TMTc_peak, interDir, saveFile):
    
    pdf = PdfPages(saveFile)
    
    files = dfId_TMTc["Files"].unique()
    
    for file in files:
        dfId_TMTc_sub = dfId_TMTc[dfId_TMTc['Files'] == file]
        
        subDir = os.path.join(interDir, 'Extracted_TMTc_peaks')
        os.makedirs(subDir, exist_ok=True)
        dfId_TMTc_sub.to_csv(os.path.join(subDir, os.path.basename(file)+"_extracted_TMTc_peaks.txt"), sep="\t", index=False)
        
        fig, ax = plt.subplots(1, 1, sharey=False, figsize = (6,4))
        n_TMTc_plot(dfId_TMTc_sub, TMTc_peak, os.path.basename(file), ax=ax)
        
        fig.tight_layout()
        pdf.savefig(fig)
        
    pdf.close()



def opt_metrics_plots(tmtcQuan, saveFile):
    
    # 
    tmtcQuan['n_mono_TMTc_peaks_grp'] = tmtcQuan['n_mono_TMTc_peaks'].astype(str)
    tmtcQuan['-log10(SSD)'] = -np.log10(tmtcQuan['Diff_square_sum'].to_numpy())
    tmtcQuan['log2(total.TMTc.intensity)'] = np.log2(tmtcQuan['total.TMTc.intensity'].to_numpy())
    
    fig, axs = plt.subplots(1, 3, sharey=False, figsize = (15,5))
    axs = axs.flatten()
    
    # SSD hist
    axs[0].hist(tmtcQuan['-log10(SSD)'], bins=50, density=True)
    y_hist, _ = np.histogram(tmtcQuan['-log10(SSD)'], bins=50, density=True)
    axs[0].vlines(x=-np.log10(0.005), ymin=0, ymax=y_hist.max(), linewidth=1, color='r')
    axs[0].text(-np.log10(0.005)*1.1, y_hist.max()*0.95, 'SSD < 0.005', fontsize=14)
    axs[0].text(-np.log10(0.005)*1.1, y_hist.max()*0.8, '(Sum of squared differences)', fontsize=14)
    axs[0].set_xlabel('-log10(SSD)', fontsize=16)
    axs[0].set_ylabel('Density', fontsize=16)
    axs[0].set_title('Distribution of SSD', fontsize=16)
    axs[0].tick_params(axis='both', which='major', labelsize=14)
    
    # relation between SSD and TMTc intensity sum
    axs[1].scatter(tmtcQuan['log2(total.TMTc.intensity)'], tmtcQuan['-log10(SSD)'])
    axs[1].hlines(y=-np.log10(0.005), xmin=tmtcQuan['log2(total.TMTc.intensity)'].min(), xmax=tmtcQuan['log2(total.TMTc.intensity)'].max(), linewidth=1, color='r')
    axs[1].set_ylabel('-log10(SSD)', fontsize=16)
    axs[1].set_xlabel('log2(total.TMTc.intensity)', fontsize=16)
    axs[1].set_title('SSD vs. TMTc peak intensity', fontsize=16)
    axs[1].tick_params(axis='both', which='major', labelsize=14)
       
    # relation between SSD and # of TMTc peaks
    grp_order = list(map(str, np.sort(tmtcQuan['n_mono_TMTc_peaks'].astype(int).unique())))
    #sns.violinplot(y='-log10(SSD)', x= 'n_mono_TMTc_peaks_grp', data=tmtcQuan, order=grp_order, orient='v', ax=axs[2])
    axs[2].set_xlabel('# of mono TMTc peaks', fontsize=16)
    axs[2].set_ylabel('-log10(SSD)', fontsize=16)
    axs[2].set_title('SSD vs. # of mono TMTc peaks', fontsize=16)
    axs[2].tick_params(axis='both', which='major', labelsize=14)
    axs[2].hlines(y=-np.log10(0.005), xmin=-0.5, xmax=4.5, linewidth=1, color='r')

    # save figure 
    fig.tight_layout()
    if saveFile != 0:
        fig.savefig(saveFile)
    
    
# sample ratio histogram 
def TMTc_quan_histogram(mtx, species, theoretical_sample_ratios, params, level, saveFile):
    
    # read parameters
    TMT_ver = params["tmt_version"]
    used_channels = params["tmt_channels_used"].split(";")
    
    # calculate theoretical values
    TMTc_group_ratio = get_combined_ratios(theoretical_sample_ratios, used_channels, TMT_ver)
    TMTc_group_ratio_num = np.array(TMTc_group_ratio['Ratio'])
        
    # normalize sample ratios
    mtx_norm = np.apply_along_axis(lambda x: x/sum(x), 1, mtx)
    theoretical_sample_ratios_norm = theoretical_sample_ratios / sum(theoretical_sample_ratios)
    #
    r1_index = theoretical_sample_ratios == 1
    normalize_factor = theoretical_sample_ratios_norm[r1_index].mean()
    #
    ratios_norm = mtx_norm / normalize_factor
    
    # define No. of rows and columns for plot
    n_TMTc_groups = TMTc_group_ratio.shape[0]
    n_row = math.ceil(math.sqrt(n_TMTc_groups))
    n_col = math.ceil(n_TMTc_groups / n_row)
   
    # calculate CV
    sd_TMTc_groups = np.around(np.std(ratios_norm, axis=0), 2)
    sd_texts = ["SD = " + sd.astype(str) for sd in sd_TMTc_groups]
    
    # generate figure
    figure_size = (4*n_col, 4*n_row)
    fig, axs = plt.subplots(n_row, n_col, sharey=False, figsize = figure_size)
    axs = axs.flatten()
    for x in range(n_TMTc_groups):
        axs[x].hist(ratios_norm[:,x], bins=30, density=False)
        y_hist, _ = np.histogram(ratios_norm[:,x], bins=30, density=False)
        axs[x].text(TMTc_group_ratio_num[x], y_hist.max()+10, sd_texts[x], fontsize=14)
        axs[x].vlines(x=TMTc_group_ratio_num[x], ymin=0, ymax=y_hist.max(), linewidth=2, color='r')
        axs[x].set_xlabel('Normalized ratio', fontsize=16)
        axs[x].set_ylabel('# of '+level+'s', fontsize=16)
        axs[x].set_title(TMTc_group_ratio['TMT_channels'][x], fontsize=16)
        axs[x].set_ylim(0, y_hist.max()+30)
        for tick in axs[x].yaxis.get_major_ticks():
            tick.label.set_fontsize(14) 
        for tick in axs[x].xaxis.get_major_ticks():
            tick.label.set_fontsize(14) 
    fig.suptitle('TMTc-based quantification (' + species + ', ' + level + '-level)', fontsize = 18)
    
    # save figure 
    fig.tight_layout()
    if saveFile != 0:
        fig.savefig(saveFile)
        

# sample ratio barplot
def TMTc_quan_barplot(mtx, species, theoretical_sample_ratios, params, level, saveFile):
    
    # read parameters
    TMT_ver = params["tmt_version"]
    used_channels = params["tmt_channels_used"].split(";")
    
    # calculate theoretical values
    TMTc_group_ratio = get_combined_ratios(theoretical_sample_ratios, used_channels, TMT_ver)
    
    # normalize sample ratios
    mtx_norm = np.apply_along_axis(lambda x: x/sum(x), 1, mtx)
    theoretical_sample_ratios_norm = theoretical_sample_ratios / sum(theoretical_sample_ratios)
    #
    #r1_index = TMTc_group_ratio["Ratio"] == 1
    r1_index = theoretical_sample_ratios == 1
    normalize_factor = theoretical_sample_ratios_norm[r1_index].mean()
    #normalize_factor = 1 / mtx_norm[:,r1_index].mean().mean()
    #
    ratios_norm = mtx_norm / normalize_factor
    #ratios_norm = mtx_norm * normalize_factor
    
    # calculate mean and sd
    mean_TMTc_groups = np.mean(ratios_norm, axis=0)
    sd_TMTc_groups = np.std(ratios_norm, axis=0)
    
    # generate figure
    fig, ax = plt.subplots(figsize = (8, 7))
    n_TMTc_groups = TMTc_group_ratio.shape[0]
    bars = ax.bar(np.arange(n_TMTc_groups), mean_TMTc_groups, yerr=[[0]*n_TMTc_groups,sd_TMTc_groups], edgecolor='black', align='center', alpha=0.5, ecolor='gray', capsize=5)
    ax.bar_label(bars, labels=np.around(mean_TMTc_groups,2), fontsize = 14)
    ax.set_ylabel('Sample Ratio', fontsize = 14)
    #ax.set_yticks(np.arange(0, 20, 3))
    ax.set_ylim(0, max(mean_TMTc_groups)+max(sd_TMTc_groups)+1)
    for tick in ax.yaxis.get_major_ticks():
        tick.label.set_fontsize(14) 
    ax.set_xticks(np.arange(n_TMTc_groups))
    x_labels = TMTc_group_ratio['TMT_channels'] + '_' + TMTc_group_ratio['Ratio'].astype(str)
    ax.set_xticklabels(x_labels, rotation='vertical', fontsize = 14)
    ax.set_title('TMTc-based quantification (' + species + ', ' + level + '-level)', fontsize = 18)
    ax.yaxis.grid(False)
    
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)

    # save figure 
    fig.tight_layout()
    if saveFile != 0:
        fig.savefig(saveFile)
        