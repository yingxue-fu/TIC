#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep 16 15:20:50 2022

@author: yfu
"""

import os
import numpy as np, pandas as pd
from tmtc_normalization import TMTc_normalization
from visualization import SDcomparePlots
from scipy.stats import norm


def rept_z_correct(reptQuan, tmtcQuan, interDir, params):
    
    print("\n  Calculate log2FC between each channel and the average")
    # For reporter quan
    # log2 transformation 
    reptQuan_log = np.log2(reptQuan)
    # calculate the average
    rept_avg = reptQuan_log.apply(np.mean, axis=1)
    # log2FC between each channel and the average
    reptQuan_logFC = reptQuan_log.sub(rept_avg, axis='index')

    print("\n  Calculate mean and standard deviation (SD) of log2FC for each reporter channel \n")
    rept_logFC_mean_SD = reptQuan_logFC.apply(lambda x: normTrimed_MeanSD(x, 2), axis=0).transpose().reset_index().rename(columns={"index": "Reporter channel", 0: "Mean", 1: "SD"})
    print(rept_logFC_mean_SD.to_string(index=False, col_space=(31,10,10)))
    
    # For TMTc quan
    print("\n  Normalize TMTc quantification by the number of overlapped reporter channels")
    n_channels_TMTc = [tmtc_channels.count('sig') for tmtc_channels in tmtcQuan.columns.to_list()]
    tmtcQuan = tmtcQuan.divide(n_channels_TMTc, axis=1)
    tmtcQuan = TMTc_normalization(tmtcQuan)

    # log2 transformation
    tmtcQuan_log = np.log2(tmtcQuan)
    # calculate the average
    TMTc_avg = tmtcQuan_log.apply(np.mean, axis=1)
    # log2FC between each channel and the average
    tmtcQuan_logFC = tmtcQuan_log.sub(TMTc_avg, axis='index')

    print("\n  Calculate mean and standard deviation (SD) of log2FC for each TMTc channel \n")
    tmtc_logFC_mean_SD = tmtcQuan_logFC.apply(lambda x: normTrimed_MeanSD(x, 2), axis=0).transpose().reset_index().rename(columns={"index": "TMTc channel", 0: "Mean", 1: "SD"})
    print(tmtc_logFC_mean_SD.to_string(index=False, col_space=(31,10,10)))

    # Calculate SD correction factor
    reptSD = rept_logFC_mean_SD['SD'].mean()
    tmtcSD = tmtc_logFC_mean_SD['SD'].mean()
    print("\n  average SD of log2FC (Reporter): {}".format(reptSD))
    print("  average SD of log2FC (TMTc):     {}".format(tmtcSD))
    SD_corrc_factor = tmtcSD / reptSD
    print("\n  SD correction factor: {}".format(SD_corrc_factor))
    
    # draw SD comparison plots
    SDcomparePlots(reptQuan_logFC, tmtcQuan_logFC, reptSD, tmtcSD, saveFile=os.path.join(interDir, 'log2FC_density_plots.pdf'))
    
    # correct reporter SD and log2FC
    rept_logFC_SD_corrc = rept_logFC_mean_SD['SD'] * SD_corrc_factor
    reptQuan_logFC_z = (reptQuan_logFC.sub(rept_logFC_mean_SD['Mean'].tolist(), axis='columns')).div(rept_logFC_mean_SD['SD'].tolist(), axis='columns')
    reptQuan_logFC_corrc = (reptQuan_logFC_z.mul(rept_logFC_SD_corrc.tolist(), axis=1)).add(rept_logFC_mean_SD['Mean'].tolist(),axis=1)
    # correct reporter intensity 
    reptQuan_log_corrc = reptQuan_logFC_corrc.add(rept_avg, axis=0)
    reptQuan_corrc = reptQuan_log_corrc.apply(lambda x: np.power(2,x), axis=1)
    
    # calculate a unified noise level 
    if params['unify_noise_level'] == '1':
        print("\n  Unify the noise level among channels")
        # Calculate the noise
        noise_pct = pd.DataFrame({'rept_max': reptQuan.apply(np.max, axis=1),
                                  'rept_min': reptQuan.apply(np.min, axis=1),
                                  'rept_mean': reptQuan.apply(np.mean, axis=1),
                                  'rept_max_corrc': reptQuan_corrc.apply(np.max, axis=1),
                                  'rept_min_corrc': reptQuan_corrc.apply(np.min, axis=1)})
        noise_pct['ratio_MaxMin_corrc'] = noise_pct['rept_max_corrc'] / noise_pct['rept_min_corrc']
        noise_pct['noise'] = ((noise_pct['ratio_MaxMin_corrc'] * noise_pct['rept_min']) - noise_pct['rept_max']) / (noise_pct['ratio_MaxMin_corrc'] - 1)
        
        # determin whether to use a cap for maximum noise level
        if params['use_noise_cap'] == '1':
            print("\n  Set maximum noise intensity as {}% of minimum reporter intensity".format(float(params['max_noise_pct'])*100))
            noise_pct['noise_reset'] = noise_pct['noise']
            # get index
            idx = noise_pct['noise'] > noise_pct['rept_min'] * float(params['max_noise_pct'])
            noise_pct['noise_reset'][idx] = noise_pct['rept_min'][idx] * float(params['max_noise_pct'])
            # subtract the noise intensity from original reporter intensities
            reptQuan_recorrc = reptQuan.sub(noise_pct['noise_reset'], axis='index')
            # define noise level
            noise_pct['noise_level'] = noise_pct['noise_reset'] / noise_pct['rept_mean']
            
        elif params['use_noise_cap'] == '0':
            reptQuan_recorrc = reptQuan.sub(noise_pct['noise'], axis='index')
            # define noise level
            noise_pct['noise_level'] = noise_pct['noise'] / noise_pct['rept_mean']
        
        return reptQuan_recorrc, noise_pct
    
    elif params['unify_noise_level'] == '0':
        
        return reptQuan_corrc


# 
def normTrimed_MeanSD(logfc, SDfold):
    # calculate mean and sd
    mu, std = norm.fit(logfc)
    
    # get a subset
    logfc_sub = logfc[(logfc > (mu-SDfold*std)) & (logfc < (mu+SDfold*std))]
    
    # recalculate mean and sd
    mu, std = norm.fit(logfc_sub)
    
    return mu, std


