# -*- coding: utf-8 -*-
"""
Created on Sat Mar 19 09:40:21 2022

@author: yfu
"""

import numpy as np
from get_tmtc_groups import get_TMTc_groups
from pre_defined_variables import TMTtag_whole_impurity_dict, TMTtag_balancer_impurity_dict


# infer TMTc relaive abundance based on arbitrary sample ratio 

# 1. infer the relative abundance of a single channel
def TMTc_iso_dist_single_channel(pep_iso_vector, whole_TMTtag_impurity_vec, balan_TMTtag_impurity_vec, n_TMT):
    
    # 1) convolve the impurity of all TMT(pro) tags in TMTc ions
    integrated_TMTtags_impurity = np.array([1])
    if n_TMT == 1: # when the peptide is only labeled with one TMT(pro) tag
        integrated_TMTtags_impurity = balan_TMTtag_impurity_vec    
    else:
        # When n_TMT >= 2, first convolve the impurity of whole TMT(pro) tags
        for x in range(n_TMT-1):
            integrated_TMTtags_impurity = np.convolve(integrated_TMTtags_impurity, whole_TMTtag_impurity_vec)
        # Then convolve the above isotopic evelope with one TMT(pro) balancer tag
        integrated_TMTtags_impurity = np.convolve(integrated_TMTtags_impurity, balan_TMTtag_impurity_vec)
    
    # 2) convolve peptide isotope distribution with the integrated TMT(pro) tag impurity
    TMTc_single_channel = np.convolve(pep_iso_vector, integrated_TMTtags_impurity)
    
    # 3) keep position P(-1) ~ P(+10)
    p_neg_one = n_TMT-1 # position of P(-1)
    TMTc_single_channel_final = TMTc_single_channel[p_neg_one:p_neg_one+12]

    return TMTc_single_channel_final


# 2. multiply sample ratio with the relaive abundance of single channel to form TMTc peaks
def infer_TMTc_peaks_abun(sample_ratios, used_channels, pep_iso_vector, n_TMT, TMT_ver, nTMTc):
    

    # 1) get unique TMTc groups   
    n_channels = len(used_channels)
    used_channel_TMTc_group = get_TMTc_groups(used_channels, TMT_ver)
    all_TMTc_groups = np.array(list(used_channel_TMTc_group['TMTc_group']))
    uni_TMTc_groups = np.unique(all_TMTc_groups)
    n_TMTc_groups = len(uni_TMTc_groups)
    
    # 2) calculate the isotopic distribution for each channel 
    TMTc_iso_dist_channels = np.zeros((n_channels, 100))
    for i, channel in zip(range(n_channels), used_channels):
        TMTc_iso_single = TMTc_iso_dist_single_channel(pep_iso_vector, TMTtag_whole_impurity_dict[TMT_ver][channel], TMTtag_balancer_impurity_dict[TMT_ver][channel], n_TMT)
        #TMTc_iso_dist_channels[i, 0:len(TMTc_iso_single)] = TMTc_iso_single 
        TMTc_iso_dist_channels[i, 0:len(TMTc_iso_single)] = TMTc_iso_single * sample_ratios[i]
    
    # 3) multiply with sample ratios and combine signal from all channels
    TMTc_peaks_all = np.zeros((n_TMTc_groups, 150))
    for j, TMTc_group in zip(range(n_TMTc_groups), uni_TMTc_groups):
        #merged_TMTc_peaks = np.mean(TMTc_iso_dist_channels[all_TMTc_groups == TMTc_group, :], axis=0)
        merged_TMTc_peaks = np.sum(TMTc_iso_dist_channels[all_TMTc_groups == TMTc_group, :], axis=0)
        #TMTc_peaks_all[j, j:j+len(merged_TMTc_peaks)] = merged_sample_ratios_num[j] * merged_TMTc_peaks 
        TMTc_peaks_all[j, j:j+len(merged_TMTc_peaks)] = merged_TMTc_peaks 
    
    
    if nTMTc <= 10:
        # 4) keep the C(0) to C(+X) 
        TMTc_peaks_all = TMTc_peaks_all[:,1:nTMTc]
        TMTc_peaks_colsum = TMTc_peaks_all.sum(axis=0)
    elif nTMTc > 10:    
        # 4) keep the C(-1) to C(+X) 
        TMTc_peaks_all = TMTc_peaks_all[:,0:nTMTc]
        TMTc_peaks_colsum = TMTc_peaks_all.sum(axis=0)
    
    # 5) Normalize the total abundance to 1
    TMTc_relative_abun = TMTc_peaks_colsum / sum(TMTc_peaks_colsum)
    
    return TMTc_relative_abun


# Calculate theoretical_TMTc_peaks for equal sample mixing ratios. Use positions 
# that have more than 1% abundance for fitting 
def theoretical_abundant_TMTc_positions(used_channels, pep_iso_vector, n_TMT, TMT_ver, threshold=0.01):
    #
    equal_sample_ratios = np.ones(len(used_channels)) / len(used_channels)
    #
    #merged_sample_ratios = calc_merged_sample_ratios(equal_sample_ratios, used_channels, TMT_ver)
    #merged_sample_ratios_num = np.array(merged_sample_ratios['merged_sample_ratio'])
    #
    #TMTc_abundant_positions = threshold < infer_TMTc_peaks_abun(merged_sample_ratios_num, used_channels, pep_iso_vector, n_TMT, TMT_ver)
    TMTc_abundant_positions = threshold < infer_TMTc_peaks_abun(equal_sample_ratios, used_channels, pep_iso_vector, n_TMT, TMT_ver)
    return TMTc_abundant_positions