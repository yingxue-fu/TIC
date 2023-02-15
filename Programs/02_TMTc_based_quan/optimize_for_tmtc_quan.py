# -*- coding: utf-8 -*-
"""
Created on Sat Mar 19 09:40:21 2022

@author: yfu
"""


import numpy as np
from scipy.optimize import least_squares
from get_tmtc_groups import get_TMTc_groups
from pre_defined_variables import TMTtag_whole_impurity_dict, TMTtag_balancer_impurity_dict


# Optimization for TMTc-based quantification 

# Define the difference function
# make difference between inferred and real TMT peaks abundance as the y-variable 
# and sample ratio as x-variable

def TMTc_Diff(ratio_x, used_channels, pep_iso_vector, n_TMT, TMT_ver, real_TMTc_relative_abun, positions_for_fitting):
    
    # inferred TMTc relative abundance
    inferred_TMTc_abun = infer_TMTc_peaks_abun(ratio_x, used_channels, pep_iso_vector, n_TMT, TMT_ver)  
    # 
    theoretical_TMTc_abun_norm = inferred_TMTc_abun[positions_for_fitting] / sum(inferred_TMTc_abun[positions_for_fitting])
    
    # difference function
    Diff = theoretical_TMTc_abun_norm - real_TMTc_relative_abun
    
    return Diff

# 
def optimize_for_ratios(used_channels, pep_iso_vector, n_TMT, TMT_ver, real_TMTc_relative_abun, positions_for_fitting):
    #equal_sample_ratios = np.ones(len(used_channels)) / len(used_channels)
    #
    #TMTc_group_ratio = calc_merged_sample_ratios(equal_sample_ratios, used_channels, TMT_ver)
    #initial_x = np.array(TMTc_group_ratio['merged_sample_ratio'])
    
    # intial ratio
    initial_x = np.ones(len(used_channels)) / len(used_channels)
    # lowest possible ratio
    lower_bound = initial_x.mean()*0.1
    # use least squares for optimization
    res = least_squares(TMTc_Diff, initial_x, bounds=(lower_bound, 1), ftol=1e-07, xtol=1e-06, args=(used_channels, 
                        pep_iso_vector, n_TMT, TMT_ver, real_TMTc_relative_abun, positions_for_fitting))
    
    return [res.fun, res.success, res.x]


###############################################################################
###### infer relaive abundance of TMTc peaks with arbitrary sample ratio ######
###############################################################################

# 1. get the isotopic distribution of TMTc ions for a single channel
# TMTc ions = TMT balancer region (N-terminal) + peptide + additional TMT tags (K)

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


# 2. multiply sample ratio with the isotopic distribution of single channel 
#    and combine all channels together

def infer_TMTc_peaks_abun(sample_ratios, used_channels, pep_iso_vector, n_TMT, TMT_ver, nTMTc = 18):
    
    # 1) calculate the isotopic distribution for each channel and multiply with sample ratios 
    n_channels = len(used_channels)
    TMTc_iso_dist_channels = np.zeros((n_channels, 20))
    for i, channel in zip(range(n_channels), used_channels):
        # calculate the isotopic distribution for each channel
        TMTc_iso_single = TMTc_iso_dist_single_channel(pep_iso_vector, TMTtag_whole_impurity_dict[TMT_ver][channel], TMTtag_balancer_impurity_dict[TMT_ver][channel], n_TMT)
        # multiply with sample ratios
        TMTc_iso_dist_channels[i, 0:len(TMTc_iso_single)] = TMTc_iso_single * sample_ratios[i]
    
    # 2) merge TMTc ions that are overlapped and shift positions 
    
    ## get unique TMTc groups   
    used_channel_TMTc_grp = get_TMTc_groups(used_channels, TMT_ver)
    all_TMTc_groups = np.array(list(used_channel_TMTc_grp['TMTc_grp']))
    n_TMTc_groups = len(np.unique(all_TMTc_groups))
    
    # 
    TMTc_peaks_all = np.zeros((n_TMTc_groups, 50))
    for j, TMTc_group in zip(range(n_TMTc_groups), np.unique(all_TMTc_groups)):
        # merge TMTc ions that are overlapped (sum together)
        merged_TMTc_peaks = np.sum(TMTc_iso_dist_channels[all_TMTc_groups == TMTc_group, :], axis=0)
        # shift positions
        TMTc_peaks_all[j, j:j+len(merged_TMTc_peaks)] = merged_TMTc_peaks
    
    # 3) sum the intensities of the same position
    if nTMTc <= 10:
        # keep the C(0) to C(+X) 
        TMTc_peaks_all = TMTc_peaks_all[:,1:nTMTc+1]
        TMTc_peaks_colsum = TMTc_peaks_all.sum(axis=0)
    elif nTMTc > 10:
        # keep the C(-1) to C(+X) 
        TMTc_peaks_all = TMTc_peaks_all[:,0:nTMTc]
        TMTc_peaks_colsum = TMTc_peaks_all.sum(axis=0)
    
    # 4) Normalize the total abundance to 1
    TMTc_relative_abun = TMTc_peaks_colsum / sum(TMTc_peaks_colsum)
    
    return TMTc_relative_abun


# Calculate theoretical_TMTc_peaks for equal sample mixing ratios. Use positions 
# that have more than 1% abundance for fitting 
def theoretical_abundant_TMTc_positions(used_channels, pep_iso_vector, n_TMT, TMT_ver, threshold=0.01):
    
    # use equal ratios among samples
    equal_sample_ratios = np.ones(len(used_channels)) / len(used_channels)
    
    # define abundant positions
    TMTc_abundant_positions = threshold < infer_TMTc_peaks_abun(equal_sample_ratios, used_channels, pep_iso_vector, n_TMT, TMT_ver)
    
    return TMTc_abundant_positions

