# -*- coding: utf-8 -*-

import numpy as np

# 1. calculate the isotopic composition of a TMT(pro)-labeled peptide

def calc_pep_tmt_iso_compos(pep_iso_vector, n_TMT, TMT_ver, TMTtag_whole_impurity_mean):
    # pep_iso_vector: the natural isotopic distribution of a peptide calculated based on its sequence
    # TMT_ver: 'TMT' or 'TMTpro'

    # 1) when the number of TMT(pro) tags on peptide > 1, first convolve the impurity of all TMT(pro) tags 
    integrated_TMTtags_impurity = np.array([1])
    for x in range(n_TMT):
        integrated_TMTtags_impurity = np.convolve(integrated_TMTtags_impurity, TMTtag_whole_impurity_mean[TMT_ver])
    
    # 2) convolve natural peptide isotope distribution with the integrated TMT tag impurity
    pep_tmt_iso_composition_temp = np.outer(pep_iso_vector, integrated_TMTtags_impurity)
    ## shift position based on peptide isotopes
    pep_tmt_iso_composition = np.zeros((pep_tmt_iso_composition_temp.shape[0], 100))
    for x in range(pep_tmt_iso_composition_temp.shape[0]):
        pep_tmt_iso_composition[x,x:x+pep_tmt_iso_composition_temp.shape[1]] = pep_tmt_iso_composition_temp[x]

    # 3) determine M+0 position 
    zero_position = pep_tmt_iso_composition[0].argmax()
    ## keep the M-1 to M+10
    pep_tmt_iso_composition = pep_tmt_iso_composition[:,(zero_position-1):(zero_position+11)]

    return pep_tmt_iso_composition


# 2. recalculate peptide isotopic distribution based on isolation window width

def pep_iso_dist_ms2(pep_iso_vector, n_TMT, TMT_ver, TMTtag_whole_impurity_mean, prec_charge, params):
    
    pep_tmt_iso_composition = calc_pep_tmt_iso_compos(pep_iso_vector, n_TMT, TMT_ver, TMTtag_whole_impurity_mean)
    
    # m/z difference between two isotopic peaks of a precursor 
    delta_mz = 1 / prec_charge 
    
    # calculate the relative m/z of all precursor peaks M-1 ~ M+10, with monoisotopic peak as M+0
    relative_mz = np.append(np.array([0-delta_mz, 0]), delta_mz*np.arange(1,11,1))
    
    # calculate the isolation window 
    isolation_width = float(params['isolation_width'])
    isolation_offset = float(params['isolation_offset'])
    right_bound = 0 + (isolation_width/2) + isolation_offset
    left_bound = 0 - (isolation_width/2) + isolation_offset

    # extract precursor peaks inside the window
    isolated_prec_peaks = pep_tmt_iso_composition[:,(relative_mz > left_bound) & (relative_mz < right_bound)]
    
    # renormalize the sum to 1
    isolated_prec_peaks_norm = isolated_prec_peaks / isolated_prec_peaks.sum()
    pep_iso_vector_ms2 = isolated_prec_peaks_norm.sum(axis = 1)

    return pep_iso_vector_ms2
