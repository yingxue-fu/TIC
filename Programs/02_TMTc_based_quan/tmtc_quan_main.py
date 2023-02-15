# -*- coding: utf-8 -*-
"""
Created on Sat Mar 19 09:40:21 2022

@author: yfu
"""

import numpy as np, pandas as pd
from pep_iso_dist import pepSeq_to_chemComp, std_aa_comp, iso_distri_neutral, iso_mass_inten_dict
from pre_defined_variables import TMTtag_whole_impurity_mean
from pep_tmt_iso_dist_ms2 import pep_iso_dist_ms2
from get_tmtc_groups import get_combined_channels, get_combined_ratios
from optimize_for_tmtc_quan import theoretical_abundant_TMTc_positions, optimize_for_ratios
from utils import progressBar
import multiprocessing as mp


# TMTc based quantification for each peptide in JUMPf output
def TMTc_based_quan(df, params):
    # input file: dfId_TMTc generated after the TMTc peaks extraction step    
    
    # read parameters
    TMT_ver = params["tmt_version"]
    used_channels = params["tmt_channels_used"].split(";")
    
    # combine channels that combine together in TMTc
    comb_channels = get_combined_channels(used_channels, TMT_ver)
  
    # for each precursor in one run #
    n_Prec = df.shape[0]

    # 
    res = pd.DataFrame(columns = comb_channels + ['Diff_square_sum','total.TMTc.intensity'])
    
    #print("\n  Performing TMTc-based quantification")
    progress = progressBar(n_Prec)
    for x in range(n_Prec):
        progress.increment()  
    
        # get peptide sequence
        pep_seq = df['Stripped.Sequence'][x]
        prec_charge = df["Precursor.Charge"][x]
        
        # Calculate natural peptide isotopic distribution 
        pep_chem_com = pepSeq_to_chemComp(pep_seq, Charge=0, aa_comp=std_aa_comp, TMT_ver='None')
        pep_iso_distri = iso_distri_neutral(iso_mass_inten_dict, pep_chem_com, isotope_cutoff=1E-6, 
                                            mass_tolerance=50, mass_calculation_method=1)
        pep_iso_vector = pep_iso_distri.isotope_inten.values
        
        # Get the number of TMT tags on the peptide
        n_TMT = pep_seq.count('K') + 1
        
        # when using a narrow MS1 isolation window, re-calculate the peptide isotopic distribution
        if float(params['isolation_width']) <= 5:
            pep_iso_vector = pep_iso_dist_ms2(pep_iso_vector, n_TMT, TMT_ver, TMTtag_whole_impurity_mean, prec_charge, params)
        
        
        # convert intensity to relative ratio for real TMTc peaks
        real_TMTc_intensity = df.loc[x, df.columns.str.contains(pat = 'sig_C')]
        
        # calculate theoretical abundant TMTc positions
        positions_for_fitting = theoretical_abundant_TMTc_positions(used_channels, pep_iso_vector, n_TMT, TMT_ver)
        real_TMTc_intensity = real_TMTc_intensity[positions_for_fitting]
        
        # calculate relative abundance of TMTc peaks
        real_TMTc_intensity_sum = sum(real_TMTc_intensity)
        real_TMTc_relative_abun = (real_TMTc_intensity / real_TMTc_intensity_sum).to_numpy(dtype=float)
             
        # optimization for minimum Diff between experimental and theoritical
        opt_result = optimize_for_ratios(used_channels, pep_iso_vector, n_TMT, TMT_ver, real_TMTc_relative_abun, positions_for_fitting)
        
        # merge ratios 
        sample_ratios = opt_result[2] / sum(opt_result[2])
        merged_ratio = get_combined_ratios(sample_ratios, used_channels, TMT_ver)
        
        # convert ratio to intensity
        temp_df = pd.DataFrame([merged_ratio.Ratio.values*real_TMTc_intensity_sum], columns = comb_channels, index = [x])
        temp_df['Diff_square_sum'] = sum(opt_result[0]**2)
        temp_df['total.TMTc.intensity'] = real_TMTc_intensity_sum
        
        res = pd.concat([res, temp_df], ignore_index=False)
    
    # drop TMTc peaks info
    df = df.loc[:, ~df.columns.str.contains(pat = 'mz_C|sig_C')]
    
    # Combine results
    res = pd.concat([df, res], axis=1)
    
    return res


# 
def parallelize_dataframe(fn, df, args):
    
    n_cores = round(mp.cpu_count() / 2)
    pool = mp.Pool(n_cores)
    
    list_of_dict = []
    df_split = np.array_split(df, n_cores)
    
    for x in df_split:
        x.reset_index(drop=True, inplace=True)
        ls1 = pool.apply_async(fn, (x, args))
        list_of_dict.append(ls1)
        
    pool.close()
    
    return list_of_dict

