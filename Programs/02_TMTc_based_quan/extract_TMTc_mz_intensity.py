#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Dec 19 16:30:30 2022

@author: yfu
"""

import os, sys
import numpy as np, pandas as pd
from utils import progressBar
from pyteomics import mzxml
from pre_defined_variables import TMTtag_reportor_mass, TMTtag_CO_mass


#
def extract_TMTc_mz_intensity(files, df, params):
    # Input arguments
    # files: mzXML files 
    # df: dataframe of ID.txt file
    # params: parameters
    
    print("\n  Extracting TMTc m/z and intensity for identified PSMs...")
    
    dictTMTc = {}
    # for each mzXML file
    for file in files:
        print("\n    Working on {}".format(os.path.basename(file)))
        
        df_sub = df[df['Files'] == file]
        
        # read mzXML file
        ext = os.path.splitext(file)[-1]
        if ext == ".mzXML": 
            reader = mzxml.MzXML(file) 
            
            scans = list(df_sub['Scan'].unique())
            # Extraction of TMTc ions in each scan
            progress = progressBar(len(scans))
            for scan in scans:
                progress.increment()
                spec = reader[str(scan)]
                df_scan = df_sub[df_sub['Scan'] == scan].reset_index(drop=True)
                res = get_TMTc_intensity(df_scan, spec, mass_tolerance=20)
                key = file + "_" + str(scan)
                dictTMTc[key] = res
            
        else:
            sys.exit("   Only .mzXML file is supported.")
    
    # Combine the results
    res = pd.concat(dictTMTc, axis=0, ignore_index=True)
    
    # summary of extracted TMTc peaks
    TMTc_summary = {}
    TMTc_summary['n_mono_TMTc_peaks'] = res['n_mono_TMTc_peaks'].value_counts().sort_index(ascending=False).to_frame(name='n_PSMs').reset_index()
    TMTc_summary['n_mono_TMTc_peaks'].rename(columns = {'index':'# of mono TMTc peaks', 'n_PSMs':'# of PSMs'}, inplace = True)
    
    return res, TMTc_summary


# Calculate theoretical TMTc ion m/z for a precursor 
# Inputs: 
    # precursor m/z and charge;
    # used TMT channels;
    # TMT version(TMT-11 or TMT-18).

def get_TMTc_mz(dfId, params, nTMTc = 18):
    
    # No. of identified precursors
    n_Prec = dfId.shape[0]
    
    # read paramters
    TMT_ver = params["tmt_version"]
    used_channels = params["tmt_channels_used"].split(";")
    
    # create a new dataframe
    dfId['theoretical.TMTc.mz'] = ""
    
    print("\n  Caculating theoretical TMTc m/z for identified PSMs...")
    progress = progressBar(n_Prec)
    for x in range(n_Prec):
        progress.increment()
        
        # required precursor information
        precMz = dfId['Precursor.Mz'][x]
        precCharge = dfId['Precursor.Charge'][x]
        
        # calculate theoretical TMTc ion mass for each channel
        theoretical_TMTc_mass = np.zeros(len(used_channels))
        for y in range(len(used_channels)):
            # TMTc ion mass = precursor ion mass - reporter ion mass - CO mass
            theoretical_TMTc_mass[y] = precMz*precCharge - TMTtag_reportor_mass[used_channels[y]] - TMTtag_CO_mass[TMT_ver][used_channels[y]] 
        
        # calculate the TMTc ions m/z 
        theoretical_TMTc_mz = theoretical_TMTc_mass / (precCharge-1)
            
        # combine peaks that are at the same position
        TMTc_positions = np.around(theoretical_TMTc_mz - theoretical_TMTc_mz.min())
        merged_TMTc_mz = pd.DataFrame({'TMTc.position': list(TMTc_positions.astype(str)),
                                          'mz': list(theoretical_TMTc_mz)})    
        merged_TMTc_mz = merged_TMTc_mz.groupby('TMTc.position', as_index=False).agg(mzMean=('mz', 'mean'))
        mono_TMTc_mz = np.array(merged_TMTc_mz['mzMean'])
        
        #
        mz_diff = merged_TMTc_mz['mzMean'].diff().mean()
        if nTMTc == len(mono_TMTc_mz):
            dfId['theoretical.TMTc.mz'][x] = mono_TMTc_mz
        elif nTMTc == len(mono_TMTc_mz)+1:
            # add only one additional isotope peak on the right
            mono_added_TMTc_mz = np.append(mono_TMTc_mz, mono_TMTc_mz.max()+mz_diff)      
            dfId['theoretical.TMTc.mz'][x] = mono_added_TMTc_mz
        elif nTMTc >= len(mono_TMTc_mz)+2:
            # add extra TMTc positions: from C(-1) to C(X) 
            n_added_TMTc = nTMTc - len(mono_TMTc_mz)
            added_TMTc_mz = np.append(mono_TMTc_mz.max()+mz_diff*np.arange(1,n_added_TMTc,1), mono_TMTc_mz.min()-mz_diff)
            mono_added_TMTc_mz = np.append(mono_TMTc_mz, added_TMTc_mz)
            mono_added_TMTc_mz.sort()
            dfId['theoretical.TMTc.mz'][x] = mono_added_TMTc_mz
        
    return dfId


# Extract real TMTc peaks' intensities 
# extract peak intensities from MS2 spectrum based on input m/z array
def get_TMTc_intensity(df_scan, spec, mass_tolerance=20): # 20 ppm
    
    # No. of precursors
    n_prec = df_scan.shape[0]
    nTMTc = len(df_scan['theoretical.TMTc.mz'][0])
    # create a dataframe to store m/z and intensities 
    res = pd.DataFrame(index=range(n_prec), columns=range(nTMTc*2 + 2))
    # add column names
    if nTMTc <= 10:
        TMTc_positions = ["C+"+str(num) for num in range(0, nTMTc)]
        res.columns = ['mz_' + s for s in TMTc_positions] + ['sig_' + s for s in TMTc_positions] + ['n_TMTc_peaks', 'n_mono_TMTc_peaks']
    if nTMTc > 10:
        TMTc_positions = ["C-1"] + ["C+"+str(num) for num in range(0, nTMTc-1)]
        res.columns = ['mz_' + s for s in TMTc_positions] + ['sig_' + s for s in TMTc_positions] + ['n_TMTc_peaks', 'n_mono_TMTc_peaks']
    # for each identified precursor in a MS2 scan
    for i in range(n_prec):
        input_mz = df_scan['theoretical.TMTc.mz'][i]
        # for each position of TMTc
        for x in range(len(input_mz)):
            # calculat ppm 
            ppm_all_peaks = ((abs(input_mz[x] - spec['m/z array'])) * 1E6) / input_mz[x]
            # match the peak
            peak_index = list(ppm_all_peaks < mass_tolerance)
            # extraction
            if sum(peak_index) == 0: # when there is no peak
                res.iloc[i, x] = np.nan # m/z
                res.iloc[i, x+nTMTc] = 0     # intensity
            elif sum(peak_index) == 1: # when there is one peak
                res.iloc[i, x] = spec['m/z array'][peak_index][0]
                res.iloc[i, x+nTMTc] = spec['intensity array'][peak_index][0]
            else: # when there are multiple peaks
                # calculate the diff with the input_mz and choose peak with the smallest diff
                mz_diff = abs(spec['m/z array'][peak_index] - input_mz[x])
                res.iloc[i, x] = spec['m/z array'][peak_index][np.argmin(mz_diff)]
                res.iloc[i, x+nTMTc] = spec['intensity array'][peak_index][np.argmin(mz_diff)]
        
        res['n_TMTc_peaks'][i] = sum(~(res.loc[i, res.columns.str.contains("sig_")] == 0))
        res['n_mono_TMTc_peaks'][i] = sum(~(res.loc[i, ['sig_C+0','sig_C+1','sig_C+2','sig_C+3','sig_C+4','sig_C+5','sig_C+6','sig_C+7','sig_C+8']] == 0))
        
    return pd.concat([df_scan, res], axis=1)
