#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 12 17:13:17 2023

@author: yfu
"""


import os, sys, shutil
import numpy as np, pandas as pd
from datetime import datetime
from utils import getParams, Tee
from parse_ID_report import parse_ID_report
#from tmtc_parse_jumpf_report import parse_jumpf_report
from extract_TMTc_mz_intensity import get_TMTc_mz, extract_TMTc_mz_intensity
from tmtc_quan_main import TMTc_based_quan, parallelize_dataframe
from tmtc_visualization import TMTc_quan_barplot, n_TMTc_fractions, opt_metrics_plots#, TMTc_quan_histogram
from tmtc_quan_filters import TMTc_filterPSMs
from tmtc_summarization import TMTc_summarization
from publication import generateTables


if __name__ == "__main__":

    print("\n  Initializing the program...")
    
    startTime = datetime.now()
    startTimeString = startTime.strftime("%Y/%m/%d %H:%M:%S")
    
    ##################
    # Initialization #
    ##################
    paramFile = sys.argv[1]
    params = getParams(paramFile)
    # make the output directory 
    saveDir = os.path.join(os.getcwd(), params["save_dir"], "02_TMTc_quan")
    os.makedirs(saveDir, exist_ok=True)
    # copy parameter file
    shutil.copy(paramFile, os.path.join(saveDir,'jump_q_02_tmtc.params'))

    with Tee(os.path.join(saveDir,'log.txt')):
        
        print("\n  JUMPq for TMTc-based quantification")
        print("  " + startTimeString)
        
        # make the intermediate directory 
        interDir = os.path.join(saveDir, 'intermediate')
        os.makedirs(interDir, exist_ok=True)
        # make the publication directory 
        pubDir = os.path.join(saveDir, 'publication')
        os.makedirs(pubDir, exist_ok=True)
        
        ####################
        # Read ID.txt file #
        ####################
        
        print("\n  Loading ID.txt file")
    
        # Note that this part may need to be revised according to the Identification result format
        dfId = pd.read_table(params["idtxt"], sep="\t", header=0)
        dfId, pep2psm, prot2psm = parse_ID_report(dfId)
        
        # the following two lines are specific to jump-f ID.txt
        #dfId = pd.read_table(params["idtxt"], sep=";", skiprows=1, header=0)
        #dfId, pep2psm, prot2psm = parse_jumpf_report(dfId)
        # 
        files = dfId["Files"].unique()
        
        ##########################
        # Extract TMTc ion peaks #
        ##########################
        
        # calculate theoretical TMTc m/z
        dfId = get_TMTc_mz(dfId, params, nTMTc = 18)
        
        # extract TMTc m/z and intensity
        dfId_TMTc, TMTc_summary = extract_TMTc_mz_intensity(files, dfId, params)
                
        print("\n  Number of extracted mono TMTc peaks for all PSMs")
        print(TMTc_summary['n_mono_TMTc_peaks'].to_string(index=False, col_space=(24,15)))
        
        # plot # of mono TMTc peaks and save extracted TMTc peaks
        #n_TMTc_fractions(dfId_TMTc, "mono", interDir, saveFile=os.path.join(interDir,'n_extracted_TMTc_peaks_in_each_run.pdf'))
    
        # filter by number of TMTc peaks (mono)
        dfId_TMTc = dfId_TMTc[dfId_TMTc['n_mono_TMTc_peaks'] >= int(params['min_n_TMTc_peaks'])].reset_index(drop=True)
        print("\n  Filtering PSMs based on # of mono TMTc peaks (n >= {})".format(int(params['min_n_TMTc_peaks'])))
        print("    Hereafter, {} PSMs will be used for TMTc-based quantification".format(len(dfId_TMTc)))
    
        #############################
        # TMTc-based quantification #
        #############################
        
        print("\n  Performing TMTc-based quantification")
        tmtcQuan = pd.DataFrame()
        result = parallelize_dataframe(TMTc_based_quan, dfId_TMTc, params)
        for x in result:
            tmtcQuan = pd.concat([tmtcQuan, x.get()]) 
        
        # write the report to file
        tmtcQuan.set_index('key', inplace=True)
        tmtcQuan.to_csv(os.path.join(interDir, "TMTc_quan_PSMs.txt"), sep="\t", index=False)
        
        # optimization metrics plot
        #opt_metrics_plots(tmtcQuan, saveFile=os.path.join(interDir,'Sum of squared differences of optimization.pdf'))
        
        # filtering by SSD cutoff
        tmtcQuan = tmtcQuan[tmtcQuan['Diff_square_sum'] < float(params['SSD_cutoff'])]
        
        print("\n  Filtering PSMs based on SSD of optimization (< {})".format(float(params['SSD_cutoff'])))
        print("    Hereafter, {} PSMs will be used for further analysis".format(len(tmtcQuan)))
        
        # filtering based on intensity
        tmtcQuan = TMTc_filterPSMs(tmtcQuan, params)
        #tmtcQuan.to_csv(os.path.join(interDir, "TMTc_quan_PSMs_filtered.txt"), sep="\t", index=False)
        
        # get the quantification matrix
        tmtcQuan = tmtcQuan.loc[:, tmtcQuan.columns.str.contains(pat = 'sig1')]
        
        ###############################################
        # Summarization for TMTc-based quantification #
        ###############################################
        # 1. Peptide-level summarization
        print("\n  Protein-level summarization (TMTc) is being performed")
        tmtc_dfPep = TMTc_summarization(pep2psm, tmtcQuan, params, 'peptide')
        # 2. Protein-level summarization
        print("\n  Protein-level summarization (TMTc) is being performed")
        tmtc_dfProt = TMTc_summarization(prot2psm, tmtcQuan, params, 'protein')
            
        ###################################
        # Draw barplot of relative ratios #
        ###################################
        if params['draw_ratio_plot'] == '1':
            if params['theoretical_ecoli_ratio'] != '0':
                # extract Ecoli proteins
                tmtc_dfProt_ecoli =  tmtc_dfProt[['ECOLI' in s for s in tmtc_dfProt.index]]
                theoretical_ecoli_ratios = np.array([int(i) for i in params["theoretical_ecoli_ratio"].split(",")])
                TMTc_quan_barplot(tmtc_dfProt_ecoli, 'E.coli', theoretical_ecoli_ratios, params, level='Protein', saveFile=os.path.join(interDir,'tmtc_quan_Ecoli_barplot.pdf'))
                #TMTc_quan_histogram(tmtc_dfProt_ecoli, 'E.coli', theoretical_ecoli_ratios, params, level='Protein', saveFile=os.path.join(interDir,'tmtc_quan_Ecoli_histogram.pdf'))
            if params['theoretical_human_ratio'] != '0':
                # extract Human proteins
                tmtc_dfProt_human =  tmtc_dfProt[['HUMAN' in s for s in tmtc_dfProt.index]]
                theoretical_human_ratios = np.array([int(i) for i in params["theoretical_human_ratio"].split(",")])
                TMTc_quan_barplot(tmtc_dfProt_human, 'Human', theoretical_human_ratios, params, level='Protein', saveFile=os.path.join(interDir,'tmtc_quan_Human_barplot.pdf'))
                #TMTc_quan_histogram(tmtc_dfProt_human, 'Human', theoretical_human_ratios, params, level='Protein', saveFile=os.path.join(interDir,'tmtc_quan_Human_histogram.pdf'))
        
        ######################
        # Publication tables #
        ######################
        dfUniPep_tmtc, dfAllPep_tmtc, dfUniProt_tmtc, dfAllProt_tmtc = generateTables(tmtc_dfPep, tmtc_dfProt, params)
        dfUniPep_tmtc.to_csv(os.path.join(pubDir, "id_uni_pep_quan_tmtc.txt"), sep="\t", index=False)
        dfAllPep_tmtc.to_csv(os.path.join(pubDir, "id_all_pep_quan_tmtc.txt"), sep="\t", index=False)
        dfUniProt_tmtc.to_csv(os.path.join(pubDir, "id_uni_prot_quan_tmtc.txt"), sep="\t", index=False)
        dfAllProt_tmtc.to_csv(os.path.join(pubDir, "id_all_prot_quan_tmtc.txt"), sep="\t", index=False)
        
        
        endTime = datetime.now()
        endTimeString = endTime.strftime("%Y/%m/%d %H:%M:%S")
        print("\n  " + endTimeString)
        elapsed = (endTime - startTime).total_seconds()
        print("  Finished in {} seconds \n".format(int(elapsed)))
