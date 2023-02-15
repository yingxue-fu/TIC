#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 12 17:19:17 2023

@author: yfu
"""

import os, sys, shutil
import numpy as np, pandas as pd
from datetime import datetime
from utils import getParams, Tee
from rept_z_correction import rept_z_correct
from rept_corrc_publication import generate_protTables
import matplotlib.pyplot as plt
from visualization import quan_mtx_barplot#, quan_mtx_histogram

if __name__ == "__main__":
    
    print("\n  Initializing the program...")
    
    startTime = datetime.now()
    startTimeString = startTime.strftime("%Y/%m/%d %H:%M:%S")
    
    ##################
    # Initialization #
    ##################
    paramFile = sys.argv[1]
    params = getParams(paramFile)
    saveDir = os.path.join(os.getcwd(), params["save_dir"], "03_reporter_quan_corrc")
    os.makedirs(saveDir, exist_ok=True)
    # copy parameter file
    shutil.copy(paramFile, os.path.join(saveDir,'jump_q_03_rept_corrc.params'))

    with Tee(os.path.join(saveDir,'log.txt')):
        
        print("\n  JUMPq for reporter correction based on TMTc")
        print("  " + startTimeString)
        
        # make the intermediate directory 
        interDir = os.path.join(saveDir,'intermediate')
        os.makedirs(interDir, exist_ok=True)
        # make the publication directory 
        pubDir = os.path.join(saveDir,'publication')
        os.makedirs(pubDir, exist_ok=True)
        
        #################################################
        # Read reporter and TMTc quantification results #
        #################################################
        ## reporter
        print("\n  Reading reporter-based quantification data")
        print("    from {}".format(params['reptQuan']))
        reptQuan_df = pd.read_table(params['reptQuan'], sep="\t")
        reptQuan = reptQuan_df.loc[:, reptQuan_df.columns.str.contains(pat = 'sig1')]
        reptQuan.set_index(reptQuan_df['Protein Accession #'], inplace=True)
        
        ## TMTc 
        print("\n  Reading TMTc-based quantification data")
        print("    from {}".format(params['tmtcQuan']))
        tmtcQuan_df = pd.read_table(params['tmtcQuan'], sep="\t")
        tmtcQuan = tmtcQuan_df.loc[:, tmtcQuan_df.columns.str.contains(pat = 'sig1')]
        
        ###############################
        # TMTc-based noise correction #
        ###############################
        if params['unify_noise_level'] == '1': # unify noise level among channels
            rept_dfProt_corrc, noise_pct = rept_z_correct(reptQuan, tmtcQuan, interDir, params)
            
            #noise_pct.to_csv(os.path.join(interDir, "noise_percentage.txt"), sep="\t", index=True)
            # draw plot of noise percentage
            #ax = noise_pct.hist(column='noise_level', bins=30, grid=False)[0][0]
            #ax.set_xlabel('Noise intensity / Mean reporter intensity', fontsize=14)
            #ax.set_ylabel('# of proteins', fontsize=14)
            #ax.set_title('Noise percentage', fontsize=16)
            #ax.tick_params(axis='both', which='major', labelsize=12)
            #plt.tight_layout()
            #plt.savefig(os.path.join(interDir, 'noise_percentage_histogram.pdf'))
            
        elif params['unify_noise_level'] == '0': # don't unify noise level among channels
            rept_dfProt_corrc = rept_z_correct(reptQuan, tmtcQuan, interDir, params)
        
        ###################################
        # Draw barplot of relative ratios #
        ###################################
        if params['draw_ratio_plot'] == '1':
            if params['theoretical_ecoli_ratio'] != '0':
                # extract Ecoli proteins
                dfProt_ecoli =  rept_dfProt_corrc[['ECOLI' in s for s in rept_dfProt_corrc.index]]
                theoretical_ecoli_ratios = np.array([int(i) for i in params["theoretical_ecoli_ratio"].split(",")])
                quan_mtx_barplot(dfProt_ecoli, 'E.coli', theoretical_ecoli_ratios, params, level='Protein', saveFile=os.path.join(interDir,'reporter_quan_corrc_Ecoli_barplot.pdf'))
                #quan_mtx_histogram(dfProt_ecoli, 'E.coli', theoretical_ecoli_ratios, params, level='Protein', saveFile=os.path.join(interDir,'reporter_quan_corrc_Ecoli_histogram.pdf'))
            if params['theoretical_human_ratio'] != '0':
                # extract Human proteins
                dfProt_human =  rept_dfProt_corrc[['HUMAN' in s for s in rept_dfProt_corrc.index]]
                theoretical_human_ratios = np.array([int(i) for i in params["theoretical_human_ratio"].split(",")])
                quan_mtx_barplot(dfProt_human, 'Human', theoretical_human_ratios, params, level='Protein', saveFile=os.path.join(interDir,'reporter_quan_corrc_Human_barplot.pdf'))
                #quan_mtx_histogram(dfProt_human, 'Human', theoretical_human_ratios, params, level='Protein', saveFile=os.path.join(interDir,'reporter_quan_corrc_Human_histogram.pdf'))
            
        ###############
        # Publication #
        ###############
        print("\n  Generating publication tables")
        dfUniProt_corrc, dfAllProt_corrc = generate_protTables(rept_dfProt_corrc, params)
        dfUniProt_corrc.to_csv(os.path.join(pubDir, "id_uni_prot_quan_corrc.txt"), sep="\t", index=False)
        dfAllProt_corrc.to_csv(os.path.join(pubDir, "id_all_prot_quan_corrc.txt"), sep="\t", index=False)
        
        
        endTime = datetime.now()
        endTimeString = endTime.strftime("%Y/%m/%d %H:%M:%S")
        print("\n  " + endTimeString)
        elapsed = (endTime - startTime).total_seconds()
        print("  Finished in {} seconds \n".format(int(elapsed)))