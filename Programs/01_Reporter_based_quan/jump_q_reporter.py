import os, sys, shutil
import numpy as np
import pandas as pd
from datetime import datetime
from utils import getParams, Tee
from reporter import extractReporters
from impurity import correctImpurity
from filters import filterPSMs
from normalization import getLoadingBias, normalization
from visualization import quan_mtx_barplot#, quan_mtx_histogram
from summarization import summarization
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
    saveDir = os.path.join(os.getcwd(), params["save_dir"], "01_reporter_quan")
    os.makedirs(saveDir, exist_ok=True)
    # copy parameter file
    shutil.copy(paramFile, os.path.join(saveDir,'jump_q_01_rept.params'))
    
    with Tee(os.path.join(saveDir,'log.txt')):

        print("\n  JUMPq for TMT reporter-based quantification")
        print("  " + startTimeString)
        
        # make the intermediate directory 
        interDir = os.path.join(saveDir, 'intermediate')
        os.makedirs(interDir, exist_ok=True)
        # make the publication directory 
        pubDir = os.path.join(saveDir, 'publication')
        os.makedirs(pubDir, exist_ok=True)
        
        ##################
        # Parsing ID.txt #
        ##################
    
        print("\n  Loading ID.txt file")
    
        # Note that this part may need to be revised according to the identification result format
        # the following three lines are specific to jump-f ID.txt
        #dfId = pd.read_table(params["idtxt"], sep=";", skiprows=1, header=0)
        #dfId["Files"] = dfId["Outfile"].apply(lambda x: os.path.dirname(x).rsplit(".", 1)[0] + ".mzXML")
        #dfId["Scan"] = dfId["Outfile"].apply(lambda x: os.path.basename(x).split(".")[1])
        
        dfId = pd.read_table(params["idtxt"], sep="\t", header=0)

        # define a key for each PSM
        dfId["key"] = dfId["Files"] + "_" + dfId["Scan"].astype(str)
        files = dfId["Files"].unique()
    
        ##################################
        # Extract TMT reporter ion peaks #
        ##################################
        # 1st round of reporter ion extraction
        dfQuan, reporterSummary = extractReporters(files, dfId, params)
    
        # Before 2nd round of TMT reporter extraction, m/z-shifts of reporters are summarized
        print("\n  m/z-shift in each TMT reporter")
        reporters = params["tmt_channels_used"].split(";")
        for reporter in reporters:
            m = reporterSummary[reporter]["meanMzShift"]
            s = reporterSummary[reporter]["sdMzShift"]
            print("    %s\tm/z-shift = %.4f [ppm]\tsd = %.4f" % (reporter, m, s))
    
        # 2nd round of reporter ion extraction
        dfQuan, reporterSummary = extractReporters(files, dfId, params, **reporterSummary)
    
        ###########################
        # TMT impurity correction #
        ###########################
        dfQuan = correctImpurity(dfQuan, params)
    
        #####################
        # Filtering of PSMs #
        #####################
        dfQuan = filterPSMs(dfQuan, params)
    
        #####################################
        # Show the loading-bias information #
        #####################################
        avgLb, sdLb, semLb, nn = getLoadingBias(dfQuan, params)
        print("\n  Loading bias (before correction)")
        print("    Reporter\tMean[%]\tSD[%]\tSEM[%]\t#PSMs")
        for i in range(len(reporters)):
            print("    %s\t%.2f\t%.2f\t%.2f\t%d" % (reporters[i], avgLb[i], sdLb[i], semLb[i], nn))
    
        #################
        # Normalization #
        #################
        dfNorm = normalization(dfQuan, params)
        dfNorm.to_csv(os.path.join(interDir, "normalized_quan_psm_nonzero.txt"), sep="\t")
        
        #################
        # Summarization #
        #################
        # 1. Peptide-level summarization
        print("\n  Peptide-level summarization is being performed")
        pep2psm = dfId.groupby("Peptide")["key"].apply(lambda x: list(np.unique(x))).to_dict()
        dfPep = summarization(pep2psm, dfNorm, params, 'peptide')
        # 2. Protein-level summarization
        print("\n  Protein-level summarization is being performed")
        prot2psm = dfId.groupby("Protein")["key"].apply(lambda x: list(np.unique(x))).to_dict()
        dfProt = summarization(prot2psm, dfNorm, params, 'protein')
        
        ###################################
        # Draw barplot of relative ratios #
        ###################################
        if params['draw_ratio_plot'] == '1':
            if params['theoretical_ecoli_ratio'] != '0':
                # extract Ecoli proteins
                dfProt_ecoli =  dfProt[['ECOLI' in s for s in dfProt.index]]
                theoretical_ecoli_ratios = np.array([int(i) for i in params["theoretical_ecoli_ratio"].split(",")])
                quan_mtx_barplot(dfProt_ecoli, 'E.coli', theoretical_ecoli_ratios, params, level='Protein', saveFile=os.path.join(interDir,'reporter_quan_Ecoli_barplot.pdf'))
                #quan_mtx_histogram(dfProt_ecoli, 'E.coli', theoretical_ecoli_ratios, params, level='Protein', saveFile=os.path.join(interDir,'reporter_quan_Ecoli_histogram.pdf'))
            if params['theoretical_human_ratio'] != '0':
                # extract Human proteins
                dfProt_human =  dfProt[['HUMAN' in s for s in dfProt.index]]
                theoretical_human_ratios = np.array([int(i) for i in params["theoretical_human_ratio"].split(",")])
                quan_mtx_barplot(dfProt_human, 'Human', theoretical_human_ratios, params, level='Protein', saveFile=os.path.join(interDir,'reporter_quan_Human_barplot.pdf'))
                #quan_mtx_histogram(dfProt_human, 'Human', theoretical_human_ratios, params, level='Protein', saveFile=os.path.join(interDir,'reporter_quan_Human_histogram.pdf'))
            
        ######################
        # Publication tables #
        ######################
        dfUniPep, dfAllPep, dfUniProt, dfAllProt = generateTables(dfPep, dfProt, params)
        dfUniPep.to_csv(os.path.join(pubDir, "id_uni_pep_quan.txt"), sep="\t", index=False)
        dfAllPep.to_csv(os.path.join(pubDir, "id_all_pep_quan.txt"), sep="\t", index=False)
        dfUniProt.to_csv(os.path.join(pubDir, "id_uni_prot_quan.txt"), sep="\t", index=False)
        dfAllProt.to_csv(os.path.join(pubDir, "id_all_prot_quan.txt"), sep="\t", index=False)
    
        endTime = datetime.now()
        endTimeString = endTime.strftime("%Y/%m/%d %H:%M:%S")
        print("\n  " + endTimeString)
        elapsed = (endTime - startTime).total_seconds()
        print("  Finished in {} seconds \n".format(int(elapsed)))
    
    