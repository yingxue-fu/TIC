#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 10 10:31:54 2023

@author: yfu
"""

import os, re
import numpy as np


# parse jump-f report (ID.txt file)
# Note that this part may need to be revised according to the Jump -f result format
def parse_jumpf_report(df):
    
    # extract information from 'Outfile' column
    df["Files"] = df["Outfile"].apply(lambda x: os.path.dirname(x).rsplit(".", 1)[0] + ".mzXML")
    df["Scan"] = df["Outfile"].apply(lambda x: os.path.basename(x).split(".")[1])
    df["Precursor.Charge"] = df["Outfile"].apply(lambda x: os.path.basename(x).split(".")[3]).astype(int)
    
    # define other key variables 
    df["key"] = df["Files"] + "_" + df["Scan"]
    df['Run'] = df["Files"].apply(lambda x: os.path.basename(x).rsplit(".", 1)[0])
    df['Precursor.Mz'] = (df['calcMH'] + (df['Precursor.Charge'] - 1) * 1.00727646662) / df['Precursor.Charge']
    df['Stripped.Sequence'] = df["Peptide"].apply(lambda x: re.sub('[@#$]','',x.split(".")[1]))
    
    # extract psm to peptide/protein relationship 
    pep2psm = df.groupby("Peptide")["key"].apply(lambda x: list(np.unique(x))).to_dict()
    prot2psm = df.groupby("Protein")["key"].apply(lambda x: list(np.unique(x))).to_dict()

    # get unique PSMs/Scans 
    ## 1) filter out comtaminants 
    print("    # of PSMs in ID.txt = {}".format(len(df)))
    print("\n  Remove PSMs that belong to contaminants...")
    df = df[~df['Protein'].str.contains('co|', regex=False)]
    print("    # of remaining PSMs = {}".format(len(df)))
    
    ## 2) remove overlapped PSMs -- remove PSM-Protein relationship
    print("\n  Remove overlapped PSMs mapping to multiple proteins...")
    df = df.drop(columns=['Protein','Ions','red','group','subgroup','unique','tryptic','pos'])
    df = df.drop_duplicates(ignore_index = True)
    print("    # of unique PSMs = {}".format(len(df)))
    
    ## 3) remove precursors with charge of 1
    df = rm_prec_z1(df)
    
    return df, pep2psm, prot2psm


def rm_prec_z1(df):
    
    # get PSM charge distribution 
    precCharge_stat = df['Precursor.Charge'].value_counts().to_frame()
    precCharge_stat.reset_index(inplace=True)
    precCharge_stat = precCharge_stat.rename({'index': 'Precursor Charge', 'Precursor.Charge': '# of Precursors'}, axis=1)
    print('\n  Precursor charge distribution (all runs)')
    print(precCharge_stat.to_string(index=False, col_space=(20,20)))
    
    # remove precursors with charge of 1
    if precCharge_stat['Precursor Charge'].min() == 1:
        print('\n  Remove precursors with charge of 1')
        df = df[df['Precursor.Charge'] != 1].reset_index(drop=True)
    
    return df