#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan 20 17:18:25 2023

@author: yfu
"""

import os
import pandas as pd

def generate_protTables(dfProt, params):
    dirPub = os.path.dirname(params['idtxt']) + "/publications/"

    ########################
    # Protein-level tables #
    ########################
    dfProt = dfProt.reset_index()
    dfProt = dfProt.rename({"index": "Protein Accession #"}, axis=1)

    # Unique proteins
    dfJumpf = pd.read_table(os.path.join(dirPub, "id_uni_prot.txt"), sep="\t", skiprows=1, header=0)
    dfUniProt = dfJumpf.merge(dfProt, on="Protein Accession #")

    # All peptides
    dfJumpf = pd.read_table(os.path.join(dirPub, "id_all_prot.txt"), sep="\t", skiprows=1, header=0)
    dfAllProt = dfJumpf.merge(dfProt, on="Protein Accession #")

    return dfUniProt, dfAllProt