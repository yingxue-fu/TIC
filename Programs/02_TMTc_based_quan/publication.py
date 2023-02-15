import os
import pandas as pd


def generateTables(dfPep, dfProt, params):
    dirPub = os.path.dirname(params['idtxt']) + "/publications/"

    ########################
    # Peptide-level tables #
    ########################
    dfPep = dfPep.reset_index()
    dfPep = dfPep.rename({"index": "Peptides"}, axis=1)

    # Unique peptides
    dfJumpf = pd.read_table(os.path.join(dirPub, "id_uni_pep.txt"), sep="\t", skiprows=3, header=0)
    dfUniPep = dfJumpf.merge(dfPep, on="Peptides")

    # All peptides
    dfJumpf = pd.read_table(os.path.join(dirPub, "id_all_pep.txt"), sep="\t", skiprows=3, header=0)
    dfAllPep = dfJumpf.merge(dfPep, on="Peptides")

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

    return dfUniPep, dfAllPep, dfUniProt, dfAllProt