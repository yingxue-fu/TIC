import os
import pandas as pd
import numpy as np


def correctImpurity(df, params):
    if params['impurity_correction'] == "1":
        reporters = params["tmt_channels_used"].split(";")
        # read impurity matrix 
        dir_path = os.path.join(os.path.dirname(os.path.realpath(__file__)), 'ImpurityTables')
        dfImpurity = pd.read_table(os.path.join(dir_path, params["impurity_matrix"]+'.ini'), sep="\t", skiprows=1, header=None, index_col=0)
        dfImpurity = pd.DataFrame(np.linalg.pinv(dfImpurity.values), dfImpurity.columns, dfImpurity.index)
        dfImpurity.columns = reporters
        dfCorrected = df[reporters].dot(dfImpurity.T)
        dfCorrected.columns = reporters
        df[reporters] = pd.concat([df[reporters]/2, dfCorrected]).groupby(level=0).max()

    return df

