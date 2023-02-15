import numpy as np
import pandas as pd
from utils import progressBar
from filters import filterByIntensity
from outliers import outlierRemoval


def summarization_1_2(df, params):
    # Summarization of proteins (and peptides?) mapped by 1 or 2 PSMs
    methods = params["min_intensity_method_1_2_psm"].split(",")
    thresholds = params["min_intensity_value_1_2_psm"].split(",")
    reporters = params["tmt_channels_used"].split(";")
    idx = []
    if len(df) == 1:
        # Proteins mapped by only one PSM
        # If the PSM is filtered by the intensity-based filter, it will not be used for the quantification
        idx = filterByIntensity(df, methods, thresholds, reporters, 0)
        if len(idx) > 0:
            df.drop(idx, inplace=True)
    elif len(df) == 2:
        # Proteins mapped by two PSMs
        # Apply the intensity-based filter first
        # - If both PSMs are filtered out, they will not be used for the quantification
        # - If one of PSMs is filtered out, the PSM will not be used for the quantification
        # - If none of PSMs is filtered out, go to the next step (two PSMs can be used for the quantification)
        #   - For each PSM, check the variation (i.e., stdev) across the reporters (in log2-space)
        #   - One with smaller variation will be used for the quantification
        #   - If both PSMs have the same variation, the one with higher mean intensity will be used
        idx = filterByIntensity(df, methods, thresholds, reporters, 0)
        if len(idx) > 0:
            df.drop(idx, inplace=True)
        elif len(idx) == 0:
            psmStd = np.log2(df[reporters]).std(axis=1)
            psmMean = np.log2(df[reporters]).mean(axis=1)
            if psmStd[0] == psmStd[1]:
                ii = np.argmin(psmMean)
            else:
                ii = np.argmax(psmStd)
            df.drop(df.index[ii], inplace=True)

    return df


def summarization(inputDict, df, params, level):
    # Input arguments
    # df: a dataframe containing PSM-level quantification information
    # inputDict: a dictionary containing the relationship between protein (or peptide) and PSMs
    #            e.g., prot2Ppsm: key = each protein, value = list of PSMs corresponding to the protein
    # params: parameters from the .param file
    # level: protein or peptide (using a string)

    resDict = {}
    nEntries = len(inputDict)
    nRemoved1, nRemoved2 = 0, 0
    reporters = params["tmt_channels_used"].split(";")
    progress = progressBar(nEntries)
    for entry, psms in inputDict.items():
        progress.increment()
        psms = df.index.join(psms, how="inner")
        if len(psms) == 0:
            nRemoved1 += 1
            continue
        else:
            subDf = df.loc[psms][reporters]

            # Limit the number of PSMs when there are TOO many PSMs(?)
            # Choose top-n PSMs according to the PSM-wise total intensity
            threshold = 100
            if len(subDf) > threshold:
                psmSum = subDf.sum(axis=1)
                topIndex = psmSum.sort_values(ascending=False).index[:threshold]
                subDf = subDf.loc[topIndex]

            # Summarization (both peptide and protein)
            if 0 < len(subDf) <= 2:
                subDf = summarization_1_2(subDf, params)
            elif len(subDf) >= 3:
                # Preprocessing for outlier removal
                # 1. Log2-transformation
                # 2. PSM-wise mean calculation
                # 2.1. Representative protein abundance by the mean of top3 PSM-wise means
                #      (equivalent to the grand mean of top3 PSMs)
                # 3. Mean-centering (using the PSM-wise mean obtained at step2)
                # 4. Outlier removal (using either Q-test or ESD test)
                # 5. Scale-back to the raw-scale
                subDf = np.log2(subDf)
                psmMeans = subDf.mean(axis=1)
                repAbundance = np.mean(sorted(psmMeans, reverse=True)[0:3])
                subDf = subDf.sub(psmMeans, axis=0)
                subDf = outlierRemoval(subDf, 0.05)  # Can I make it faster?
                if len(subDf) > 0:
                    subDf = 2 ** (subDf.mean(axis=0).to_frame().T + repAbundance)

            if len(subDf) > 0:
                resDict[entry] = subDf.iloc[0].to_dict()
            else:
                nRemoved2 += 1

    print("    {} (out of {}) {}s are quantified".format(nEntries - nRemoved1 - nRemoved2, nEntries, level))
    print("    {} (out of {}) {}s are NOT quantified due to the filtering of PSMs".format(nRemoved1 + nRemoved2, nEntries, level))
    # print("    {} (out of {}) {}s are NOT quantified due to the intensity-based filtering of PSMs".format(nRemoved1, nEntries, level))
    # print("    {} (out of {}) {}s are NOT quantified due to the further filtering of PSMs".format(nRemoved2, nEntries, level))
    res = pd.DataFrame.from_dict(resDict, orient="index")
    return res
