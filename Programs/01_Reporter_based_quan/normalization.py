import numpy as np


def getSubset(df, params):
    # Get a subset of a dataframe to calculate loading-bias information
    # 1. Filter out PSMs based on the intensity level
    reporters = params["tmt_channels_used"].split(";")
    noiseLevel = 1000
    snRatio = float(params["SNratio_for_correction"])
    subDf = df[reporters][(df[reporters] > noiseLevel * snRatio).prod(axis=1).astype(bool)]  # Zero-intensity PSMs are excluded

    # 2. Filter out highly variant PSMs in each column (reporter)
    psmMean = subDf.mean(axis=1)
    subDf = np.log2(subDf.divide(psmMean, axis=0))
    pctTrimmed = float(params["percentage_trimmed"])
    n = 0
    for reporter in reporters:
        if n == 0:
            ind = ((subDf[reporter] > subDf[reporter].quantile(pctTrimmed / 200)) &
                   (subDf[reporter] < subDf[reporter].quantile(1 - pctTrimmed / 200)))
        else:
            ind = ind & ((subDf[reporter] > subDf[reporter].quantile(pctTrimmed / 200)) &
                         (subDf[reporter] < subDf[reporter].quantile(1 - pctTrimmed / 200)))
        n += 1

    subDf = subDf.loc[ind]
    return subDf


def getLoadingBias(df, params):
    ###########################
    # Loading-bias evaluation #
    ###########################
    subDf = getSubset(df, params)
    n = len(subDf)
    sm = 2 ** subDf.mean(axis=0)    # Sample-mean values
    msm = np.mean(sm)    # Mean of sample-mean values
    avg = sm / msm * 100
    sdVal = subDf.std(axis=0)
    sd = ((2 ** sdVal - 1) + (1 - 2 ** (-sdVal))) / 2 * 100
    sem = sd / np.sqrt(n)

    return avg, sd, sem, n


def normalization(df, params):
    ################################################
    # Normalization (i.e. loading-bias correction) #
    ################################################
    doNormalization = params["loading_bias_correction"]
    normalizationMethod = params["loading_bias_correction_method"]
    reporters = params["tmt_channels_used"].split(";")
    res = df.copy()

    if doNormalization == "1":
        # First, get a subset for calculating normalization factors (same as loading-bias calculation)
        # Note that this subset is 1) divided by row-wise mean (i.e. PSM-wise mean) and then 2) log2-transformed
        subDf = getSubset(df, params)

        # Calculate normalization factors for samples (reporters)
        if normalizationMethod == "1":  # Trimmed-mean
            sm = subDf.mean(axis=0)
        elif normalizationMethod == "2":  # Trimmed-median
            sm = subDf.median(axis=0)
        target = np.mean(sm)
        normFactors = sm - target

        # Normalize the input dataframe, df (in log2-scale and then scale-back)
        psmMeans = res[reporters].mean(axis=1)
        res[reporters] = np.log2(res[reporters].divide(psmMeans, axis=0).replace(0, np.nan))
        res[reporters] = res[reporters] - normFactors
        res[reporters] = 2 ** res[reporters]
        res[reporters] = res[reporters].multiply(psmMeans, axis=0)
        # After the normalization, no need to show loading-bias again (should be 100% for all samples)
        print("\n  Normalization is finished")
    else:
        print("\n  Normalization is skipped according to the parameter")

    return res
