import numpy as np, pandas as pd


# Get a subset of a dataframe to calculate loading-bias information
def TMTc_getSubset(df):
    
    # Filter out highly variant PSMs in each column (reporter)
    rowMean = df.mean(axis=1)
    df_log = np.log2(df.divide(rowMean, axis=0))
    pctTrimmed = 20
    n = 0
    for i in range(df_log.shape[1]):
        if n == 0:
            ind = ((df_log.iloc[:, i] > df_log.iloc[:, i].quantile(pctTrimmed / 200)) &
                   (df_log.iloc[:, i] < df_log.iloc[:, i].quantile(1 - pctTrimmed / 200)))
        else:
            ind = ind & ((df_log.iloc[:, i] > df_log.iloc[:, i].quantile(pctTrimmed / 200)) &
                         (df_log.iloc[:, i] < df_log.iloc[:, i].quantile(1 - pctTrimmed / 200)))
        n += 1

    df_log = df_log.loc[ind]
    
    return df_log



def TMTc_getLoadingBias(df):
    ###########################
    # Loading-bias evaluation #
    ###########################
    subDf = TMTc_getSubset(df)
    n = len(subDf)
    sm = 2 ** subDf.mean(axis=0)    # Sample-mean values
    msm = np.mean(sm)    # Mean of sample-mean values
    avg = sm / msm * 100
    sdVal = subDf.std(axis=0)
    sd = ((2 ** sdVal - 1) + (1 - 2 ** (-sdVal))) / 2 * 100
    sem = sd / np.sqrt(n)
    
    tmtc_LoadingBias = pd.DataFrame({'TMTc_group': list(df.columns),
                                     'Mean[%]': list(avg),
                                     'SD[%]': list(sd),
                                     'SEM[%]': list(sem),
                                     '#Proteins': [n]*len(df.columns)})
    
    return tmtc_LoadingBias



def TMTc_normalization(df):
    ################################################
    # Normalization (i.e. loading-bias correction) #
    ################################################

    res = df.copy()

    # First, get a subset for calculating normalization factors (same as loading-bias calculation)
    # Note that this subset is 1) divided by row-wise mean and then 2) log2-transformed
    subDf = TMTc_getSubset(df)

    # Calculate normalization factors for samples (reporters)
    sm = subDf.mean(axis=0)

    target = np.mean(sm)
    normFactors = sm - target

    # Normalize the input dataframe, df (in log2-scale and then scale-back)
    rowMeans = res.mean(axis=1)
    res = np.log2(res.divide(rowMeans, axis=0).replace(0, np.nan))
    res = res - normFactors
    res = 2 ** res
    res = res.multiply(rowMeans, axis=0)
    # After the normalization, no need to show loading-bias again (should be 100% for all samples)

    return res
