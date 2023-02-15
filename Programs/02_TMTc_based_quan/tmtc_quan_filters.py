import sys


def TMTc_getFileteredIndexes(df, method, threshold):
    # This module produces indexes to be removed (i.e., filtered indexes)
    
    df = df.loc[:, df.columns.str.contains(pat = 'sig1')]
    
    if method == '1':  # Minimum-based filter
        idx = df[df.min(axis=1) < threshold].index
    elif method == '2':  # Maximum-based filter
        idx = df[df.max(axis=1) < threshold].index
    elif method == '3':  # Mean-based filter
        idx = df[df.mean(axis=1) < threshold].index
    elif method == '4':  # Median-based filter
        idx = df[df.median(axis=1) < threshold].index
    else:
        sys.exit("  Please check 'min_intensity_method' parameter. It should be 0, 1, 2, 3, or 4")

    return idx


def TMTc_filterByIntensity(df, methods, thresholds, verbose=1):
    methodsStr = ["none", "minimum intensity", "maximum intensity", "mean intensity", "median intensity"]
    res = []
    n = 0
    for i in range(len(methods)):
        #
        if methods[i] == "0":
            pass
        else:
            idx = TMTc_getFileteredIndexes(df, methods[i], float(thresholds[i]))    # "idx" is the index to be removed
            res.extend(idx.values)
        #
        if verbose:
            res = list(set(res))
            print("    Removed {} PSMs based on the intensity-based filter ({})".format(len(res) - n, methodsStr[int(methods[i])]))
            n = len(res)

    return res


def TMTc_filterPSMs(df, params):
    
    print("\n  Examining the TMTc quantification results of PSMs (n = {})".format(len(df)))
    
    methods = params["tmtc_min_intensity_method"].split(",")
    thresholds = params["tmtc_min_intensity_value"].split(",")
    
    idx = TMTc_filterByIntensity(df, methods, thresholds, 1)  # Indexes to be removed
    idx = list(set(idx))
    df = df[~df.index.isin(idx)]
    print("    Hereafter, {} PSMs will be used for the quantification".format(len(df)))

    return df
