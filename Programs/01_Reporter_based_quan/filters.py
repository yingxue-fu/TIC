import sys


def getFileteredIndexes(df, method, threshold, reporters):
    # This module produces indexes to be removed (i.e., filtered indexes)
    if method == '1':  # Minimum-based filter
        idx = df[(df[reporters] < threshold).any(axis=1)].index
    elif method == '2':  # Maximum-based filter
        idx = df[(df[reporters] > threshold).any(axis=1)].index
    elif method == '3':  # Mean-based filter
        idx = df[df[reporters].mean(axis=1) < threshold].index
    elif method == '4':  # Median-based filter
        idx = df[df[reporters].median(axis=1) < threshold].index
    else:
        sys.exit("  Please check 'min_intensity_method' parameter. It should be 0, 1, 2, 3, or 4")

    return idx


def filterByIntensity(df, methods, thresholds, reporters, verbose=1):
    methodsStr = ["none", "minimum intensity", "maximum intensity", "mean intensity", "median intensity"]
    res = []
    n = 0
    for i in range(len(methods)):
        if methods[i] == "0":
            pass
        else:
            idx = getFileteredIndexes(df, methods[i], float(thresholds[i]), reporters)    # "idx" is the index to be removed
            res.extend(idx.values)
        if verbose:
            res = list(set(res))
            print("    Removed {} PSMs based on the intensity-based filter ({})".format(len(res) - n, methodsStr[int(methods[i])]))
            n = len(res)

    return res


def filterPSMs(df, params):
    print("\n  Examining the extracted TMT reporter ions in PSMs")
    reporters = params["tmt_channels_used"].split(";")

    # 0. Zero-intensity filter
    n = len(df)
    df = df[(df[reporters] > 0).all(axis=1)]
    print("    Removed {} PSMs due to zero intensity at least one channel".format(n - len(df)))

    # 1. Intensity-based filtering of all PSMs
    methods = params["min_intensity_method"].split(",")
    thresholds = params["min_intensity_value"].split(",")
    idx = filterByIntensity(df, methods, thresholds, reporters, 1)  # Indexes to be removed
    idx = list(set(idx))
    df = df[~df.index.isin(idx)]
    print("    Hereafter, {} PSMs will be used for the quantification".format(len(df)))

    """
    2021/12/28
    The following filter was used in the old implementation to filter out PSMs mapped to a protein.
    However, rather than using the filter, it would be more reasonable to consider it in the summarization step.
    Therefore, the filter is commented and not going to be used

    # 2. Further filtering when only 1 or 2 PSMs are mapped to a protein
    print("    Further filtering of 1 or 2 PSMs mapped to a protein")
    methods = params["min_intensity_method_1_2_psm"].split(",")
    thresholds = params["min_intensity_value_1_2_psm"].split(",")
    idx = []
    progress = progressBar(len(prot2psm))
    for prot, psms in prot2psm.items():
        progress.increment()
        psms = df.index.join(psms, how="inner")

        if len(psms) == 0:
            continue
        elif len(psms) == 1:
            # Proteins mapped by only one PSM
            # If the PSM is filtered by the intensity-based filter, it will not be used for the quantification
            idxProt = filterByIntensity(df.loc[psms], methods, thresholds, reporters, 0)
            if len(idxProt) > 0:
                idx.extend(idxProt)
        elif len(psms) == 2:
            # Proteins mapped by two PSMs
            # Apply the intensity-based filter first
            # - If both PSMs are filtered out, they will not be used for the quantification
            # - If one of PSMs is filtered out, the PSM will not be used for the quantification
            # - If none of PSMs is filtered out, go to the next step (two PSMs can be used for the quantification)
            #   - For each PSM, check the variation (i.e., stdev) across the reporters (in log2-space)
            #   - One with smaller variation will be used for the quantification
            #   - If both PSMs have the same variation, the one with higher mean intensity will be used
            #         subDf = filterPSM1(subDf, methods, thresholds, reporters, 0)
            idxProt = filterByIntensity(df.loc[psms], methods, thresholds, reporters, 0)
            if len(idxProt) > 0:
                idx.extend(idxProt)
            else:
                psmStd = np.log2(df.loc[psms][reporters]).std(axis=1)
                psmMean = np.log2(df.loc[psms][reporters]).mean(axis=1)
                if psmStd[0] == psmStd[1]:
                    ii = np.argmin(psmMean)
                else:
                    ii = np.argmax(psmStd)
                idx.extend([psms[ii]])

    idx = list(set(idx))
    print("    Removed {} PSMs due to the larger variation than the other PSM mapped to the same protein".format(len(idx)))
    df = df[~df.index.isin(idx)]
    """

    return df
