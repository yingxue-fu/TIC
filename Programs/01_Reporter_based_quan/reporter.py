import numpy as np
import os, sys, re
import pandas as pd
from pyteomics import ms2, mzxml
from utils import progressBar


def extractReporters(files, df, params, **kwargs):
    # Input arguments
    # files: mzXML or ms2 files to be quantified
    # df: dataframe of ID.txt file
    # params: parameters

    if "sig126" in kwargs:
        print("\n  Refined extraction of TMT reporter ion peaks")
    else:
        print("\n  Extraction of TMT reporter ion peaks")

    dictQuan = {}
    for file in files:
        print("    Working on {}".format(os.path.basename(file)))
        ext = os.path.splitext(file)[-1]
        if ext == ".mzXML":
            reader = mzxml.MzXML(file)  # mzXML file reader
        elif ext == ".ms2":
            reader = ms2.IndexedMS2(file)  # MS2 file reader
        else:
            sys.exit(" Currently, either .mzXML or .ms2 file is supported")

        # Extraction of TMT reporter ions in each fraction
        scans = list(df['Scan'][df['Files'] == file].unique())
        progress = progressBar(len(scans))
        for scan in scans:
            progress.increment()
            spec = reader[str(scan)]
            res = getReporterIntensity(spec, params, **kwargs)  # Array of reporter m/z and intensity values
            key = file + "_" + str(scan)
            dictQuan[key] = res

    # Create a dataframe of quantification data
    reporters = params["tmt_channels_used"].split(";")
    colNames = [re.sub("sig", "mz", i) for i in reporters] + reporters
    res = pd.DataFrame.from_dict(dictQuan, orient='index', columns=colNames)

    # Summary of quantified TMT reporter ions
    print()
    reporterSummary = getReporterSummary(res, reporters)
    nTot = len(res)
    for reporter in reporters:
        n = reporterSummary[reporter]["nPSMs"]
        print("    %s\t%d (%.2f%%) matched" % (reporter, n, n / nTot * 100))

    return res, reporterSummary


def getReporterIntensity(spec, params, **kwargs):
    tol = 10
    reporterNames = params["tmt_channels_used"].split(";")
    mzArray = []
    intensityArray = []

    for reporter in reporterNames:
        if reporter in kwargs:
            mz = getReporterMz(reporter) * (1 + kwargs[reporter]["meanMzShift"] / 1e6)
            tol = kwargs[reporter]["sdMzShift"] * float(params['tmt_peak_extraction_second_sd'])
        else:
            mz = getReporterMz(reporter)

        lL = mz - mz * tol / 1e6
        uL = mz + mz * tol / 1e6
        ind = np.where((spec["m/z array"] >= lL) & (spec["m/z array"] <= uL))[0]
        if len(ind) == 0:
            mz = 0
        elif len(ind) == 1:
            ind = ind[0]
            mz = spec["m/z array"][ind]
        elif len(ind) > 1:
            if params['tmt_peak_extraction_method'] == '2':
                ind2 = np.argmin(abs(mz - spec["m/z array"][ind]))
                ind = ind[ind2]
                mz = spec["m/z array"][ind]
            else:
                ind2 = np.argmax(spec["intensity array"][ind])
                ind = ind[ind2]
                mz = spec["m/z array"][ind]
        if lL <= mz < uL:
            intensity = spec["intensity array"][ind]
        else:
            intensity = 0
        mzArray.append(mz)
        intensityArray.append(intensity)

    outArray = mzArray + intensityArray
    return outArray


def getReporterMz(name):
    if name == "sig126":
        return 126.127726
    elif name == "sig127" or name == "sig127N":
        return 127.124761
    elif name == "sig127C":
        return 127.131081
    elif name == "sig128N":
        return 128.128116
    elif name == "sig128" or name == "sig128C":
        return 128.134436
    elif name == "sig129" or name == "sig129N":
        return 129.131471
    elif name == "sig129C":
        return 129.137790
    elif name == "sig130N":
        return 130.134825
    elif name == "sig130" or name == "sig130C":
        return 130.141145
    elif name == "sig131" or name == "sig131N":
        return 131.138180
    elif name == "sig131C":
        return 131.144500
    elif name == "sig132N":
        return 132.141535
    elif name == "sig132C":
        return 132.147855
    elif name == "sig133N":
        return 133.144890
    elif name == "sig133C":
        return 133.151210
    elif name == "sig134N":
        return 134.148245
    elif name == "sig134C":
        return 134.1545646
    elif name == "sig135N":
        return 135.1515995
    else:
        sys.exit("  {} is an incorrect reporter name. Please check the parameter file again".format(name))


def getReporterSummary(df, reporters):
    print("  Summary of quantified TMT reporter ions")
    res = {}
    for reporter in reporters:
        res[reporter] = {}
        reporterMz = getReporterMz(reporter)
        measuredMz = df[reporter.replace("sig", "mz")]
        measuredMz = measuredMz[measuredMz > 0]
        n = len(measuredMz)
        meanMzShift = ((measuredMz - reporterMz) / reporterMz * 1e6).mean()
        sdMzShift = ((measuredMz - reporterMz) / reporterMz * 1e6).std(ddof=1)
        res[reporter]['nPSMs'] = n
        res[reporter]['meanMzShift'] = meanMzShift
        res[reporter]['sdMzShift'] = sdMzShift

    return res
