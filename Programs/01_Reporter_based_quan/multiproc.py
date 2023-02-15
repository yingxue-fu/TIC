import os, sys, re, numpy as np, pandas as pd
from pyteomics import ms2, mzxml
import multiprocessing, tqdm
from contextlib import contextmanager
from functools import partial
from itertools import islice


@contextmanager
def poolcontext(*args, **kwargs):
    pool = multiprocessing.Pool(*args, **kwargs)
    yield pool
    pool.terminate()


def parExtractReporters(files, df, params, nCores, **kwargs):
    if "sig126" in kwargs:
        print("\n  Refined extraction of TMT reporter ion peaks")
    else:
        print("\n  Extraction of TMT reporter ion peaks")

    results = []
    with poolcontext(processes=nCores) as pool:
        for result in tqdm.tqdm(pool.imap(partial(singleExtractReporters, df=df, params=params, **kwargs), files),
                                total=len(files),
                                bar_format="    Progress: [{bar:20}] {percentage:3.0f}%"):
            results.append(result)

    # Processing the output
    res = pd.concat(results, axis=0)

    # Summary of quantified TMT reporter ions
    print()
    reporters = params["tmt_reporters_used"].split(";")
    reporterSummary = getReporterSummary(res, reporters)
    nTot = len(res)
    for reporter in reporters:
        n = reporterSummary[reporter]["nPSMs"]
        print("    %s\t%d (%.2f%%) matched" % (reporter, n, n / nTot * 100))

    return res, reporterSummary


def singleExtractReporters(file, df, params, **kwargs):
    # Input arguments
    # files: mzXML or ms2 files to be quantified
    # df: dataframe of ID.txt file
    # params: parameters

    dictQuan = {}
    ext = os.path.splitext(file)[-1]
    if ext == ".mzXML":
        reader = mzxml.MzXML(file)  # mzXML file reader
    elif ext == ".ms2":
        reader = ms2.IndexedMS2(file)  # MS2 file reader
    else:
        sys.exit(" Currently, either .mzXML or .ms2 file is supported")

    # Extraction of TMT reporter ions in each fraction
    scans = list(df['scan'][df['frac'] == file].unique())
    for scan in scans:
        spec = reader[str(scan)]
        res = getReporterIntensity(spec, params, **kwargs)  # Array of reporter m/z and intensity values
        key = file + "_" + str(scan)
        dictQuan[key] = res

    # Create a dataframe of quantification data
    reporters = params["tmt_reporters_used"].split(";")
    colNames = [re.sub("sig", "mz", i) for i in reporters] + reporters
    res = pd.DataFrame.from_dict(dictQuan, orient='index', columns=colNames)
    return res


# 2021/9/12
# While testing, "pool" function produces an error in a node of HPC
# It is still not clear why such an error occurs. I will discuss it with Karthik in RIS later
def chunks(data, nCPU):
    size = round(len(data) / nCPU)
    it = iter(data)
    res = []
    for i in range(0, len(data), size):
        res.append({k: data[k] for k in islice(it, size)})

    return res


def parSummarization(inputDict, df, params):
    nProc = round(multiprocessing.cpu_count() / 2)    # For safety, the half of available cores will be used
    listDict = chunks(inputDict, nProc)
    results = []
    with poolcontext(processes=nProc) as pool:
        for result in tqdm.tqdm(pool.imap(partial(summarization, df=df, params=params), listDict), total=len(listDict),
                                bar_format="    Progress: [{bar:20}] {percentage:3.0f}%"):
            results.append(result)
    # Processing the output
    res = pd.concat(results, axis=0)
    return res


# singleSummarization needs to be defined and implemented