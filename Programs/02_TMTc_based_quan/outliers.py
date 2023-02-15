import numpy as np
import numpy.ma as ma
from scipy.stats import t


def outlierRemoval(df, alpha):
    n = len(df)
    nOutliers = int(np.round(n * 0.2))

    # e.g.,   n = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10], then
    # nOutliers = [0, 0, 1, 1, 1, 1, 1, 2, 2, 2]
    # So, the rule becomes as follows,
    # 1. n < 3, no outlier removal
    # 2. 3 <= n <= 7, Dixon's Q-test for a single outlier removal
    # 3. n >= 8, ESD test for the removal of (potentially) multiple outliers

    indArray = []
    if nOutliers > n - 2:
        nOutliers = n - 2
    if nOutliers > 1:
        for i in range(df.shape[1]):
            ind = ESDtest(df.iloc[:, i], alpha, nOutliers)
            indArray.extend(ind)
    else:
        for i in range(df.shape[1]):
            ind = Qtest(df.iloc[:, i], alpha, nOutliers)
            indArray.extend(ind)

    # PSMs including one or more outliers will not be considered for the subsequent quantification
    indArray = list(set(indArray))    # Indices of outliers across all reporters
    df.drop(df.index[indArray], axis=0, inplace=True)
    return df


def ESDtest(x, alpha, maxOLs):
    xm = ma.array(x)
    n = len(xm)
    R, L, minds = [], [], []
    for i in range(maxOLs):
        # Compute mean and std of x
        xmean = xm.mean()
        xstd = xm.std(ddof=1)

        # Find maximum deviation
        rr = np.abs((xm - xmean) / xstd)
        minds.append(np.argmax(rr))
        R.append(rr[minds[-1]])
        p = 1.0 - alpha / (2.0 * (n - i))
        perPoint = t.ppf(p, n - i - 2)
        L.append((n - i - 1) * perPoint / np.sqrt((n - i - 2 + perPoint ** 2) * (n - i)))

        # Mask that value and proceed
        xm[minds[-1]] = ma.masked

    # Find the number of outliers
    ofound = False
    for i in range(maxOLs - 1, -1, -1):
        if R[i] > L[i]:
            ofound = True
            break

    # Prepare return value
    if ofound:
        return minds[0: i + 1]    # There are outliers
    else:
        return []    # No outliers could be detected


# Dictionary used for Dixon's Q-test
def Qtest(data, left=True, right=True, alpha=0.05):
    """
    From https://sebastianraschka.com/Articles/2014_dixon_test.html#implementing-a-dixon-q-test-function

    Keyword arguments:
        data = A ordered or unordered list of data points (int or float).
        left = Q-test of minimum value in the ordered list if True.
        right = Q-test of maximum value in the ordered list if True.
        q_dict = A dictionary of Q-values for a given confidence level,
            where the dict. keys are sample sizes N, and the associated values
            are the corresponding critical Q values. E.g.,
            {3: 0.97, 4: 0.829, 5: 0.71, 6: 0.625, ...}
    Returns a list of 2 values for the outliers, or None.
    E.g.,
       for [1,1,1] -> [None, None]
       for [5,1,1] -> [None, 5]
       for [5,1,5] -> [1, None]

    """

    q90 = [0.941, 0.765, 0.642, 0.56, 0.507, 0.468, 0.437,
           0.412, 0.392, 0.376, 0.361, 0.349, 0.338, 0.329,
           0.32, 0.313, 0.306, 0.3, 0.295, 0.29, 0.285, 0.281,
           0.277, 0.273, 0.269, 0.266, 0.263, 0.26
           ]
    Q90 = {n: q for n, q in zip(range(3, len(q90) + 1), q90)}
    q95 = [0.97, 0.829, 0.71, 0.625, 0.568, 0.526, 0.493, 0.466,
           0.444, 0.426, 0.41, 0.396, 0.384, 0.374, 0.365, 0.356,
           0.349, 0.342, 0.337, 0.331, 0.326, 0.321, 0.317, 0.312,
           0.308, 0.305, 0.301, 0.29
           ]
    Q95 = {n: q for n, q in zip(range(3, len(q95) + 1), q95)}
    q99 = [0.994, 0.926, 0.821, 0.74, 0.68, 0.634, 0.598, 0.568,
           0.542, 0.522, 0.503, 0.488, 0.475, 0.463, 0.452, 0.442,
           0.433, 0.425, 0.418, 0.411, 0.404, 0.399, 0.393, 0.388,
           0.384, 0.38, 0.376, 0.372
           ]
    Q99 = {n: q for n, q in zip(range(3, len(q99) + 1), q99)}

    if isinstance(data, list):
        pass
    else:
        x = list(data)

    if alpha == 0.1:
        q_dict = Q90
    elif alpha == 0.05:
        q_dict = Q95
    elif alpha == 0.01:
        q_dict = Q99

    assert(left or right), 'At least one of the variables, `left` or `right`, must be True.'
    assert(len(data) >= 3), 'At least 3 data points are required'
    assert(len(data) <= max(q_dict.keys())), 'Sample size too large'

    sdata = sorted(data)
    Q_mindiff, Q_maxdiff = (0,0), (0,0)

    if left:
        Q_min = (sdata[1] - sdata[0])
        try:
            Q_min /= (sdata[-1] - sdata[0])
        except ZeroDivisionError:
            pass
        Q_mindiff = (Q_min - q_dict[len(data)], sdata[0])

    if right:
        Q_max = abs((sdata[-2] - sdata[-1]))
        try:
            Q_max /= abs((sdata[0] - sdata[-1]))
        except ZeroDivisionError:
            pass
        Q_maxdiff = (Q_max - q_dict[len(data)], sdata[-1])

    if not Q_mindiff[0] > 0 and not Q_maxdiff[0] > 0:
        outliers = []
    elif Q_mindiff[0] == Q_maxdiff[0]:
        outliers = [Q_mindiff[1], Q_maxdiff[1]]
    elif Q_mindiff[0] > Q_maxdiff[0]:
        outliers = [Q_mindiff[1]]
    else:
        outliers = [Q_maxdiff[1]]

    outlierInd = [i for i, v in enumerate(data) if v in outliers]
    # survivedInd = np.setdiff1d(range(len(data)), outlierInd)

    return outlierInd
