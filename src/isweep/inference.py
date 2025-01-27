# imports
import gzip
import numpy as np
import pandas as pd
from random import randint
from math import floor, exp, log, ceil
from numpy.random import exponential as Exp
from scipy.stats import norm, binom
from scipy.optimize import minimize, minimize_scalar
from .coalescent import walk_variant_backward, probability_quasi_geometric, probability_erlang_segments
from .utilities import bin_ibd_segments

# read .ibd file

def read_ibd_file(ibd_file, header = 1, include_length = 1):
    '''Create list of IBD segments

    Parameters
    ----------
    ibd_file : str
        Name of text file with IBD segments
    header : int
        0 for no header
    include_length : int
        0 for no length

    Returns
    -------
    list
        (ID1, ID2, cM length) pairs
    '''
    segments = []
    with gzip.open(ibd_file, 'rt') as f:
        if header != 0:
            f.readline()
        for line in f:
            row = line.strip().split('\t')
            id1 = str(row[0]) + '_' + str(row[1])
            id2 = str(row[2]) + '_' + str(row[3])
            if include_length != 0:
                length = float(row[7])
                segments.append((id1, id2, length))
            else:
                segments.append((id1, id2))
    return segments

# goodness-of-fit statistics (estimate selection coeffficent)

def chi2_labeled_isweep(s, p0, Ne, n, obs1, obs0, ab, one_step_model = 'a', tau0 = 0, sv=-0.01, ploidy = 2):
    '''Chi-squared statistic for sweep model (labeled)

    Parameters
    ----------
    s : float
        Selection coefficient
    p0 : float
        Variant frequency at generation 0
    Ne : dict
        Effective population sizes
    n : int
        Pairs sample size
    obs1 : numpy.array
        Observed counts for ibd segment bins (labeled 1)
    obs0 : numpy.array
        Observed counts for ibd segment bins (labeled 0)
    ab : numpy.array
        Increasing floats in centiMorgans
    one_step_model : str
        'm', 'a', 'd', or 'r'
    tau0 : int
        Generation at which neutrality begins
    sv: float
        Allele frequency of standing variation
        (Default -0.01 will assume de novo sweep)
    ploidy : int
        1 for haploid or 2 for diploid

    Returns
    -------
    float
        Goodness-of-fit statistic
    '''

    # assertions
    assert len(obs1) == len(obs0)
    assert len(obs1) == (len(ab) - 1)

    # probabilities
    ps, Ns, xs = walk_variant_backward(s, p0, Ne, False, one_step_model, tau0, sv, ploidy)
    mass1 = (ps[0] ** 2) * probability_quasi_geometric(ps, Ns, ploidy)
    mass0 = ((1 - ps[0]) ** 2) * probability_quasi_geometric(1 - ps, Ns, ploidy)

    # expecteds
    E0 = probability_erlang_segments(ab, mass0)
    E1 = probability_erlang_segments(ab, mass1)
    E1 = n * E1
    E0 = n * E0

    # chi-squared goodness-of-fit
    OE1 = ((obs1 - E1) ** 2) / E1
    OE0 = ((obs0 - E0) ** 2) / E0

    return np.sum(OE1) + np.sum(OE0)

def chi2_isweep(s, p0, Ne, n, obs, ab, one_step_model = 'a', tau0 = 0, sv=-0.01, ploidy = 2):
    '''Chi-squared statistic for sweep model (unlabeled)

    Parameters
    ----------
    s : float
        Selection coefficient
    p0 : float
        Variant frequency at generation 0
    Ne : dict
        Effective population sizes
    n : int
        Pairs sample size
    obs : numpy.array
        Observed counts for ibd segment bins (unlabeled)
    ab : numpy.array
        Increasing floats in centiMorgans
    one_step_model : str
        'm', 'a', 'd', or 'r'
    tau0 : int
        Generation at which neutrality begins
    sv: float
        Allele frequency of standing variation
        (Default -0.01 will assume de novo sweep)
    ploidy : int
        1 for haploid or 2 for diploid

    Returns
    -------
    float
        Goodness-of-fit statistic
    '''

    assert len(obs) == (len(ab) - 1)
    ps, Ns, xs = walk_variant_backward(s, p0, Ne, False, one_step_model, tau0, sv, ploidy)
    mass1 = (ps[0] ** 2) * probability_quasi_geometric(ps, Ns, ploidy)
    mass0 = ((1 - ps[0]) ** 2) * probability_quasi_geometric(1 - ps, Ns, ploidy)
    E0 = probability_erlang_segments(ab, mass0)
    E1 = probability_erlang_segments(ab, mass1)
    E = E1 + E0
    E = n * E
    OE = ((obs - E) ** 2) / E

    return np.sum(OE)

### bootstrap (study uncertainty) ###

def bootstrap_standard(val, boot, alpha1 = 0.025, alpha2 = 0.975):
    '''Implements the standard bootstrap interval estimator

    Parameters
    ----------
    val : float
        Parameter estimate
    boot : numpy.array
        Bootstraps
    alpha1 : float
        Percentile
    alpha2 : float
        Percentile

    Returns
    -------
    tuple
        (lower, middle, upper) interval estimator

    '''
    za = norm.ppf(alpha1)
    zb = norm.ppf(alpha2)
    sd = np.std(boot)

    return (val + za * sd, val, val + zb * sd)

def bootstrap_standard_bc(val, boot, alpha1 = 0.025, alpha2 = 0.975):
    '''Implements the standard bootstrap interval estimator (w/ bias-correction)

    Parameters
    ----------
    val : float
        Parameter estimate
    boot : numpy.array
        Bootstraps
    alpha1 : float
        Percentile
    alpha2 : float
        Percentile

    Returns
    -------
    tuple
        (lower, middle, upper) interval estimator
    '''

    za = norm.ppf(alpha1)
    zb = norm.ppf(alpha2)
    bc = np.mean(boot) - val
    sd = np.std(boot)

    return (val - bc + za * sd, val - bc, val - bc + zb * sd)

def bootstrap_percentile_bc(val, boot, alpha1 = 0.025, alpha2 = 0.975):
    '''Implement percentile-based interval estimator with bias correction

    Parameters
    ----------
    val : float
        Parameter estimate
    boot : numpy.array
        Bootstraps
    alpha1 : float
        Percentile
    alpha2 : float
        Percentile

    Returns
    -------
    tuple
        (lower, middle, upper) interval estimator
    '''
    boot = np.array(boot)
    err = boot - val
    return (
            val - np.quantile(err, alpha2),
            val - np.quantile(err, 0.5),
            val - np.quantile(err, alpha1)
           )

def bootstrap_percentile(val, boot, alpha1 = 0.025, alpha2 = 0.975):
    '''Implements percentile-based interval estimator

    Parameters
    ----------
    val : float
        Parameter estimate
    boot : numpy.array
        Bootstraps
    alpha1 : float
        Percentile
    alpha2 : float
        Percentile

    Returns
    -------
    tuple
        (lower, middle, upper) interval estimator
    '''

    vl = np.quantile(boot, alpha1)
    vu = np.quantile(boot, alpha2)
    vm = np.quantile(boot, 0.5)

    return (vl, vm, vu)

##### time to event #####

def when_count(ct, s, p0, Ne, random_walk = True, one_step_model = 'a', tau0 = 0, sv=-0.01, ploidy=2):
    '''Report when variant count reaches set value

    Parameters
    ----------
    ct : int
        Variant count
    s : float
        Selection coefficient
    p0 : float
        Variant frequency at generation 0
    Ne : dict
        Effective population sizes
    random_walk : bool
        True for random walk
    one_step_model : str
        'm', 'a', 'd', or 'r'
    tau0 : int
        Generation when neutrality begins
    sv: float
        Allele frequency of standing variation
        (Default -0.01 will assume de novo sweep)
    ploidy : int
        1 for haploid or 2 for diploid

    Returns
    -------
    int
        Generation time
    '''

    var = walk_variant_backward(s, p0, Ne, random_walk, one_step_model, tau0, sv, ploidy)[2]
    idx = np.argmax(var <= ct)
    if idx == 0: # handle variant introduction
        return len(var)

    return idx

def bootstrap_count(ct, B, boots, bootp, Ne, random_walk = True, one_step_model = 'a', tau0 = 0, sv=-0.01, ploidy=2):
    '''Parametric bootstrap for variant count generation time

    Parameters
    ----------
    ct : int
        Variant count
    B : int
        Number of bootstraps
    boots : array-like
        NumPy array of bootstraps for selection coefficient
    bootp : array-like
        NumPy array of bootstraps for variant frequency
    Ne : dict
        Effective population sizes
    random_walk : bool
        True for random walk
    one_step_model : str
        'm', 'a', 'd', or 'r'
    tau0 : int
        Generation when neutrality begins
    sv: float
        Allele frequency of standing variation
        (Default -0.01 will assume de novo sweep)
    ploidy : int
        1 for haploid or 2 for diploid

    Returns
    -------
    array-like
        NumPy array of bootstraps for variant count generation time
    '''

    assert len(boots) == len(bootp)
    boott = []
    for i in range(len(boots)):
        for b in range(B):
            boott.append(when_count(ct, boots[i], bootp[i], Ne, random_walk, one_step_model, tau0, sv, ploidy))

    return boott

def when_freq(maf, s, p0, Ne, random_walk = True, one_step_model = 'a', tau0 = 0, sv=-0.01, ploidy=2):
    '''Report when variant frequency reaches set value

    Parameters
    ----------
    maf : float
        Variant frequency
    s : float
        Selection coefficient
    p0 : float
        Variant frequency at generation 0
    Ne : dict
        Effective population sizes
    random_walk : bool
        True for random walk
    one_step_model : str
        'm', 'a', 'd', or 'r'
    tau0 : int
        Generation when neutrality begins
    sv: float
        Allele frequency of standing variation
        (Default -0.01 will assume de novo sweep)
    ploidy : int
        1 for haploid or 2 for diploid

    Returns
    -------
    int
        Generation time
    '''

    var = walk_variant_backward(s, p0, Ne, random_walk, one_step_model, tau0, sv, ploidy)[0]
    idx = np.argmax(var <= maf)
    if idx == 0: # handle variant introduction
        return len(var)

    return idx

def bootstrap_freq(maf, B, boots, bootp, Ne, random_walk = True, one_step_model = 'a', tau0 = 0, sv=-0.01, ploidy=2):
    '''Parametric bootstrap for variant frequency generation time

    Parameters
    ----------
    maf : float
        Variant frequency
    B : int
        Number of bootstraps
    boots : numpy.array
        NumPy array of bootstraps for selection coefficient
    bootp : numpy.array
        NumPy array of bootstraps for variant frequency
    Ne : dict
        Effective population sizes
    random_walk : bool
        True for random walk
    one_step_model : str
        'm', 'a', 'd', or 'r'
    tau0 : int
        Generation when neutrality begins
    sv: float
        Allele frequency of standing variation
        (Default -0.01 will assume de novo sweep)
    ploidy : int
        1 for haploid or 2 for diploid

    Returns
    -------
    numpy.array
        NumPy array of bootstraps for variant frequency generation time
    '''

    assert len(boots) == len(bootp)
    boott = []
    for i in range(len(boots)):
        for b in range(B):
            boott.append(when_freq(maf, boots[i], bootp[i], Ne, random_walk, one_step_model, tau0, sv, ploidy))

    return boott

### metrics (evaluate experimental results) ###

# def check_envelope(x, a, b):
#     '''Check if value is in an interval

#     Parameters
#     ----------
#     x : float
#         Parameter estimate
#     a : float
#         Left endpoint
#     b : float
#         Right endpoint

#     Returns
#     -------
#     bool
#         If covered
#     '''

#     if x >= a:
#         if x <= b:
#             return True

#     return False

# def make_confusion_matrix(real, pred):
#     '''Form confusion matrix for statistical classification

#     Parameters
#     ----------
#     real : array_like
#         Actual 1s and 0s
#     pred : array_like
#         Predicted 1s and 0s

#     Returns
#     -------
#     pandas DataFrame object
#         2 x 2 crosstab
#     '''

#     data = {'Actual':real, 'Predicted':pred}
#     conf = pd.DataFrame(data, columns = ['Actual', 'Predicted'])
#     matr = pd.crosstab(conf['Actual'],
#                        conf['Predicted'],
#                        rownames=['Actual'],
#                        colnames=['Predicted'])

#     return matr

# def false_omission_rate(matrix):
#     '''Compute false omission rate

#     Parameters
#     ----------
#     matrix : pandas DataFrame object
#         2 x 2 cross tab

#     Returns
#     -------
#     float
#         The false omission rate (https://en.wikipedia.org/wiki/Confusion_matrix)
#     '''

#     tn = matrix[0][0]
#     fp = matrix[1][0]
#     fn = matrix[0][1]
#     tp = matrix[1][1]

#     return fn / (fn + tn)

# def false_positive_rate(matrix):
#     '''Compute false positive rate

#     Parameters
#     ----------
#     matrix : pandas DataFrame object
#         2 x 2 cross tab

#     Returns
#     -------
#     float
#         The false positive rate (https://en.wikipedia.org/wiki/Confusion_matrix)
#     '''

#     tn = matrix[0][0]
#     fp = matrix[1][0]
#     fn = matrix[0][1]
#     tp = matrix[1][1]

#     return fp / (fp + tn)

# def true_positive_rate(matrix):
#     '''Compute true positive rate

#     Parameters
#     ----------
#     matrix : pandas DataFrame object
#         2 x 2 cross tab

#     Returns
#     -------
#     float
#         The true positive rate (https://en.wikipedia.org/wiki/Confusion_matrix)
#     '''

#     tn = matrix[0][0]
#     fp = matrix[1][0]
#     fn = matrix[0][1]
#     tp = matrix[1][1]

#     return tp / (tp + fn)
