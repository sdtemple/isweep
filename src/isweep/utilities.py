# imports
import numpy as np
import pandas as pd
from copy import deepcopy

### binning ###

def read_bins(file):
    '''Read *.bins file

    Parameters
    ----------
    file : str
        Input file name

    Returns
    -------
    array-like
        Increasing floats in centiMorgans
    '''

    bins = pd.read_csv(file, header = None)
    ab = bins[0]

    return list(ab)

def bin_ibd_segments(ell, ab):
    """Put ibd segments into bins

    Parameters
    ----------
    ell : numpy.array
        ibd segments
    ab : numpy.array
        Increasing floats in centiMorgans

    Returns
    -------
    numpy.array
        Observed counts for ibd segment bins
    """

    M = len(ab)
    obs = []
    for m in range(1, M):
        idx = ell < ab[m]
        obs.append(sum(idx))
        ell = ell[ell >= ab[m]]

    return np.array(obs)

### formatting vectors for plotting ###

def big_format_distribution(distr, counts):
    '''Reformat a vector for plotting with matplotlib, seaborn

    Parameters
    ----------
    distr : array-like
        Vector of realizations (lengths, times, etc.)
    counts : array-like
        Vector of realization multiplicities

    Returns
    -------
    numpy.array
        Adds copies of realization if multiplicity > 1
    '''
    nobs = sum(counts)
    dist = np.zeros(nobs)
    nlen = len(counts)
    itr = 0
    for i in range(nlen):
        ctr = counts[i]
        dis = distr[i]
        while ctr > 0:
            ctr -= 1
            dist[itr] = dis
            itr += 1
    return np.array(dist)

### effective population size ###

from math import floor, exp

def extend_vector(vec, val, mxn):
    '''Extend NumPy array

    Parameters
    ----------
    vec : array-like
        NumPy array
    val : float
        Impute value
    mxn : int
        Array size

    Returns
    -------
    array-like
        NumPy array
    '''
    num = mxn - (len(vec) - 1)
    ext = np.array([val for i in range(num)])
    return np.concatenate((vec, ext))

def zero_shift_Ne(Ne):
    '''Shift the Ne dictionary to start at generation 0

    Parameters
    ----------
    Ne : dict
        Effective population sizes

    Return
    ------
    dict
        Effective population sizes, with keys shifted
    '''
    minG = min(Ne.keys())
    return {(k-minG):v for k,v in Ne.items()}

def cut_Ne(Ne, cut):
    '''Shrink the Ne dictionary to end at specified generational time

    Parameters
    ----------
    Ne : dict
        Effective population sizes
    cut : float
        Generational time

    Returns
    -------
    dict
        Effective population sizes
    '''
    return {k:v for k,v in Ne.items() if k <= cut}

def fill_Ne(Ne):
    '''Fill in effective sizes for non-specified generations

    Parameters
    ----------
    Ne : dict
        Effective population sizes

    Returns
    -------
    dict
        Modified Ne dictionary
    '''
    assert min(Ne.keys())==0
    Me=deepcopy(Ne)
    values=list(Ne.values())
    keys=list(Ne.keys())
    keys=sorted(keys,reverse=True)
    maxkey=max(keys)
    k = 0
    j = 1
    while j <= maxkey:
        prev = values[k]
        try:
            Me[j]
            k += 1
        except KeyError:
            Me[j] = prev
        j += 1
    return dict(sorted(Me.items(),key=lambda x:x[0]))

def to_max_Ne(Ne, maxG = 500):
    '''Fill in effective sizes to a maximum assuming constant in deep past

    Parameters
    ----------
    Ne : dict
        Effective population sizes
    maxG : numeric
        Maximum generation

    Returns
    -------
    Modified Ne dictionary
    '''
    maxG = int(float(maxG))
    lastG = max(Ne.keys())
    Me = deepcopy(Ne)
    if lastG >= maxG:
        pass
    else:
        val = Me[lastG]
        for key in range(lastG+1,maxG+1):
            Me[key] = val
    return Me

def read_Ne(file):
    '''Read *.ne file

    Parameters
    ----------
    file : str
        Input file name

    Returns
    -------
    dict
        dict[generation] = size
    '''

    Ne = {}
    f = open(file, 'r')
    f.readline() # header
    for line in f:
        g, size = line.split('\t')[:2]
        Ne[int(g)] = int(float(size))
    f.close()

    return Ne

def write_Ne(Ne, output_file):
    '''Write Ne dictionary to .ne file

    Parameters
    ----------
    Ne : dict
        Effective population sizes
    output_file : str
        File name

    Returns
    -------
    None
    '''
    with open(output_file, 'w') as f:
        f.write('GEN\tNE\n')
        for k in Ne.keys():
            f.write(str(k)); f.write('\t'); f.write(str(Ne[k])); f.write('\n')
    return None

def extend_Ne(file, maxg, output):
    '''Extend *.ne file

    Constant size until user-specificied generation. Creates a new file.

    Parameters
    ----------
    file : str
        Input file name
    maxg : int
        Maximum generation
    output : str
        Output file name

    Returns
    -------
    None
    '''

    with open(file, 'r') as f:
        lines = f.readlines()
        g, size = lines[-1].strip().split('\t')[:2]

    g = int(g) + 1
    if maxg >= g:
        with open(output, 'w') as o:
            for line in lines:
                o.write(line)
            while maxg >= g:
                o.write(str(g) + '\t' + size + '\n')
                g += 1
    else:
        raise ValueError('*.ne file contains coalescent times greater than maxg')

    return None

def make_constant_Ne(file, size, maxg):
    '''Create *.ne file for constant size population

    Parameters
    ----------
    file: str
        Output file name
    size : float
        Effective population size
    maxg : int
        Maximum generation

    Returns
    -------
    None
    '''

    size = floor(size)
    maxg = floor(maxg)
    with open(file, 'w') as f:
        f.write('GEN\tNE\n')
        f.write(str(0)); f.write('\t'); f.write(str(size)); f.write('\n')
        for i in range(1, int(maxg) + 1):
            f.write(str(i)); f.write('\t'); f.write(str(size)); f.write('\n')

    return None

def make_exponential_Ne(file, size, maxg, rate):
    '''Create *.ne file for exponentially growing population

    Parameters
    ----------
    file : str
        Output file name
    size : float
        Effective population size at generation 0
    maxg : list (int)
        Maximum generation (s)
    rate : list (float)
        Exponential growth rate(s)

    Returns
    -------
    None
    '''

    assert len(rate) == len(maxg)
    size = floor(size)
    with open(file, 'w') as f:
        f.write('GEN\tNE\n')
        f.write(str(0)); f.write('\t'); f.write(str(size)); f.write('\n')
        g1 = 0 + 1
        for j in range(len(rate)):
            g2 = floor(maxg[j])
            r = rate[j]
            for i in range(g1, g2 + 1):
                size = floor(size * exp(r))
                f.write(str(i)); f.write('\t'); f.write(str(size)); f.write('\n')
            g1 = g2 + 1

    return None

def make_bottleneck_Ne(file, isize, bsize, rt1, rt2, rt3):
    '''Create *.ne file for bottlneck population

    Parameters
    ----------
    file : str
        Output file name
    isize : float
        Initial population size
    bsize : float
        Bottleneck population size
    rt1 : tuple
        First phase (exponential rate, time)
    rt2 : tuple
        Second phase (exponential rate, time)
    rt3 : tuple
        Third phase (exponential rate, time)

    Returns
    -------
    None
    '''

    # initialize
    size = []
    r1 = rt1[0]
    t1 = rt1[1]
    r2 = rt2[0]
    t2 = rt2[1]
    r3 = rt3[0]
    t3 = rt3[1]

    # first growth phase
    for i in range(t1+1):
        size.append(isize * exp(r1 * i))

    # second growth phase
    for j in range(1,t2):
        size.append(isize * exp(r1 * t1) * exp(r2 * j))

    # bottleneck
    size.append(bsize)

    # third growth phase
    for k in range(1,t3+1):
        size.append(bsize * exp(r3 * k))

    # output
    size.reverse()
    with open(file, 'w') as f:
        f.write('GEN\tNE\n')
        for l in range(len(size)):
            f.write(str(l))
            f.write('\t')
            f.write(str(size[l]))
            f.write('\n')

    return None

def make_two_growth_phases_Ne(file, isize, rt1, rt2):
    '''Create *.ne file for two growth phases population

    Parameters
    ----------
    file : str
        Output file name
    isize : float
        Initial population size
    rt1 : tuple
        First phase (exponential rate, time)
    rt2 : tuple
        Second phase (exponential rate, time)

    Returns
    -------
    None
    '''

    # initialize
    size = []
    r1 = rt1[0]
    t1 = rt1[1]
    r2 = rt2[0]
    t2 = rt2[1]

    # first growth phase
    for i in range(t1+1):
        size.append(isize * exp(r1 * i))

    # second growth phase
    for j in range(1,t2):
        size.append(isize * exp(r1 * t1) * exp(r2 * j))

    # output
    size.reverse()
    with open(file, 'w') as f:
        f.write('GEN\tNE\n')
        for l in range(len(size)):
            f.write(str(l))
            f.write('\t')
            f.write(str(size[l]))
            f.write('\n')

    return None

def make_three_growth_phases_Ne(file, isize, rt1, rt2, rt3):
    '''Create *.ne file three growth phases population

    Parameters
    ----------
    file : str
        Output file name
    isize : float
        Initial population size
    rt1 : tuple
        First phase (exponential rate, time)
    rt2 : tuple
        Second phase (exponential rate, time)
    rt3 : tuple
        Third phase (exponential rate, time)

    Returns
    -------
    None
    '''

    # initialize
    size = []
    r1 = rt1[0]
    t1 = rt1[1]
    r2 = rt2[0]
    t2 = rt2[1]
    r3 = rt3[0]
    t3 = rt3[1]

    # first growth phase
    for i in range(t1+1):
        size.append(isize * exp(r1 * i))

    # second growth phase
    for j in range(1,t2):
        size.append(isize * exp(r1 * t1) * exp(r2 * j))

    # third growth phase
    for k in range(1,t3+1):
        size.append(isize * exp(r1 * t1) * exp(r2 * j) * exp(r3 * k))

    # output
    size.reverse()
    with open(file, 'w') as f:
        f.write('GEN\tNE\n')
        for l in range(len(size)):
            f.write(str(l))
            f.write('\t')
            f.write(str(size[l]))
            f.write('\n')

    return None
