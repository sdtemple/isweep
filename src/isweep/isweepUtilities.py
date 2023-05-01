#!/bin/python

# imports
import numpy as np
import pandas as pd

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
    ell : array-like
        ibd segments
    ab : array-like
        Increasing floats in centiMorgans

    Returns
    -------
    NumPy array
        Observed counts for ibd segment bins
    """

    M = len(ab)
    obs = []
    for m in range(1, M):
        idx = ell < ab[m]
        obs.append(sum(idx))
        ell = ell[ell >= ab[m]]

    return np.array(obs)
