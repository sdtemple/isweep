#!/bin/python

# imports
import numpy as np
import pandas as pd
import networkx as nx
from scipy.stats import binom
from random import randint, sample
from math import floor, exp, log, ceil
from numpy.random import exponential as Exp
from copy import deepcopy
from .cis.cisUtilities import *
from .cis.coalescentIBD import *
# from .cis.cisUtilities import *
# from .cis.coalescentIBD import *

def simulate_ibd_isweep_independent(n, s, p0, Ne, long_ibd=2, random_walk=True, one_step_model='a', tau0=0, sv=-0.01, ploidy=2):
    """Simulator for independent ibd segment lengths in recent sweep scenario

    Parameters
    ----------
    n : int
        Sample size (individuals)
    s : float
        Selection coefficient
    p0 : float
        Variant frequency at generation 0
    Ne : dict
        Effective population sizes
    long_ibd : float
        cM length threshold
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
    tuple (NumPy arrays)
        (all, selected, non-selected) ibd lengths + (all, selected, non-selected) coalescent times
    """

    assert ploidy in [1,2]
    assert p0 >= 0
    assert p0 <= 1
    assert sv < 1
    if p0 == 0 or p0 == 1:
        out = simulate_ibd_independent(n, Ne, long_ibd, ploidy)
        return (out[0], np.array([]), np.array([]), out[1], np.array([]), np.array([]))
    Ne = cut_Ne(to_max_Ne(fill_Ne(Ne),500),2000)
    ps, Ns, xs = walk_variant_backward(s, p0, Ne, random_walk, one_step_model, tau0, sv, ploidy)
    n = int(n) * ploidy
    n1 = floor(n * p0)
    n0 = floor(n * (1-p0))
    m1 = n1 * (n1 - 1) / 2
    m0 = n0 * (n0 - 1) / 2
    mass1 = probability_quasi_geometric(ps, Ns, ploidy)
    mass0 = probability_quasi_geometric(1-ps, Ns, ploidy)
    geom1 = simulate_quasi_geometric(m1, mass1)
    geom0 = simulate_quasi_geometric(m0, mass0)
    geom1 = geom1[geom1 != np.Inf]
    geom0 = geom0[geom0 != np.Inf]
    ell1 = simulate_erlang_segments(geom1)
    ell0 = simulate_erlang_segments(geom0)
    geom1 = geom1[ell1 >= long_ibd]
    geom0 = geom0[ell0 >= long_ibd]
    ell1 = ell1[ell1 >= long_ibd]
    ell0 = ell0[ell0 >= long_ibd]
    ell = np.concatenate((ell1, ell0))
    geom = np.concatenate((geom1, geom0))

    return ell, ell1, ell0, geom, geom1, geom0

def probability_ibd_isweep(s, p0, Ne, long_ibd = 2, one_step_model = 'a', tau0 = 0, sv=-0.01, ploidy = 2):
    '''Approximate probability of ibd given a sweep model

    Parameters
    ----------
    s : float
        Selection coefficient
    p0 : float
        Variant frequency at generation 0
    Ne : dict
        Effective population sizes
    long_ibd : float
        cM length threshold
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
    float
        approx P(\ell > c) where \ell is ibd length
    '''
    assert ploidy in [1,2]
    assert p0 >= 0
    assert p0 <= 1
    assert sv < 1
    if p0==1 or p0==0:
        Ne = cut_Ne(to_max_Ne(fill_Ne(Ne),500),1000)
        b = np.array(list(Ne.values()))
        a = extend_vector(np.array([1]),1,max(Ne.keys()))
        d = probability_ibd(a, b, long_ibd, ploidy)
        f = 0
    else:
        Ne = cut_Ne(to_max_Ne(fill_Ne(Ne),500),1000)
        a, b, c = walk_variant_backward(s, p0, Ne, False, one_step_model, tau0, sv, ploidy)
        a = extend_vector(a, a[-1], max(Ne.keys()))
        b = extend_vector(b, b[-1], max(Ne.keys()))
        d = probability_ibd(a, b, long_ibd, ploidy)
        f = probability_ibd(1-a, b, long_ibd, ploidy)
    return d + f
