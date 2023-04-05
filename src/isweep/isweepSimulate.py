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

def simulate_ibd_isweep(n, s, p0, Ne, long_ibd=2.0, short_ibd=1.0, random_walk=True, one_step_model='m', tau0=0, ploidy=2, record_dist=True, pairwise_output=True):
    '''ibd segments from a coalescent with selection

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
    long_ibd, short_ibd : float
        cM length threshold
    random_walk : bool
        True for random walk
    one_step_model : str
        'm', 'a', 'd', or 'r'
    tau0 : int
        Generation when neutrality begins
    ploidy : int
        1 for haploid or 2 for diploid
    record_dist : bool
        To save tract length and coalescent time distributions or not to (default True)
    pairwise_output : bool
        To save pairwise segments or not to (default True)

    Returns
    -------
    tuple(s)
        (all, adaptive allele, non-adaptive allele) then pairwise segments
        Each tuple is (number of tracts, group sizes, length distr., time distr., count distr.)
    '''

    assert ploidy in [1,2]
    assert long_ibd > 0
    assert short_ibd > 0
    assert long_ibd >= short_ibd
    assert one_step_model in ['m','a','r','d']
    assert p0 <= 1
    assert p0 >= 0

    global H, H1, H0, ldist, ldist1, ldist0, tdist1, tdist0, tdist, ddist1, ddist0, ddist
    global pairwise_segments
    pairwise_segments = []
    H1 = nx.Graph()
    H0 = nx.Graph()
    ldist1 = []
    ldist0 = []
    tdist1 = []
    tdist0 = []
    ddist1 = []
    ddist0 = []

    # renaming
    p = p0
    t = tau0
    mdl = one_step_model
    stoc = random_walk
    Ne = cut_Ne(to_max_Ne(fill_Ne(Ne),500),2000)

    # should p0 be fixed
    if (p0 == 0) or (p0 == 1):
        out = simulate_ibd(n, Ne,
                           long_ibd, short_ibd,
                           ploidy, continuous_time,
                           record_dist, pairwise_output
                          )
        # numpy-ify
        return (out,
                (np.nan,np.array([]),np.array([]),np.array([]),np.array([])),
                (np.nan,np.array([]),np.array([]),np.array([]),np.array([])),
                tuple(pairwise_segments)
               )

    # calculating structured demographies
    ps, Ns, xs = walk_variant_backward(s, p, Ne, stoc, mdl, tau0, ploidy)
    n = int(float(n))
    n = n * ploidy
    G = max(Ne.keys())
    n1 = floor(n * ps[0])
    n0 = n - n1
    N1 = [ceil(i) for i in Ns * ps]
    N1 = {i:N1[i] for i in range(len(N1))}
    qs = extend_vector(1 - ps, 1, G)
    Ns = np.array(list(Ne.values()))
    N0 = [ceil(i) for i in Ns * qs]
    N0 = {i:N0[i] for i in range(len(N0))}

    # arrival times
    times1 = simulate_coalgen(n1, N1, ploidy, to_tmrca=False)
    interiors1 = [Node(i, 0, pairwise_output) for i in range(n1)]
    interiors0 = [Node(i, 0, pairwise_output) for i in range(n1,n1+n0)]

    # pairwise comparisons pop 1
    m = n1
    indxs = [i for i in range(m)]
    itr = n1 + n0
    for t in times1:
        m -= 1
        new_node = Node(itr, t, pairwise_output)
        sindx = sorted(two_randint(m))
        s1 = sindx[0]
        s2 = sindx[1]
        new_node.pair_coalesce(interiors1[s1], interiors1[s2], long_ibd, short_ibd, 1, record_dist, pairwise_output)
        interiors1.append(new_node)
        interiors1.pop(s1)
        interiors1.pop(s2-1)
        indxs.pop(s1)
        indxs.pop(s2-1)
        indxs.append(itr)
        itr += 1
    numTracts1 = sum([node._num_tracts for node in interiors1])
    ibdGroupSizes1 = sorted([len(c) for c in nx.connected_components(H1)],reverse=True)

    # pairwise comparisons p0
    maxt1=max([node._time for node in interiors1])
    N00 = {key:val for key,val in N0.items() if key <= maxt1}
    N01 = {key:val for key,val in N0.items() if key > maxt1}
    times0 = simulate_coalgen(n0, N00, ploidy, to_tmrca=False)
    m = n0
    indxs = [i for i in range(m)]
    for t in times0:
        m -= 1
        new_node = Node(itr, t, pairwise_output)
        sindx = sorted(two_randint(m))
        s1 = sindx[0]
        s2 = sindx[1]
        new_node.pair_coalesce(interiors0[s1], interiors0[s2], long_ibd, short_ibd, 0, record_dist, pairwise_output)
        interiors0.append(new_node)
        interiors0.pop(s1)
        interiors0.pop(s2-1)
        indxs.pop(s1)
        indxs.pop(s2-1)
        indxs.append(itr)
        itr += 1
    numTracts0 = sum([node._num_tracts for node in interiors0])
    ibdGroupSizes0 = sorted([len(c) for c in nx.connected_components(H0)],reverse=True)

    # arrival times, pairwise comparisons to TMRCA
    H = nx.compose(H1,H0)
    ldist = [*ldist1, *ldist0]
    tdist = [*tdist1, *tdist0]
    ddist = [*ddist1, *ddist0]
    interiors = interiors1 + interiors0
    m = len(interiors)
    indxs = [i for i in range(m)]
    if len(N01.keys()) > 0:
        times2 = simulate_coalgen(m, N01, ploidy, to_tmrca=True)
        for t in times2:
            m -= 1
            new_node = Node(itr, t, pairwise_output)
            sindx = sorted(two_randint(m))
            s1 = sindx[0]
            s2 = sindx[1]
            new_node.pair_coalesce(interiors[s1], interiors[s2], long_ibd, short_ibd, 2, record_dist, pairwise_output)
            interiors.append(new_node)
            interiors.pop(s1)
            interiors.pop(s2-1)
            indxs.pop(s1)
            indxs.pop(s2-1)
            indxs.append(itr)
            itr += 1
    numTracts = sum([node._num_tracts for node in interiors])
    ibdGroupSizes = sorted([len(c) for c in nx.connected_components(H)],reverse=True)

    # numpy-ify
    return ((numTracts,np.array(ibdGroupSizes),np.array(ldist),np.array(tdist),np.array(ddist)),
            (numTracts1,np.array(ibdGroupSizes1),np.array(ldist1),np.array(tdist1),np.array(ddist1)),
            (numTracts0,np.array(ibdGroupSizes0),np.array(ldist0),np.array(tdist0),np.array(ddist0)),
            tuple(pairwise_segments)
           )

def simulate_ibd_isweep_independent(n, s, p0, Ne, long_ibd=2, random_walk=True, one_step_model='m', tau0=0, ploidy=2):
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
    if p0 == 0 or p0 == 1:
        out = simulate_ibd_independent(n, Ne, long_ibd, random_walk, one_step_model, tau0, ploidy)
        return (out[0], np.array([]), np.array([]), out[1], np.array([]), np.array([]))
    Ne = cut_Ne(to_max_Ne(fill_Ne(Ne),500),2000)
    ps, Ns, xs = walk_variant_backward(s, p0, Ne, random_walk, one_step_model, tau0, ploidy)
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

def probability_ibd_isweep(s, p0, Ne, one_step_model = 'm', tau0 = 0, long_ibd = 2, ploidy = 2):
    '''Approximate probability of ibd given a sweep model

    Parameters
    ----------
    s : float
        Selection coefficient
    p0 : float
        Variant frequency at generation 0
    Ne : dict
        Effective population sizes
    one_step_model : str
        'm', 'a', 'd', or 'r'
    tau0 : int
        Generation when neutrality begins
    long_ibd : float
        cM length threshold
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
    if p0==1 or p0==0:
        Ne = cut_Ne(to_max_Ne(fill_Ne(Ne),500),1000)
        b = np.array(list(Ne.values()))
        a = extend_vector(np.array([1]),1,max(Ne.keys()))
        d = probability_ibd(a, b, long_ibd, ploidy)
        f = 0
    else:
        Ne = cut_Ne(to_max_Ne(fill_Ne(Ne),500),1000)
        a, b, c = walk_variant_backward(s, p0, Ne, False, one_step_model, tau0, ploidy)
        a = extend_vector(a, a[-1], max(Ne.keys()))
        b = extend_vector(b, b[-1], max(Ne.keys()))
        d = probability_ibd(a, b, long_ibd, ploidy)
        f = probability_ibd(1-a, b, long_ibd, ploidy)
    return d + f
