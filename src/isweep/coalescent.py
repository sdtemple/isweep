# imports
import numpy as np
import pandas as pd
import networkx as nx
from scipy.stats import binom
from random import randint, sample
from math import floor, exp, log, ceil, sqrt
from numpy.random import exponential as Exp
from copy import deepcopy
from .utilities import *

H1 = nx.Graph()
H0 = nx.Graph()
H = nx.Graph()
ldist = []
ldist1 = []
ldist0 = []
tdist = []
tdist1 = []
tdist0 = []
ddist1 = []
ddist0 = []
ddist = []
pairwise_segments = []

### random walks ###

def walk_variant_backward(s, p0, Ne, random_walk = False, one_step_model = 'a', tau0 = 0, sv=-0.01, ploidy = 2):
    '''Variant frequencies backward in time

    Parameters
    ----------
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
    tuple
        NumPy arrays for frequencies and sizes
    '''

    # local functions
    assert ploidy in [1,2]
    assert one_step_model in ['m','a','r','d']
    assert p0 <= 1
    assert p0 >= 0
    assert sv < 1
    def haploid_bwd(p, s): # haploid is same as multiplicative (Felsenstein, 2017)
        return p / (1 + s - s * p)
    def multiplicative_bwd(p, s):
        num = 1
        dnm = 1 + (1 - p) * s
        return p * num / dnm
    def dominant_bwd(p, s):
        if s <= 0:
            return p
        a = p * s
        b = 1 + s - 2 * p * s
        c = - p
        qf = - b + sqrt((b ** 2) - 4 * a * c)
        qf = qf / 2 / a
        return qf
    def additive_bwd(p, s):
        if s <= 0:
            return p
        a = s
        b = 1 + s - 2 * p * s
        c = - p
        qf = - b + sqrt((b ** 2) - 4 * a * c)
        qf = qf / 2 / a
        return qf
    def recessive_bwd(p, s):
        if s <= 0:
            return p
        a = (1 - p) * s
        b = 1
        c = - p
        qf = - b + sqrt((b ** 2) - 4 * a * c)
        qf = qf / 2 / a
        return qf

    # one step calculation
    if ploidy == 1:
        one_step = haploid_bwd
    else:
        if one_step_model == 'a':
            one_step = additive_bwd
        elif one_step_model == 'r':
            one_step = recessive_bwd
        elif one_step_model == 'd':
            one_step = dominant_bwd
        else:
            one_step = multiplicative_bwd

    # initialize
    ps = [] # frequencies
    xs = [] # variants
    Ns = [] # sizes
    t = floor(tau0)
    p = p0
    N = Ne[0]
    x = floor(p * ploidy * N)
    Ns.append(N)
    xs.append(x)
    ps.append(p)

    if random_walk: # random walk

        for G in range(1, max(Ne.keys())+1):
            try: # population size change
                N = Ne[G]
            except KeyError:
                pass
            if G > t:
                p = one_step(p, s)
            x = int(binom.rvs(int(ploidy * N), p))
            p = x / ploidy / N
            if x < 1:
                break
            if p >= 1:
                break
            if p <= sv:
                s = 0
            ps.append(p)
            xs.append(x)
            Ns.append(N)

        return np.array([ps, Ns, xs], dtype=float) # numpy-ify

    else: # deterministic

        for G in range(1, max(Ne.keys())+1):
            try: # population size change
                N = Ne[G]
            except KeyError:
                pass
            if G > t:
                p = one_step(p, s)
            x = floor(p * ploidy * N)
            if x < 1:
                break
            if p >= 1:
                break
            if p <= sv:
                s = 0
            Ns.append(N)
            xs.append(x)
            ps.append(p)

        return np.array([np.array(ps), np.array(Ns), np.array(xs)]) # numpy-ify

def walk_variant_forward(s, pG, Ne, random_walk = False, one_step_model = 'a', tau0 = 0, ploidy = 2):
    '''Variant frequencies forward in time

    Parameters
    ----------
    s : float
        Selection coefficient
    pG : float
        Variant frequency at maximum generation
    Ne : dict
        Effective population sizes
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
    tuple
        NumPy arrays for frequencies and sizes
    '''

    # local functions
    assert ploidy in [1,2]
    def haploid_fwd(p, s): # haploid is same as multiplicative
        return p * (1 + s) / (1 + p * s)
    def additive_fwd(p, s):
        num = 1 + s + p * s
        dnm = 1 + 2 * p * s
        return p * num / dnm
    def dominant_fwd(p, s):
        num = 1 + s
        dnm = 1 + 2 * p * s * (1 - p) + p * p * s
        return p * num / dnm
    def multiplicative_fwd(p, s):
        num = 1 + s
        dnm = 1 + p * s
        return p * num / dnm
    def recessive_fwd(p, s):
        num = 1 + p * s
        dnm = 1 + p * p * s
        return p * num / dnm

    # one step calculation
    if ploidy == 1:
        one_step = haploid_fwd
    else:
        if one_step_model == 'a':
            one_step = additive_fwd
        elif one_step_model == 'r':
            one_step = recessive_fwd
        elif one_step_model == 'd':
            one_step = dominant_fwd
        else:
            one_step = multiplicative_fwd

    # initialize
    ps = [] # frequencies
    xs = [] # variants
    Ns = [] # sizes
    t = floor(tau0)
    p = pG
    G = max(Ne.keys())
    N = Ne[G]
    x = ceil(p * ploidy * N)
    Ns.append(N)
    xs.append(x)
    ps.append(p)

    if random_walk: # random walk

        while G >= 0:
            G -= 1
            try:
                N = Ne[G]
            except KeyError:
                pass
            if G > t:
                p = one_step(p, s)
            x = int(binom.rvs(int(ploidy * N), p))
            p = x / ploidy / N
            if x < 1:
                break
            if p >= 1:
                break
            ps.append(p)
            xs.append(x)
            Ns.append(N)

        return np.array([np.array(ps), np.array(Ns), np.array(xs)]) # numpy-ify

    else: # deterministic

        while G >= 0:
            G -= 1
            try:
                N = Ne[G]
            except KeyError:
                pass
            if G > t:
                p = one_step(p, s)
            x = floor(p * ploidy * N)
            if x < 1:
                break
            if p >= 1:
                break
            Ns.append(N)
            xs.append(x)
            ps.append(p)

        return np.array([np.array(ps), np.array(Ns), np.array(xs)]) # numpy-ify

### simple coalescent simulators ###

def basic_coalescent(n):
    '''Simulate times in basic coalescent (scale post-hoc by population size)

    Parameters
    ----------
    n : int
        Sample size

    Returns
    -------
    numpy.array
        Interarrival times
    '''
    n = int(float(n))
    k = n
    t = np.zeros(n - 1)
    itr = 0
    curr = 0
    while k > 1:
        bcoef = 2 / k / (k - 1)
        k -= 1
        draw = Exp(bcoef)
        curr += draw
        t[itr] = curr
        itr += 1
    return t

def varying_Ne_coalescent(n, Ne, ploidy = 2, to_tmrca=True):
    '''Simulate times in varying population size coalescent

    Parameters
    ----------
    n : int
        Sample size
    Ne : dict
        Effective population sizes
    ploidy : int
    to_tmrca : bool
        Go to TMRCA

    Returns
    -------
    numpy.array
        Arrival times in generations
    '''

    # initialize
    Ne = zero_shift_Ne(Ne)
    maxG = max(Ne.keys())
    N = Ne[maxG]
    Nt = np.array(list(Ne.values()))
    lambdas = N / Nt # relative (inverse) coalescent rates
    loglambdas = np.log(lambdas)
    Lambdas = np.cumsum(lambdas) # piecewise integral
    Lambdas = np.insert(Lambdas, 0, 0)
    times = np.log(basic_coalescent(n)) + log(ploidy) + log(N) # convert to generations
    times = np.exp(times)
    times = np.insert(times, 0, 0)
    dtime = np.diff(times)
    K = len(dtime)
    Taus = np.zeros(K)
    k = 0
    itr = 0
    LambdaNext = 0
    vCurr = 0
    # while k <= (K-1) means until most recent common ancestor (MRCA)
    # itr < maxG controls updates with piecewise Ne
    while k <= (K-1) and itr < maxG:
        vPrev = vCurr
        LambdaPrev = LambdaCurr = LambdaNext
        RHS = dtime[k] + LambdaPrev
        LambdaNext = RHS
        steps = 0
        try:
            while LambdaCurr < RHS:
                itr += 1
                steps += 1
                LambdaCurr = Lambdas[itr]
            itr -= 1
            boole = int(steps>1)
            LambdaCurr = Lambdas[itr]*boole + LambdaPrev*(1-boole)
            if boole:
                vCurr += exp(log(LambdaCurr-LambdaPrev) - loglambdas[itr-1])
            RHS -= LambdaCurr
            logscalar = loglambdas[itr-1]
            logsoluti = log(RHS) - logscalar
            vCurr += exp(logsoluti)
            Taus[k] = (vCurr - vPrev)
            k += 1
        except IndexError:
            boole = int(steps>1)
            LambdaCurr = Lambdas[-1]*boole + LambdaPrev*(1-boole)
            RHS -= LambdaCurr
            logscalar = loglambdas[-1]
            while steps > 1:
                vCurr += 1
                steps -= 1
            logsoluti = log(RHS) - logscalar
            vCurr += exp(logsoluti)
            Taus[k] = (vCurr - vPrev)
            k += 1
    logscalar = loglambdas[-1]
    if to_tmrca:
        while k <= (K-1): # finish
            Taus[k] = exp(log(dtime[k]) - logscalar)
            k += 1
    Etas = np.cumsum(Taus)
    Etas = Etas[Etas > 0]
#     return np.floor(Etas)+1 # round to discrete time
    return Etas # continuous time

def binomial_step(k, Nei, ploidy=2):
    '''Take a Binomial(n,p) step for 1 gen of WF process

    Parameters
    ----------
    k : int
        Sample size
    Nei : float
        Effective size at time i
    ploidy : int

    Return
    ------
    int
        Count from binomial experiments
    '''
    denom = Nei * ploidy
    numer = int(k * (k - 1) / 2)
    bstep = binom.rvs(numer, 1/denom, size=1)[0]
    return bstep

def binomial_pmf_leq1(k, Nei, ploidy=2):
    '''Compute mass P(\leq 1) for Binomial(n,p) for 1 gen of WF process

    Parameters
    ----------
    k : int
        Sample size
    Nei : float
        Effective size at time i
    ploidy : int

    Return
    ------
    float
        A probability mass
    '''
    numer = int(k * (k - 1) / 2)
    denom = 1 / (Nei * ploidy)
    return binom.pmf(0, numer, denom) + binom.pmf(1, numer, denom)

def wright_fisher(n, Ne, ploidy = 2, to_tmrca=True):
    '''Simulate times in Wright Fisher model with varying population size
    (After last generation in Ne, assume constant size
    and apply scaled basic_coalescent)

    Parameters
    ----------
    n : int
        Sample size
    Ne : dict
        Effective population sizes
    ploidy : int
    to_tmrca : bool
        Continue to TMRCA

    Returns
    -------
    numpy.array
        Arrival times in generations
    '''
    # initalize
    n = int(float(n))
    k = n
    itr = 0
    gentimes = np.zeros(n-1)
    Me = deepcopy(Ne)
    timer = min(Me.keys())
    lastG = max(Me.keys())
    lastN = Me[lastG]
    cuml = 0
    while k > 1 and timer <= lastG:
        # record coalescent events
        size = Me[timer] * ploidy
        draw = [randint(0, size - 1) for i in range(k)]
        table = {}
        for i in range(len(draw)):
            try:
                table[draw[i]] += 1
            except KeyError:
                table[draw[i]] = 1
        events = [val for key, val in table.items() if val >= 2]
        l = sum(events) - len(events)
        cuml += l
        timer += 1
        gentimes[itr:cuml] = timer
        k -= l
        itr += l
    if k > 1: # finish with coalescent
        if to_tmrca:
            finish = basic_coalescent(k) * lastN + lastG
            finish = np.floor(finish) + 1 # round to discrete time
            gentimes[cuml:] = finish
    return gentimes[gentimes > 0]

def simulate_coalgen(n, Ne, ploidy, binom_cut=0.01, to_tmrca=True):
    '''Simulate times in varying size population with large sample.
    Use binomial walk, Wright Fisher process, and the coalescent.

    Parameters
    ----------
    n : int
        Sample size
    Ne : dict
        Effective population sizes
    ploidy : int
    to_tmrca : bool
        Go to TMRCA
    binom_cut : float
        Threshold to conduct binomial update [P(X <= 1) <= \alpha]

    Returns
    -------
    numpy.array
        Arrival times in generations
    '''

    assert n > 1
    assert binom_cut > 0
    assert binom_cut < 1
    n = int(float(n))
    j = n
    Taus = np.zeros(n-1)
    itr = min(Ne.keys()) + 1
    k = 0
    try: # while loop ends, i.e. condition to switch from binomial walk achieved

        # binomial walk process
        while binomial_pmf_leq1(j, Ne[itr], ploidy) <= binom_cut:
            # i is binomial count
            # j is nodes in wf process
            # k is index iterating up
            i = binomial_step(j, Ne[itr], ploidy)
            j -= i
            if j > 1:
                itr += 1
                Taus[k:(k+i)] = itr
                k += i
            else:
                itr += 1
                Taus[k:] = itr
        Taus = Taus[Taus>0]

        # wright fisher / coalescent process
        # wright fisher runs a basic coalescent at end of Ne dictionary
        sNe = {k:v for k,v in Ne.items() if k >= itr}
        Etas = wright_fisher(j, sNe, ploidy, to_tmrca)
        Taus = np.concatenate((Taus,Etas))

    except KeyError: # while loop does not end, i.e. reach end of defined Ne
        Taus = Taus[Taus>0]

    return Taus

### coalescent network simulator ###

import numpy as np
import pandas as pd
import networkx as nx
from scipy.stats import binom
from random import randint, sample
from math import floor, exp, log, ceil
from numpy.random import exponential as Exp

### coalescent network simulator ###

class DirectedPaths:

    def __init__(self, _id, _pairwise_output):
        self._left_length = np.inf
        self._right_length = np.inf
        self._num_leaves = 1
        self._leading_leaf = _id
        self._pairwise_output = _pairwise_output
        if _pairwise_output:
            self._leaves = [_id]
        return None

    def add_edge(self, _left_length, _right_length):
        self._left_length = min(self._left_length, _left_length)
        self._right_length = min(self._right_length, _right_length)
        return None

class Node:

    def __init__(self, _id, _time, _pairwise_output):
        self._id = _id
        self._time = _time
        if _time <= 0:
            self._paths_list = [DirectedPaths(_id, _pairwise_output)]
            self._num_leaves = 1
        else:
            self._paths_list = []
            self._num_leaves = 0
        self._num_tracts = 0
        return None
    def __repr__(self):
        return 'Node(_id = ' + str(self._id) + ')'
    def __str__(self):
        return 'Node(_id = ' + str(self._id) + ')'

    def add_incoming_node(self, node):
        self._num_leaves += node._num_leaves
        self._num_tracts += node._num_tracts
        return None

    def add_incoming_edge(self, scale, scales=True):
        if scales: # for multiple compare
            _left_length = Exp(scale=scale, size=1)[0]
            _right_length = Exp(scale=scale, size=1)[0]
            for paths in self._paths_list:
                paths.add_edge(_left_length, _right_length)
        else: # for pair compare
            if scale != 0:
                scale = 1 / scale
                _left_length = Exp(scale=scale, size=1)[0]
                _right_length = Exp(scale=scale, size=1)[0]
                for paths in self._paths_list:
                    paths.add_edge(_left_length, _right_length)
            else:
                pass
        return None

    def merge(self):
        paths_list = self._paths_list
        merged = set() # merge set
        doskip = set() # skip set
        k = len(paths_list)
        for i in range(k):
            if i not in doskip:
                merged.add(i)
                for j in range(i+1,k):
                    if j not in doskip:
                        paths1 = paths_list[i]
                        paths2 = paths_list[j]
                        if paths1._left_length == paths2._left_length:
                            if paths1._right_length == paths2._right_length:
                                m1 = paths1._num_leaves
                                m2 = paths2._num_leaves
                                m = m1 + m2
                                paths_list[i]._num_leaves = m
                                paths_list[i]._leading_leaf = min(paths1._leading_leaf, paths2._leading_leaf)
                                if paths_list[i]._pairwise_output:
                                    paths_list[i]._leaves += paths2._leaves
                                doskip.add(j)
        # merging
        self._paths_list = [paths_list[i] for i in merged]
        return None

    def prune_paths(self, cutoff):
        self._paths_list = [paths for paths in self._paths_list if paths._left_length + paths._right_length >= cutoff]
        return None

    def pair_coalesce(self, node1, node2, cutoff1, cutoff2, subgraph_id=2, record_dist=True, pairwise_output=True):
        deltas = [self._time - node1._time, self._time - node2._time]
        rates = [delta / 100 for delta in deltas]
        node1.add_incoming_edge(rates[0], scales=False)
        node1.prune_paths(cutoff2)
        node1.merge()
        node2.add_incoming_edge(rates[1], scales=False)
        node2.prune_paths(cutoff2)
        node2.merge()
        _num_tracts, n1, n2 = pairwise_compare(node1, node2, cutoff1, cutoff2, self._time, subgraph_id, record_dist, pairwise_output)
        node1 = n1
        node2 = n2
        self._num_tracts += _num_tracts
        self.add_incoming_node(node1)
        for paths in node1._paths_list:
            self._paths_list.append(paths)
        self.add_incoming_node(node2)
        for paths in node2._paths_list:
            self._paths_list.append(paths)
        return None

def pairwise_compare(node1, node2, cutoff1, cutoff2, _time, subgraph_id=2, record_dist=True, pairwise_output=True):
    '''Compare paths up coalescent tree based in IBD

    Parameters
    ----------
    cutoff1 : float
        Long IBD cM threshold
    cutoff2 : float
        Long IBD cM threshold
    _time : int
        Generation time

    Returns
    -------
    tuple
        (int # of long IBD segments, Node class instance, Node class instance)
    '''
    # initializing
    global H1, H0, H, ldist1, ldist0, ldist, tdist1, tdist0, tdist, ddist1, ddist0, ddist
    global pairwise_segments
    _num_tracts = 0
    K = len(node1._paths_list)
    L = len(node2._paths_list)
    for k in range(K):
        i = node1._paths_list[k]
        for l in range(L):
            j = node2._paths_list[l]
            _left_length = min(i._left_length, j._left_length)
            _right_length = min(i._right_length, j._right_length)
            _length = _right_length + _left_length
            # networking
            if _length >= cutoff2:
                if subgraph_id == 1:
                    H1.add_edge(i._leading_leaf, j._leading_leaf)
                elif subgraph_id == 0:
                    H0.add_edge(i._leading_leaf, j._leading_leaf)
                else:
                    H.add_edge(i._leading_leaf, j._leading_leaf)
                if i._leading_leaf <= j._leading_leaf:
                    node2._paths_list[l]._leading_leaf = node1._paths_list[k]._leading_leaf
                else:
                    node1._paths_list[k]._leading_leaf = node2._paths_list[l]._leading_leaf
                # segment counting
                if _length >= cutoff1:
                    m1 = i._num_leaves
                    m2 = j._num_leaves
                    density = m1 * m2
                    _num_tracts += density
                    if record_dist:
                        if subgraph_id == 1:
                            ldist1.append(_length)
                            tdist1.append(_time)
                            ddist1.append(density)
                        elif subgraph_id == 0:
                            ldist0.append(_length)
                            tdist0.append(_time)
                            ddist0.append(density)
                        else:
                            ldist.append(_length)
                            tdist.append(_time)
                            ddist.append(density)
                    if i._pairwise_output:
                        for a in i._leaves:
                            for b in j._leaves:
                                pairwise_segments.append((a,b,_length,_time))

    return _num_tracts, node1, node2

def two_randint(n):
    x = randint(0,n)
    y = randint(0,n)
    while x == y:
        y = randint(0,n)
    return sorted([x,y])

def simulate_ibd(n, Ne, long_ibd=2.0, short_ibd=1.0, ploidy=2, record_dist=True, pairwise_output=True):
    '''ibd segments from a coalescent

    Parameters
    ----------
    n : int
        Sample size (individuals)
    Ne : dict
        Effective population sizes
    long_ibd: float
        cM length threshold
    short_ibd : float
        cM length threshold
    ploidy : int
        1 for haploid or 2 for diploid
    record_dist : bool
        To save tract length and coalescent time distributions or not to (default True)
    pairwise_output : bool
        To save pairwise segments or not to (default True)

    Returns
    -------
    tuple
        (number of tracts, group sizes, length distr., time distr., count distr., pairwise segments)
    '''

    # checks
    assert long_ibd > 0
    assert short_ibd > 0
    assert long_ibd >= short_ibd
    n = int(float(n))
    ploidy = int(float(ploidy))

    # initialize graph network
    global H, H1, H0, ldist, ldist1, ldist0, tdist1, tdist0, tdist, ddist1, ddist0, ddist
    global pairwise_segments
    pairwise_segments = []
    H = nx.Graph()
    ldist = []
    ddist = []
    tdist = []

    # initalize
    m = ploidy*n
    interiors = [Node(i, 0, pairwise_output) for i in range(m)]
    itr = m

    # simulate times
    Ne=cut_Ne(to_max_Ne(fill_Ne(Ne),500),2000)
    times = simulate_coalgen(m, Ne, ploidy, to_tmrca=True)

    indxs = [i for i in range(m)]
    for t in times:
        m -= 1
        new_node = Node(itr, t, pairwise_output)
        sindx = sorted(two_randint(m))
        s1 = sindx[0]
        s2 = sindx[1]
        new_node.pair_coalesce(interiors[s1], interiors[s2], long_ibd, short_ibd, 2, record_dist)
        interiors.append(new_node)
        interiors.pop(s1)
        interiors.pop(s2-1)
        indxs.pop(s1)
        indxs.pop(s2-1)
        indxs.append(itr)
        itr += 1

    numTracts = sum([node._num_tracts for node in interiors])
    ibdGroupSizes = sorted([len(c) for c in nx.connected_components(H)],reverse=True)

    # tract count, group sizes, distributions: lengths, times, density
    return (numTracts,
            np.array(ibdGroupSizes),
            np.array(ldist),
            np.array(tdist),
            np.array(ddist),
            tuple(pairwise_segments)
           )

def simulate_ibd_isweep(n, s, p0, Ne, long_ibd=2.0, short_ibd=1.0, random_walk=True, one_step_model='a', tau0=0, sv=-0.01, ploidy=2, record_dist=True, pairwise_output=True):
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
    sv: float
        Allele frequency of standing variation
        (Default -0.01 will assume de novo sweep)
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
    assert sv < 1

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
                           long_ibd, 
                           short_ibd,
                           ploidy,
                           record_dist, 
                           pairwise_output
                          )
        # numpy-ify
        return (out,
                (np.nan,np.array([]),np.array([]),np.array([]),np.array([])),
                (np.nan,np.array([]),np.array([]),np.array([]),np.array([])),
                tuple(pairwise_segments)
               )

    # calculating structured demographies
    ps, Ns, xs = walk_variant_backward(s, p, Ne, stoc, mdl, tau0, sv, ploidy)
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

### segments ###

# quasi geometric random variables

def probability_quasi_geometric(ps, Ns, ploidy = 2):
    """Probability masses for quasi-geometric coalescent

    Parameters
    ----------
    ps : numpy.array
        Variant frequencies (back in time)
    Ns : numpy.array
        Effective population sizes
    ploidy : int
        1 for haploid or 2 for diploid

    Returns
    -------
    numpy.array
        P(G = g) for generation g; zero element is P(G > g) for some large g
    """

    assert ploidy in [1,2]
    ones = np.array([1 for i in range(len(ps))])
    marg_mass = 1 / (ploidy * ps[:-1] * Ns[:-1])
    marg_mass = np.insert(marg_mass, 0, 0)
    marg_mass = np.minimum(ones, marg_mass)
    surv_mass = ones - marg_mass
    surv_mass = np.cumprod(surv_mass)
    marg_mass = 1 / (ploidy * ps * Ns)
    marg_mass = np.minimum(ones, marg_mass)
    full_mass = surv_mass * marg_mass
    full_mass = np.insert(full_mass, 0, np.maximum(1 - np.sum(full_mass), 0))

    return full_mass

def simulate_quasi_geometric(n, mass):
    """Simulator for quasi-geometric coalescent

    Parameters
    ----------
    n : int
        Simulation size
    mass : numpy.array
        Probability masses

    Returns
    -------
    numpy.array
        Coalescent generation times; np.inf for P(G > g) for some large g
    """

    n = int(float(n))
    time = np.array([np.inf] + [i for i in range(1, len(mass))])
    coal = np.random.choice(time, size = n, p = mass)

    return coal

# erlang random variables

def simulate_erlang_segments(geom):
    """Simulator for Erlang (shape=2) ibd segment lengths

    Parameters
    ----------
    geom : numpy.array
        Coalescent generation times (quasi-geometric)

    Returns
    -------
    numpy.array
        Independent ibd segment lengths
    """

    geom = geom[geom != np.inf]
    m = len(geom)
    left = np.random.exponential(scale = 50 / geom, size = m)
    right = np.random.exponential(scale = 50 / geom, size = m)

    return left + right

def probability_erlang_segments(ab, mass):
    '''Interval probabilities for Erlang (shape=2) ibd lengths

    Parameters
    ----------
    ab : numpy.array
        Increasing floats in centiMorgans
    mass : numpy.array
        P(G = g) for generation g

    Returns
    -------
    numpy.array
        P(\ell \in [u,v)) for many [u,v) intervals
    '''

    # local function
    def puv(u, v, g):

        # error handling
        if u < 0:
            raise ValueError("Left endpoint less than zero")
        if u >= v:
            raise ValueError("Left endpoint exceeds right")

        # integral values
        U = np.exp(- g / 50 * u) * (u * g + 50)
        if np.isinf(v):
            V = 0
        else:
            V = np.exp(- g / 50 * v) * (v * g + 50)

        return (U - V) / 50

    # initialize
    mass = mass[1:]
    geom = np.array([i for i in range(1, len(mass) + 1)])
    M = len(ab)

    # double integral
    puvs = []
    for m in range(1, M):
        pab = puv(ab[m-1], ab[m], geom)
        puvs.append(sum(mass * pab))

    return np.array(puvs)

# probability of ibd

def probability_ibd(ps, Ns, long_ibd = 2, ploidy = 2):
    """Approximate probability of ibd

    Parameters
    ----------
    ps : numpy.array
        Variant frequencies (back in time)
    Ns : numpy.array
        Effective population sizes
    long_ibd : float
        cM length threshold
    ploidy : int
        1 for haploid or 2 for diploid

    Returns
    -------
    float
        approx P(\ell > c) where \ell is ibd length
    """

    assert ploidy in [1,2]

    # initialize
    p0 = ps[0]

    # cumulative products
    with np.errstate(divide = 'ignore'):
        dnm1 = 1 / (ploidy * ps[:-1] * Ns[:-1])
    dnm1 = np.insert(dnm1, 0, 0)
    one = np.array([1 for i in range(len(dnm1))])
    sub1 = np.minimum(one, dnm1)
    onesub1 = one - sub1
    cp1 = np.cumprod(onesub1)

    # products
    with np.errstate(divide = 'ignore'):
        dnm1 = 1 / (ploidy * ps * Ns)
    dnm1 = np.minimum(one, dnm1)

    # partial sum
    pgs = (p0 ** 2) * cp1 * dnm1
    gstar = len(pgs)

    masses = []
    for g in range(1, gstar + 1):
        masses.append((long_ibd * g / 50 + 1) * exp(- g * long_ibd / 50))
    masses = np.array(masses)
    partialp = sum(pgs * masses)

    # approximate coalescent
    if dnm1[-1] >= 0: # check this
        FRp = - np.inf
    else:
        hp = (p0 ** 2) * cp1[-1] * (1 - dnm1[-1])
        dp = 1 / ploidy / Ns[-1] / ps[-1]
        Cp = log(1 - dp) if not np.isinf(dp) else 0
        CCp = Cp - long_ibd / 50
        FRp = log(-(CCp + long_ibd / 50 * ((gstar + 1) * CCp - 1)))
        FRp = FRp + CCp * (gstar + 1) - log(CCp ** 2)
        FRp = FRp + log(hp) + log(dp) - Cp * (gstar + 1)

    return partialp + exp(FRp)

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

### for asymptotic tests ###

def simulate_ibd_constant(n, Ne, long_ibd=2.0, short_ibd=1.0, ploidy=2, record_dist=True, pairwise_output=True):
    '''ibd segments from a coalescent w/ constant Ne

    Parameters
    ----------
    n : int
        Sample size (individuals)
    Ne : int
        Constant effective population size
    long_ibd, short_ibd : float
        cM length threshold
    ploidy : int
        1 for haploid or 2 for diploid
    record_dist : bool
        To save tract length and coalescent time distributions or not to (default True)
    pairwise_output : bool
        To save pairwise segments or not to (default True)

    Returns
    -------
    tuple
        (number of tracts, group sizes, length distr., time distr., count distr., pairwise segments)
    '''

    # checks
    assert long_ibd > 0
    assert short_ibd > 0
    assert long_ibd >= short_ibd
    n = int(float(n))
    ploidy = int(float(ploidy))

    # initialize graph network
    global H, H1, H0, ldist, ldist1, ldist0, tdist1, tdist0, tdist, ddist1, ddist0, ddist
    global pairwise_segments
    pairwise_segments = []
    H = nx.Graph()
    ldist = []
    ddist = []
    tdist = []

    # initalize
    m = ploidy*n
    interiors = [Node(i, 0, pairwise_output) for i in range(m)]
    itr = m

    # this is the simple change
    # do basic coalescent instead of wright-fisher
    times = basic_coalescent(m) * float(int(Ne)*ploidy)

    indxs = [i for i in range(m)]
    for t in times:
        m -= 1
        new_node = Node(itr, t, pairwise_output)
        sindx = sorted(two_randint(m))
        s1 = sindx[0]
        s2 = sindx[1]
        new_node.pair_coalesce(interiors[s1], interiors[s2], long_ibd, short_ibd, 2, record_dist)
        interiors.append(new_node)
        interiors.pop(s1)
        interiors.pop(s2-1)
        indxs.pop(s1)
        indxs.pop(s2-1)
        indxs.append(itr)
        itr += 1

    numTracts = sum([node._num_tracts for node in interiors])
    ibdGroupSizes = sorted([len(c) for c in nx.connected_components(H)],reverse=True)

    # tract count, group sizes, distributions: lengths, times, density
    return (numTracts,
            np.array(ibdGroupSizes),
            np.array(ldist),
            np.array(tdist),
            np.array(ddist),
            tuple(pairwise_segments)
           )


### time-varying selection ###

def walk_variant_forward_tv(s_s, g_s, pG, Ne, random_walk = False, ploidy = 2):
    '''Variant frequencies forward in time

    Parameters
    ----------
    s_s : list
        List of selection coefficients
    g_s : list
        List of transition times for selection coefficients
        This should be aligned with the s_s parameter
            e.g., s_s=[0.01,0.03,0.02] and g_s=[50,100] means 
            s=0.02 between 100-, s=0.03 between 100-50 and s=0.01 between 50-0
        Should be length one less than s_s
    pG : float
        Variant frequency at maximum generation
    Ne : dict
        Effective population sizes
    random_walk : bool
        True for random walk
    ploidy : int
        1 for haploid or 2 for diploid

    Additive genic selection is assumed

    Returns
    -------
    tuple
        NumPy arrays for frequencies and sizes
    '''

    # local functions
    assert ploidy in [1,2]
    def haploid_fwd(p, s): # haploid is same as multiplicative
        return p * (1 + s) / (1 + p * s)
    def additive_fwd(p, s):
        num = 1 + s + p * s
        dnm = 1 + 2 * p * s
        return p * num / dnm

    # one step calculation
    if ploidy == 1:
        one_step = haploid_fwd
    else:
        one_step = additive_fwd

    # initialize
    ps = [] # frequencies
    xs = [] # variants
    Ns = [] # sizes
    p = pG
    G = max(Ne.keys())
    assert len(s_s) == (len(g_s) + 1)
    s_s2 = deepcopy(s_s)
    g_s2 = deepcopy(g_s)
    s2 = s_s2[-1]
    s_s2 = s_s2[:-1]
    sdict = {g_s2[i]:s_s2[i] for i in range(len(g_s))}
    N = Ne[G]
    x = ceil(p * ploidy * N)
    Ns.append(N)
    xs.append(x)
    ps.append(p)

    if random_walk: # random walk

        while G >= 0:
            G -= 1
            try:
                N = Ne[G]
            except KeyError:
                pass
            try:
                s2 = sdict[G]
            except KeyError:
                pass
            p = one_step(p, s2)
            x = int(binom.rvs(int(ploidy * N), p))
            p = x / ploidy / N
            if x < 1:
                break
            if p >= 1:
                break
            ps.append(p)
            xs.append(x)
            Ns.append(N)

        return np.array([np.array(ps), np.array(Ns), np.array(xs)]) # numpy-ify

    else: # deterministic

        while G >= 0:
            G -= 1
            try:
                N = Ne[G]
            except KeyError:
                pass
            try:
                s2 = sdict[G]
            except KeyError:
                pass
            p = one_step(p, s2)
            x = floor(p * ploidy * N)
            if x < 1:
                break
            if p >= 1:
                break
            Ns.append(N)
            xs.append(x)
            ps.append(p)

        return np.array([np.array(ps), np.array(Ns), np.array(xs)]) # numpy-ify

def walk_variant_backward_tv(s_s, g_s, p0, Ne, random_walk = False, ploidy = 2):
    '''Variant frequencies backward in time

    Parameters
    ----------
    s_s : list
        List of selection coefficients
    g_s : list
        List of transition times for selection coefficients
        This should be aligned with the s_s parameter
            e.g., s_s=[0.01,0.03,0.02] and g_s=[50,100] means 
            s=0.01 between 0-50, s=0.03 between 50-100 and s=0.02 between 100-
        Should be length one less than s_s
    p0 : float
        Variant frequency at generation 0
    Ne : dict
        Effective population sizes
    random_walk : bool
        True for random walk
    ploidy : int
        1 for haploid or 2 for diploid


    Returns
    -------
    tuple
        NumPy arrays for frequencies and sizes
    '''

    # local functions
    assert ploidy in [1,2]
    assert p0 <= 1
    assert p0 >= 0
    def haploid_bwd(p, s): # haploid is same as multiplicative (Felsenstein, 2017)
        return p / (1 + s - s * p)
    def additive_bwd(p, s):
        if s <= 0:
            return p
        a = s
        b = 1 + s - 2 * p * s
        c = - p
        qf = - b + sqrt((b ** 2) - 4 * a * c)
        qf = qf / 2 / a
        return qf

    # one step calculation
    if ploidy == 1:
        one_step = haploid_bwd
    else:
        one_step = additive_bwd

    # initialize
    ps = [] # frequencies
    xs = [] # variants
    Ns = [] # sizes
    p = p0
    N = Ne[0]
    assert len(s_s) == (len(g_s)+1)
    s_s2 = deepcopy(s_s)
    g_s2 = deepcopy(g_s) 
    s2 = deepcopy(s_s2[0])
    s_s2 = s_s2[1:]
    sdict = {g_s2[i]:s_s2[i] for i in range(len(g_s))}
    x = floor(p * ploidy * N)
    Ns.append(N)
    xs.append(x)
    ps.append(p)

    if random_walk: # random walk

        for G in range(1, max(Ne.keys())+1):
            try: # population size change
                N = Ne[G]
            except KeyError:
                pass
            try:
                s2 = sdict[G]
            except KeyError:
                pass
            p = one_step(p, s2)
            x = int(binom.rvs(int(ploidy * N), p))
            p = x / ploidy / N
            if x < 1:
                break
            if p >= 1:
                break
            ps.append(p)
            xs.append(x)
            Ns.append(N)

        return np.array([ps, Ns, xs], dtype=float) # numpy-ify

    else: # deterministic

        for G in range(1, max(Ne.keys())+1):
            try: # population size change
                N = Ne[G]
            except KeyError:
                pass
            try:
                s2 = sdict[G]
            except KeyError:
                pass
            p = one_step(p, s2)
            x = floor(p * ploidy * N)
            if x < 1:
                break
            if p >= 1:
                break
            Ns.append(N)
            xs.append(x)
            ps.append(p)

        return np.array([np.array(ps), np.array(Ns), np.array(xs)]) # numpy-ify
    
def simulate_ibd_isweep_tv(n, 
                           s_s, 
                           g_s, 
                           p0, 
                           Ne, 
                           long_ibd=2.0, 
                           short_ibd=1.0, 
                           random_walk=True, 
                           ploidy=2, 
                           record_dist=True, 
                           pairwise_output=True
                           ):
    '''ibd segments from a coalescent with selection

    Parameters
    ----------
    n : int
        Sample size (individuals)
    s_s : list
        List of selection coefficients
    g_s : list
        List of transition times for selection coefficients
        This should be aligned with the s_s parameter
            e.g., s_s=[0.01,0.03,0.02] and g_s=[50,100] means 
            s=0.01 between 0-50, s=0.03 between 50-100 and s=0.02 between 100-
        Should be length one less than s_s
    p0 : float
        Variant frequency at generation 0
    Ne : dict
        Effective population sizes
    long_ibd, short_ibd : float
        cM length threshold
    random_walk : bool
        True for random walk
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
    stoc = random_walk
    Ne = cut_Ne(to_max_Ne(fill_Ne(Ne),500),2000)

    # should p0 be fixed
    if (p0 == 0) or (p0 == 1):
        out = simulate_ibd(n, Ne,
                           long_ibd, 
                           short_ibd,
                           ploidy,
                           record_dist, 
                           pairwise_output
                          )
        # numpy-ify
        return (out,
                (np.nan,np.array([]),np.array([]),np.array([]),np.array([])),
                (np.nan,np.array([]),np.array([]),np.array([]),np.array([])),
                tuple(pairwise_segments)
               )

    # calculating structured demographies
    ps, Ns, xs = walk_variant_backward_tv(s_s, g_s, p, Ne, stoc, ploidy)
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

