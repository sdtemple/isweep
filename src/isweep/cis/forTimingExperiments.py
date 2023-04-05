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
        self._left_length = np.Inf
        self._right_length = np.Inf
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
    
    def pair_coalesce(self, node1, node2, cutoff1, cutoff2, subgraph_id=2, record_dist=True, execute_merge=True, execute_prune=True):
        deltas = [self._time - node1._time, self._time - node2._time]
        rates = [delta / 100 for delta in deltas]
        node1.add_incoming_edge(rates[0], scales=False)
        if execute_prune:
            node1.prune_paths(cutoff2)
        if execute_merge:
            node1.merge()
        node2.add_incoming_edge(rates[1], scales=False)
        if execute_prune:
            node2.prune_paths(cutoff2)
        if execute_merge:
            node2.merge()
        _num_tracts, n1, n2 = pairwise_compare(node1, node2, cutoff1, cutoff2, self._time, subgraph_id, record_dist)
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
    
def pairwise_compare(node1, node2, cutoff1, cutoff2, _time, subgraph_id=2, record_dist=True, pairwise_output=False):
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

def simulate_ibd(n, Ne, long_ibd=2.0, short_ibd=1.0, ploidy=2, record_dist=False, pairwise_output=False, execute_merge=True, execute_prune=True):
    '''ibd segments from a coalescent
    
    Parameters
    ----------
    n : int
        Sample size (individuals)
    Ne : dict
        Effective population sizes
    long_ibd, short_ibd : float
        cM length threshold
    ploidy : int
        1 for haploid or 2 for diploid
    continuous_time : bool
        To implement coalescent or Wright-Fisher (default True is coalescent)
    record_dist : bool
        To save tract length and coalescent time distributions or not to (default False)
    pairwise_output : bool
        To save pairwise segments or not to (default False)
    execute_merge : bool
        To implement algorithm with merging or not to (default True)
    execute_prune : bool
        To implement algorithm with pruning or not to (default True)
        
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
    global H, ldist, tdist, ddist
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
    times = simulate_coalescent_generations(m, Ne, ploidy, to_tmrca=True)
        
    indxs = [i for i in range(m)]
    for t in times:
        m -= 1
        new_node = Node(itr, t, pairwise_output)
        sindx = sorted(two_randint(m))
        s1 = sindx[0]
        s2 = sindx[1]
        new_node.pair_coalesce(interiors[s1], interiors[s2], long_ibd, short_ibd, 2, record_dist, execute_merge, execute_prune)
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

def simulate_ibd_from_selective_sweep(n, s, p0, Ne, long_ibd=2.0, short_ibd=1.0, random_walk=True, one_step_model='m', tau0=0, ploidy=2, record_dist=False, pairwise_output=False, execute_merge=True, execute_prune=True):
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
    continuous_time : bool
        To implement coalescent or Wright-Fisher (default True is coalescent)
    record_dist : bool
        To save tract length and coalescent time distributions or not to (default False)
    pairwise_output : bool
        To save pairwise segments or not to (default False)
    execute_merge : bool
        To implement algorithm with merging or not to (default True)
    execute_prune : bool
        To implement algorithm with pruning or not to (default True)
    
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
    Ne = cut_Ne(to_max_Ne(fill_Ne(Ne),500),1000)
    
    # should p0 be fixed
    if (p0 == 0) or (p0 == 1):
        out = simulate_ibd(n, Ne, 
                           long_ibd, short_ibd, 
                           ploidy, continuous_time, 
                           record_dist, pairwise_output,
                           execute_merge, execute_prune
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
    times1 = simulate_coalescent_generations(n1, N1, ploidy, to_tmrca=False)
    
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
        new_node.pair_coalesce(interiors1[s1], interiors1[s2], long_ibd, short_ibd, 1, record_dist, execute_merge, execute_prune)
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
    N01 = {key:val for key,val in N0.items() if key <= maxt1}
    N00 = {key:val for key,val in N0.items() if key > maxt1}
    times0 = simulate_coalescent_generations(n0, N01, ploidy, to_tmrca=False)
    m = n0
    indxs = [i for i in range(m)]
    for t in times0:
        m -= 1
        new_node = Node(itr, t, pairwise_output)
        sindx = sorted(two_randint(m))
        s1 = sindx[0]
        s2 = sindx[1]
        new_node.pair_coalesce(interiors0[s1], interiors0[s2], long_ibd, short_ibd, 0, record_dist, execute_merge, execute_prune)
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
    if len(N00.keys()) > 0:
        times2 = simulate_coalescent_generations(m, N00, ploidy, to_tmrca=True)
        for t in times2:
            m -= 1
            new_node = Node(itr, t, pairwise_output)
            sindx = sorted(two_randint(m))
            s1 = sindx[0]
            s2 = sindx[1]
            new_node.pair_coalesce(interiors[s1], interiors[s2], long_ibd, short_ibd, 2, record_dist, execute_merge, execute_prune)
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
