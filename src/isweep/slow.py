# coalescent ibd simulator for paper
# slow multiple coalesce, wf process based on
# one with on/off buttons for merge, prune
# these are much slower, hence the name slow.py

# importing relevant packages
import numpy as np
import pandas as pd
import networkx as nx
from scipy.stats import binom
from random import randint, sample
from math import floor, exp, log, ceil
from numpy.random import exponential as Exp

from .coalescent import *

def empty_function():
    return None

def two_randint(n):
    x = randint(0,n)
    y = randint(0,n)
    while x == y:
        y = randint(0,n)
    return sorted([x,y])

### multiple coalesce, wf process only ###

class WfDirectedPaths:

    def __init__(self, _id):
        self._left_length = np.inf
        self._right_length = np.inf
        self._num_leaves = 1
        self._leading_leaf = _id
        return None

    def add_edge(self, _left_length, _right_length):
        self._left_length = min(self._left_length, _left_length)
        self._right_length = min(self._right_length, _right_length)
        return None

class WfNode:

    def __init__(self, _id, _time):
        self._id = _id
        self._time = _time
        if _time <= 0:
            self._paths_list = [WfDirectedPaths(_id)]
            self._num_leaves = 1
        else:
            self._paths_list = []
            self._num_leaves = 0
        self._num_tracts = 0
        return None
    def __repr__(self):
        return 'WfNode(_id = ' + str(self._id) + ')'
    def __str__(self):
        return 'WfNode(_id = ' + str(self._id) + ')'

    def add_incoming_node(self, node):
        self._num_leaves += node._num_leaves
        self._num_tracts += node._num_tracts
        return None

    def add_incoming_edge(self, scale):
        _left_length = Exp(scale=scale, size=1)[0]
        _right_length = Exp(scale=scale, size=1)[0]
        for paths in self._paths_list:
            paths.add_edge(_left_length, _right_length)
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
                                doskip.add(j)
        # merging
        self._paths_list = [paths_list[i] for i in merged]
        return None

    def prune_paths(self, cutoff):
        self._paths_list = [paths for paths in self._paths_list if paths._left_length + paths._right_length >= cutoff]
        return None

    def multiple_coalesce(self, nodes, cutoff1, cutoff2):
        # initializing
        nodes = [node for node in nodes]
        deltas = [self._time - node._time for node in nodes]
        scales = [100 / delta for delta in deltas]
        # updating node info
        for i in range(len(nodes)):
            nodes[i].add_incoming_edge(scales[i])
            nodes[i].prune_paths(cutoff2)
            nodes[i].merge()
        # pairwise comparing
        for i in range(len(nodes)-1):
            for j in range(i+1,len(nodes)):
                _num_tracts, node1, node2 = wf_pairwise_compare(nodes[i], nodes[j], cutoff1, cutoff2)
                nodes[i] = node1
                nodes[j] = node2
                self._num_tracts += _num_tracts
        # passing node info up
        for i in range(len(nodes)):
            self.add_incoming_node(nodes[i])
            for paths in nodes[i]._paths_list:
                self._paths_list.append(paths)
        return None

def wf_pairwise_compare(node1, node2, cutoff1, cutoff2):
    # initializing
    global H
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
    return _num_tracts, node1, node2

def simulate_ibd_wf(n, Ne, ibd_min_span=2.0, ibs_min_span=1.0):
    '''ibd segments from a wright fisher process (slower but exact)

    Parameters
    ----------
    n : int
        Haploid sample size
    Ne : dict
        Effective population sizes
    ibd_min_span : float
        cM length threshold
    ibs_min_span : float
        cM length threshold

    Returns
    -------
    list
        Objects of Node class
    '''

    # initialize graph network
    global H
    H = nx.Graph()

    # initalize
    n = int(float(n))
    interiors = [WfNode(i, 0) for i in range(n)]
    itr = n
    timer = 1
    end = max(Ne.keys())
    while timer <= end:
        # record coalescent events
        size = Ne[timer]
        draw = [randint(0, size - 1) for i in range(len(interiors))]
        table = {}
        for i in range(len(draw)):
            try:
                table[draw[i]].append(i)
            except KeyError:
                table[draw[i]] = [i]
        events = [val for key, val in table.items() if len(val) >= 2]
        # record pairwise ibd
        drop = []
        for merge in events:
            # incoporate incoming node data into next level interior node
            itr += 1
            new_node = WfNode(itr, timer)
            nodes = []
            for i in merge:
                nodes.append(interiors[i])
                drop.append(i)
            new_node.multiple_coalesce(nodes, ibd_min_span, ibs_min_span)
            # multiple coalescent, often 2
            interiors.append(new_node)
        # dropping old nodes
        drop = sorted(drop)
        popr = 0
        for d in drop:
            interiors.pop(d + popr)
            popr -= 1
        timer += 1
    numTracts = sum([node._num_tracts for node in interiors])
    maxIbdGroup = max([len(c) for c in nx.connected_components(H)])

    return numTracts, maxIbdGroup

def simulate_ibd_isweep_wf(n, s, p0, Ne, ibd_min_span = 2, ibs_min_span = 1.0, random = True, model = 'm', tau0 = 0, ploidy = 2):
    '''ibd segments from a wright fisher process with selection (slower but exact)

    Parameters
    ----------
    n : int
        Haploid sample size
    s : float
        Selection coefficient
    p0 : float
        Variant frequency at generation 0
    Ne : dict
        Effective population sizes
    ibd_min_span : float
        cM length threshold
    ibs_min_span : float
        cM length threshold
    random : bool
        True for random walk
    model : str
        'm', 'a', 'd', or 'r'
    tau0 : int
        Generation when neutrality begins
    ploidy : int
        Assume diploid

    Returns
    -------
    int
        Count of ibd segments
    '''

    # local function
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

    global H
    H = nx.Graph()

    # renaming
    p = p0
    t = tau0
    mdl = model
    stoc = random

    # calculating structured demographies
    ps, Ns, xs = walk_variant_backward(s, p, Ne, stoc, mdl, tau0, ploidy)
    n = int(float(n))
    n = n * ploidy
    G = max(Ne.keys())
    n1 = floor(n * ps[0])
    n0 = n - n1
    N1 = [ceil(i) for i in ploidy * Ns * ps]
    N1 = {i:N1[i] for i in range(len(N1))}
    qs = extend_vector(1 - ps, 1, G)
    Ns = np.array(list(Ne.values()))
    N0 = [ceil(i) for i in ploidy * Ns * qs]
    N0 = {i:N0[i] for i in range(len(N0))}

    # population 1
    timer = 1
    end = max(N1.keys())
    itr = n1
    interiors1 = [WfNode(i, 0) for i in range(n1)]
    while timer <= end:
        # record coalescent events
        size = N1[timer]
        draw = [randint(0, size - 1) for i in range(len(interiors1))]
        table = {}
        for i in range(len(draw)):
            try:
                table[draw[i]].append(i)
            except KeyError:
                table[draw[i]] = [i]
        events = [val for key, val in table.items() if len(val) >= 2]
        # record pairwise ibd
        drop = []
        for merge in events:
            # incoporate incoming node data into next level interior node
            itr += 1
            new_node = WfNode(itr, timer)
            nodes = []
            for i in merge:
                nodes.append(interiors1[i])
                drop.append(i)
            new_node.multiple_coalesce(nodes, ibd_min_span, ibs_min_span)
            # multiple coalescent, often 2
            interiors1.append(new_node)
        # dropping old nodes
        drop = sorted(drop)
        popr = 0
        for d in drop:
            interiors1.pop(d + popr)
            popr -= 1
        timer += 1
    numTracts1 = sum([node._num_tracts for node in interiors1])

    # population 0
    timer = 1
    itr += 1
    interiors0 = [WfNode(i, 0) for i in range(itr, itr + n0)]
    itr += n0
    while timer <= end:
        # record coalescent events
        size = N0[timer]
        draw = [randint(0, size - 1) for i in range(len(interiors0))]
        table = {}
        for i in range(len(draw)):
            try:
                table[draw[i]].append(i)
            except KeyError:
                table[draw[i]] = [i]
        events = [val for key, val in table.items() if len(val) >= 2]
        # record pairwise ibd
        drop = []
        for merge in events:
            # incoporate incoming node data into next level interior node
            itr += 1
            new_node = WfNode(itr, timer)
            nodes = []
            for i in merge:
                nodes.append(interiors0[i])
                drop.append(i)
            new_node.multiple_coalesce(nodes, ibd_min_span, ibs_min_span)
            # multiple coalescent, often 2
            interiors0.append(new_node)
        # dropping old nodes
        drop = sorted(drop)
        popr = 0
        for d in drop:
            interiors0.pop(d + popr)
            popr -= 1
        timer += 1
    numTracts0 = sum([node._num_tracts for node in interiors0])

    # pre adaptation
    itr += 1
    end = max(Ne.keys())
    interiors = interiors0 + interiors1
    while timer <= end:
        # record coalescent events
        size = Ne[timer]
        draw = [randint(0, size - 1) for i in range(len(interiors))]
        table = {}
        for i in range(len(draw)):
            try:
                table[draw[i]].append(i)
            except KeyError:
                table[draw[i]] = [i]
        events = [val for key, val in table.items() if len(val) >= 2]
        # record pairwise ibd
        drop = []
        for merge in events:
            # incoporate incoming node data into next level interior node
            itr += 1
            new_node = WfNode(itr, timer)
            nodes = []
            for i in merge:
                nodes.append(interiors[i])
                drop.append(i)
            new_node.multiple_coalesce(nodes, ibd_min_span, ibs_min_span)
            # multiple coalescent, often 2
            interiors.append(new_node)
        # dropping old nodes
        drop = sorted(drop)
        popr = 0
        for d in drop:
            interiors.pop(d + popr)
            popr -= 1
        timer += 1
    numTracts = sum([node._num_tracts for node in interiors])
    maxIbdGroup = max([len(c) for c in nx.connected_components(H)])

    return numTracts, maxIbdGroup, numTracts1, numTracts0
