import gzip
import numpy as np
import pandas as pd
import networkx as nx
import matplotlib.pyplot as plt
from random import randint
from math import floor, exp, log, ceil
from numpy.random import exponential as Exp
from scipy.stats import norm, binom
from scipy.optimize import minimize, minimize_scalar
from .cis.coalescentIBD import walk_variant_backward, probability_quasi_geometric, probability_erlang_segments
from .isweepUtilities import bin_ibd_segments

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

def chi2_labeled_isweep(s, p0, Ne, n, obs1, obs0, ab, one_step_model = 'm', tau0 = 0, ploidy = 2):
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
    obs1 : array-like
        Observed counts for ibd segment bins (labeled 1)
    obs0 : array-like
        Observed counts for ibd segment bins (labeled 0)
    ab : array-like
        Increasing floats in centiMorgans
    one_step_model : str
        'm', 'a', 'd', or 'r'
    tau0 : int
        Generation at which neutrality begins

    Returns
    -------
    float
        Goodness-of-fit statistic
    '''

    # assertions
    assert len(obs1) == len(obs0)
    assert len(obs1) == (len(ab) - 1)

    # probabilities
    ps, Ns, xs = walk_variant_backward(s, p0, Ne, False, one_step_model, tau0, ploidy)
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

def chi2_isweep(s, p0, Ne, n, obs, ab, one_step_model = 'm', tau0 = 0, ploidy = 2):
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
    obs : array-like
        Observed counts for ibd segment bins (unlabeled)
    ab : array-like
        Increasing floats in centiMorgans
    one_step_model : str
        'm', 'a', 'd', or 'r'
    tau0 : int
        Generation at which neutrality begins

    Returns
    -------
    float
        Goodness-of-fit statistic
    '''

    assert len(obs) == (len(ab) - 1)
    ps, Ns, xs = walk_variant_backward(s, p0, Ne, False, one_step_model, tau0, ploidy)
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
    boot : array-like
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
    boot : array-like
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

def bootstrap_hall(val, boot, alpha1 = 0.025, alpha2 = 0.975):
    '''Implement Hall's interval estimator

    Parameters
    ----------
    val : float
        Parameter estimate
    boot : array-like
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

def bootstrap_efron(val, boot, alpha1 = 0.025, alpha2 = 0.975):
    '''Implements Efron's interval estimator

    Parameters
    ----------
    val : float
        Parameter estimate
    boot : array-like
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

def bootstrap_efron_bc(val, boot, alpha1 = 0.025, alpha2 = 0.975):
    '''Implements Efron's interval estimator (w/ bias-correction)

    Parameters
    ----------
    val : float
        Parameter estimate
    boot : array-like
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

    prop = sum(boot > val) / len(boot)
    z0 = norm.ppf(1 - prop)
    za = norm.ppf(alpha1)
    zb = norm.ppf(alpha2)
    ql = norm.cdf(2 * z0 + za)
    qu = norm.cdf(2 * z0 + zb)
    qm = norm.cdf(2 * z0)
    vl = np.quantile(boot, ql)
    vu = np.quantile(boot, qu)
    vm = np.quantile(boot, qm)

    return (vl, vm, vu)

##### time to event #####

def when_count(ct, s, p0, Ne, random_walk = True, one_step_model = 'm', tau0 = 0):
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

    Returns
    -------
    int
        Generation time
    '''

    var = walk_variant_backward(s, p0, Ne, random_walk, one_step_model, tau0)[2]
    idx = np.argmax(var <= ct)
    if idx == 0: # handle variant introduction
        return len(var)

    return idx

def bootstrap_count(ct, B, boots, bootp, Ne, random_walk = True, one_step_model = 'm', tau0 = 0):
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

    Returns
    -------
    array-like
        NumPy array of bootstraps for variant count generation time
    '''

    assert len(boots) == len(bootp)
    boott = []
    for i in range(len(boots)):
        for b in range(B):
            boott.append(when_count(ct, boots[i], bootp[i], Ne, random_walk, one_step_model, tau0))

    return boott

def when_freq(maf, s, p0, Ne, random_walk = True, one_step_model = 'm', tau0 = 0):
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

    Returns
    -------
    int
        Generation time
    '''

    var = walk_variant_backward(s, p0, Ne, random_walk, one_step_model, tau0)[0]
    idx = np.argmax(var <= maf)
    if idx == 0: # handle variant introduction
        return len(var)

    return idx

def bootstrap_freq(maf, B, boots, bootp, Ne, random_walk = True, one_step_model = 'm', tau0 = 0):
    '''Parametric bootstrap for variant frequency generation time

    Parameters
    ----------
    maf : float
        Variant frequency
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

    Returns
    -------
    array-like
        NumPy array of bootstraps for variant frequency generation time
    '''

    assert len(boots) == len(bootp)
    boott = []
    for i in range(len(boots)):
        for b in range(B):
            boott.append(when_freq(maf, boots[i], bootp[i], Ne, random_walk, one_step_model, tau0))

    return boott

### metrics (evaluate experimental results) ###

def check_envelope(x, a, b):
    '''Check if value is in an interval

    Parameters
    ----------
    x : float
        Parameter estimate
    a : float
        Left endpoint
    b : float
        Right endpoint

    Returns
    -------
    bool
        If covered
    '''

    if x >= a:
        if x <= b:
            return True

    return False

def make_confusion_matrix(real, pred):
    '''Form confusion matrix for statistical classification

    Parameters
    ----------
    real : array_like
        Actual 1s and 0s
    pred : array_like
        Predicted 1s and 0s

    Returns
    -------
    pandas DataFrame object
        2 x 2 crosstab
    '''

    data = {'Actual':real, 'Predicted':pred}
    conf = pd.DataFrame(data, columns = ['Actual', 'Predicted'])
    matr = pd.crosstab(conf['Actual'],
                       conf['Predicted'],
                       rownames=['Actual'],
                       colnames=['Predicted'])

    return matr

def false_omission_rate(matrix):
    '''Compute false omission rate

    Parameters
    ----------
    matrix : pandas DataFrame object
        2 x 2 cross tab

    Returns
    -------
    float
        The false omission rate (https://en.wikipedia.org/wiki/Confusion_matrix)
    '''

    tn = matrix[0][0]
    fp = matrix[1][0]
    fn = matrix[0][1]
    tp = matrix[1][1]

    return fn / (fn + tn)

def false_positive_rate(matrix):
    '''Compute false positive rate

    Parameters
    ----------
    matrix : pandas DataFrame object
        2 x 2 cross tab

    Returns
    -------
    float
        The false positive rate (https://en.wikipedia.org/wiki/Confusion_matrix)
    '''

    tn = matrix[0][0]
    fp = matrix[1][0]
    fn = matrix[0][1]
    tp = matrix[1][1]

    return fp / (fp + tn)

def true_positive_rate(matrix):
    '''Compute true positive rate

    Parameters
    ----------
    matrix : pandas DataFrame object
        2 x 2 cross tab

    Returns
    -------
    float
        The true positive rate (https://en.wikipedia.org/wiki/Confusion_matrix)
    '''

    tn = matrix[0][0]
    fp = matrix[1][0]
    fn = matrix[0][1]
    tp = matrix[1][1]

    return tp / (tp + fn)

# old

# def ibs_community(graph, n, min_degree_centrality=0.25, first_frequency_cutoff=0.01, second_frequency_cutoff=0.05, distance_cutoff=4, quantile=0.05, quantile_scalar=1.0, ploidy=2):
#     '''Compute the count of a community in a graph that:
#     (1) has high minimum degree centrality (out of all degree centralities);
#     (2) most node pairings are connected by no more than 2, 4, or 6 edges;
#     (3) includes the node with the highest degree centrality

#     Algorithm (two cases)
#         Second case (highest degree centrality exceeds min_degree_centrality)
#             1. Find node with highest degree centrality, the focal node
#             2. Ignore all nodes with degree centrality less than second_frequency_cutoff
#             3. Add all neighbors of focal node to community, the first order neighbors
#             4. If distance_cutoff >= 4, add all neighbors of neighbors of focal node, the 2nd order neighbors
#             5. If distance_cutoff >= 6, add all 3rd order neighbors
#             6. Add remaining nodes if more edges w/ community than quantile_scalar * quantile(degree centrality of community)
#         First case
#             1. Find node with highest degree centrality, the focal node
#             2. Ignore all nodes with degree centrality less than first_frequency_cutoff
#             3. Add all neighbors of focal node to community
#             4. Add remaining nodes if more edges w/ community than quantile_scalar * quantile(degree centrality of community)

#     Parameters
#     ----------
#     graph : networkx Graph object
#     n : int
#         Sample size (individual)
#     min_normalized_degree_centrality : float
#         Threshold to conduct second case
#     first_frequency_cutoff : float
#         Threshold to insist component has nodes
#     second_frequency_cutoff : float
#         (Second case) Threshold to insist component has nodes
#     distance_cutoff : int
#         (Second case) Threshold for path distances connecting node pairings
#     quantile : float
#     quantile_scalar : float
#         quantile * quantile_scalar is threshold for adding tightly connected nodes back
#     ploidy : int
#         1 for haploid, 2 for diploid, 4 for tetraploid, etc.

#     Returns
#     -------
#     int
#         A community size satisfying (1-3)
#     '''
#     assert first_frequency_cutoff >= 0
#     assert second_frequency_cutoff > first_frequency_cutoff
#     assert second_frequency_cutoff < 1
#     assert min_degree_centrality >= 0
#     assert min_degree_centrality <= 1
#     assert quantile >= 0
#     assert quantile <= 1
#     assert quantile_scalar > 0
#     if distance_cutoff >= 4:
#         if distance_cutoff >= 6:
#             print('Warning: distance_cutoff >= 6 is O(N^3)')
#         else:
#             print('Warning: 6 > distance_cutoff >= 4 is O(N^2)')
#     N = ploidy * n
#     first_count_cutoff = first_frequency_cutoff * N
#     second_count_cutoff = second_frequency_cutoff * N
#     m = graph.number_of_nodes() - 1
#     degrees = nx.degree_centrality(graph)
#     degrees = {key:round(val*m,0) for key, val in degrees.items()}
#     # path length filter
#     degreesT = sorted([(key,val) for key,val in degrees.items()], key=lambda item:item[1], reverse=True)
#     friendliest_node = degreesT[0][0]
#     friendly_count = len(graph[friendliest_node])
#     within = set()
#     removenodes = [key for key, val in degrees.items() if val < first_count_cutoff] # frequency filter
#     graph.remove_nodes_from(removenodes)
#     for node1 in graph.neighbors(friendliest_node):
#         within.add(node1)
#     outnodes = [node for node in graph.nodes() if node not in within]
#     if friendly_count >= (min_degree_centrality * m): # degree centrality filter
#         removenodes = [node for node in outnodes if degrees[node] < second_count_cutoff] # frequency filter
#         outnodes = [node for node in outnodes if degrees[node] >= second_count_cutoff]
#         graph.remove_nodes_from(removenodes)
#         for node1 in graph.neighbors(friendliest_node):
#             if distance_cutoff >= 4:
#                 for node2 in graph.neighbors(node1):
#                     within.add(node2)
#                     if distance_cutoff >= 6:
#                         for node3 in graph.neighbors(node2):
#                             within.add(node3)
#     subgraph = graph.subgraph(within)
#     # add back tightly connected out nodes
#     outnodes = [node for node in graph.nodes() if node not in within]
#     subdegrees = nx.degree_centrality(subgraph)
#     M = subgraph.number_of_nodes() - 1
#     subdegrees = {key:round(val*M,0) for key, val in subdegrees.items()}
#     vals = np.array([val for key, val in subdegrees.items()])
#     quantile_sdeg = np.quantile(vals, quantile)
#     addnodes = []
#     add_thre = quantile_sdeg * quantile_scalar
#     for node1 in outnodes:
#         ctr = 0
#         studynodes = [node for node in graph[node1]]
#         for node in studynodes:
#             if subgraph.has_node(node):
#                 ctr += 1
#         if ctr >= add_thre:
#             addnodes.append(node1)
#     addwithin = within.union(set(addnodes))
#     subgraph = graph.subgraph(addwithin)
#     return subgraph

# def community_estimate_p0(ibs_file, n, min_degree_centrality=0.25, first_frequency_cutoff=0.01, second_frequency_cutoff=0.05, distance_cutoff=4, quantile=0.05, quantile_scalar=1.0, ploidy=2):
#     '''Create a graph based on IBS/HBS segments

#     Parameters
#     ----------
#     ibs_file : str
#         Name of text file with edges (IBS/HBS segments)
#     n : int
#         Sample size (individual)

#     * See ibs_community()
#     min_normalized_degree_centrality : float
#         Threshold to conduct second case
#     first_frequency_cutoff : float
#         Threshold to insist component has nodes
#     second_frequency_cutoff : float
#         (Second case) Threshold to insist component has nodes
#     distance_cutoff : int
#         (Second case) Threshold for path distances connecting node pairings
#     quantile : float
#     quantile_scalar : float
#         quantile * quantile_scalar is threshold for adding tightly connected nodes back

#     ploidy : int
#         1 for haploid, 2 for diploid, 4 for tetraploid, etc.

#     Returns
#     -------
#     float
#         Allele frequency estimate
#     '''
#     # create graph
#     graph = nx.Graph()
#     # read in edges (IBS)
#     edges_list = read_ibd_file(ibs_file, 0, 0)
#     graph.add_edges_from(edges_list)
#     # determine community of interest
#     IbsCommunitySize = len(ibs_community(graph, n, min_degree_centrality, first_frequency_cutoff, second_frequency_cutoff, distance_cutoff, quantile, quantile_scalar, ploidy).nodes())
#     return IbsCommunitySize / n / 2

# def community_segment_counts(ibd_file, ibs_file, n, min_degree_centrality=0.25, first_frequency_cutoff=0.01, second_frequency_cutoff=0.05, distance_cutoff=4, quantile=0.05, quantile_scalar=1.0, ploidy=2):
#     '''Create a graph based on IBS segments

#     Parameters
#     ----------
#     ibd_file : str
#         Name of text file with edges (IBD segments)
#     ibs_file : str
#         Name of text file with edges (IBS/HBS segments)
#     n : int
#         Sample size (individual)

#     * See ibs_community()
#     min_normalized_degree_centrality : float
#         Threshold to conduct second case
#     first_frequency_cutoff : float
#         Threshold to insist component has nodes
#     second_frequency_cutoff : float
#         (Second case) Threshold to insist component has nodes
#     distance_cutoff : int
#         (Second case) Threshold for path distances connecting node pairings
#     quantile : float
#     quantile_scalar : float
#         quantile * quantile_scalar is threshold for adding tightly connected nodes back

#     ploidy : int
#         1 for haploid, 2 for diploid, 4 for tetraploid, etc.


#     Returns
#     -------
#     tuple
#         (Total segment count, size of community 1, segment count of community 1, segment count of community 0)
#     '''
#     # create graph
#     graph = nx.Graph()
#     # read in edges (IBS/HBS)
#     edges_list = read_ibd_file(ibs_file, 0, 0)
#     graph.add_edges_from(edges_list)
#     # determine community of interest
#     IbsCommunity = ibs_community(graph, n, min_degree_centrality, first_frequency_cutoff, second_frequency_cutoff, distance_cutoff, quantile, quantile_scalar, ploidy)
#     IbsCommunityNodes = [node for node in IbsCommunity.nodes()]
#     # label segments 1 or 0
#     segments = read_ibd_file(ibd_file, 1, 0)
#     numTracts1 = 0
#     numTracts0 = 0
#     for segment in segments:
#         if segment[0] in IbsCommunityNodes:
#             numTracts1 += 1
#         elif segment[1] in IbsCommunityNodes:
#             numTracts1 += 1
#         else:
#             numTracts0 += 1
#     return numTracts1 + numTracts0, len(IbsCommunityNodes), numTracts1, numTracts0
