# imports
import numpy as np
import networkx as nx
from copy import deepcopy

def make_ibd_graph(ibd_segments):
    '''Create IBD graph from pairwise segments

    Parameters
    ----------
    ibd_segments : tuple
        Collection of pairwise segments

    Returns
    -------
    networkx.Graph
        Graph w/ edges if haplotypes have a detectable IBD segment
    '''
    ibdgraph=nx.Graph()
    edgelist=[]
    for edge in ibd_segments:
        edgelist.append((edge[0],edge[1]))
    ibdgraph.add_edges_from(edgelist)
    return ibdgraph

def diameter_communities(graph, K=3, max_communities=np.inf):
    '''Method to find connected communities with max diameter 2*K

    Parameters
    ----------
    graph : networkx.Graph
    K : int
        Default is max diameter 2*3
    max_communities : int
        Default is to find all communities

    Returns
    -------
    list
        Node sets for fully connected subgraphs with max diameter 2*K
    '''

    # initialize community detection loop
    copygraph = deepcopy(graph)
    communities = []
    stayput=True
    if max_communities >= np.inf:
        max_communities = copygraph.number_of_nodes()

    # initialize degree centrality
    M = copygraph.number_of_nodes() - 1
    degrees = nx.degree_centrality(copygraph)
    degrees = {key:round(val*M,0) for key, val in degrees.items()}
    sorted_degrees = sorted([(key, val) for key, val in degrees.items()], key=lambda x:x[1], reverse=True)

    for i in range(max_communities):

        # find community center
        community = dict()
        while stayput:
            try:
                center = sorted_degrees[0][0]
                if copygraph.has_node(center):
                    stayput=False
                else:
                    sorted_degrees.pop(0)
            except:
                return communities
        current_neighbors = [center]
        community[center] = 0

        # find up to K-order neighbors
        k = 1
        while k <= K: # up to K-order neighbors
            new_neighbors = []
            for neighbor in current_neighbors:
                adjacent_neighbors = copygraph.neighbors(neighbor)
                for adjacent_neighbor in adjacent_neighbors:
                    if adjacent_neighbor not in community:
                        community[adjacent_neighbor] = 0
                        new_neighbors.append(adjacent_neighbor)
                copygraph.remove_node(neighbor)
            current_neighbors = new_neighbors
            k += 1
        for neighbor in current_neighbors:
            if copygraph.has_node(neighbor):
                copygraph.remove_node(neighbor)
        communities.append(set(community.keys()))
        stayput=True

    return communities

def community_list(graph, communities):
    '''Return list of subgraphs

    Parameters
    ----------
    graph : networkx.Graph
    communities : list of sets
    '''
    return [graph.subgraph(comm) for comm in communities]

def largest_community(communities):
    ''' Return haploids in largest community

    Parameters
    ----------
    list
        Haplotype/sample IDs in the largest community
    '''
    return list(communities[0])

def outlier_communities(communities, scalar=5):
    '''Combine outlier (or single largest) excess IBD cluster

    Parameters
    ----------
    communities : list
    scalar : float
        Multiplier of numpy.std()

    Returns
    -------
    list
        IDs for haploids in outlier communities
    '''
    communities2 = deepcopy(communities)
    communities2 = sorted(communities2, key=len, reverse=True)
    community_sizes = [len(community) for community in communities]
    community_sizes = np.array(community_sizes)
    cutoff = community_sizes.mean() + community_sizes.std() * scalar
    oindices = (community_sizes >= cutoff).nonzero()
    outliers = set()
    for i in oindices[0]:
        community = communities[i]
        outliers = community.union(outliers)
    outliers = list(outliers)
    if len(outliers) == 0:
        return communities2[0]
    else:
        return outliers
