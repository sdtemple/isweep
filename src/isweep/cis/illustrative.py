from .cis.cisUtilities import *
from .cis.coalescentIBD import *
import matplotlib.pyplot as plt
import numpy as np

def plot_selection_Ne(s, p, Ne, xlimmax=300, colors=['tab:red','tab:blue','tab:pink','tab:cyan'], title='', loc='upper right', log_scale=False, show_neutral=True):
    '''Make an illustrative plot of how selection impacts Ne over time
    
    Parameters:
    -----------
    s : float
        Selection coefficient
    p : float
        Adaptive allele frequency
    Ne : dict
        Effective population sizes
    xlimmax : float
        Maximum generation on x-axis
    colors : list
        4 colors in matplotlib
    title : str
    loc : str
        Location of legend
    log_scale : bool
        To put y-axis on log10 scale
    show_neutral : bool
        To plot Ne without selection for comparison
    
    Returns
    -------
    None
    '''
    assert len(colors) == 4
    assert p <= 1
    assert p >= 0
    # with selection
    ps, Ns, xs = walk_variant_backward(s, p, Ne)
    G = max(Ne.keys())
    N1 = [ceil(i) for i in Ns * ps]
    N1 = {i:N1[i] for i in range(len(N1))}
    qs = extend_vector(1 - ps, 1, G)
    Ns = np.array(list(Ne.values()))
    N0 = [ceil(i) for i in Ns * qs]
    N0 = {i:N0[i] for i in range(len(N0))}
    g=min(xlimmax,G)
    y1=np.array([val for val in N1.values()])
    x1=np.array([key for key in N1.keys()])
    y1=y1[:g]
    x1=x1[:g]
    if log_scale:
        y1=np.log10(y1)
    y0=np.array([val for val in N0.values()])
    x0=np.array([key for key in N0.keys()])
    y0=y0[:g]
    x0=x0[:g]
    if log_scale:
        y0=np.log10(y0)
    plt.plot(x1,y1, color=colors[0], label='Derived (adaptive)')
    plt.plot(x0,y0, color=colors[1], label='Ancestral (adaptive)')
    
    # without selection
    if show_neutral:
        ps, Ns, xs = walk_variant_backward(0, p, Ne)
        G = max(Ne.keys())
        N1 = [ceil(i) for i in Ns * ps]
        N1 = {i:N1[i] for i in range(len(N1))}
        qs = extend_vector(1 - ps, 1, G)
        Ns = np.array(list(Ne.values()))
        N0 = [ceil(i) for i in Ns * qs]
        N0 = {i:N0[i] for i in range(len(N0))}
        g=min(xlimmax,G)
        y1=np.array([val for val in N1.values()])
        x1=np.array([key for key in N1.keys()])
        y1=y1[:g]
        x1=x1[:g]
        if log_scale:
            y1=np.log10(y1)
        y0=np.array([val for val in N0.values()])
        x0=np.array([key for key in N0.keys()])
        y0=y0[:g]
        x0=x0[:g]
        if log_scale:
            y0=np.log10(y0)
        plt.plot(x1,y1, color=colors[2], label='Derived (neutral)')
        plt.plot(x0,y0, color=colors[3], label='Ancestral (neutral)')
    plt.legend(loc=loc)
    plt.title(title)
    plt.xlim(-10,g+10)
    plt.xlabel('Generation')
    if log_scale:
        plt.ylabel('Log10 Effective size')
    else:
        plt.ylabel('Effective size')
    return None

def plot_illustrative_ibd_graph(ks,js, colors=['tab:red','tab:blue'], node_size=5):
    '''Plot an example of an IBD graph
    
    Parameters
    ----------
    ks : list
        Counts, where some counts far exceed counts in js
    js : list
        Counts
    colors : list
        matplotlib colors
    node_size : float
    
    Returns
    -------
    None
    '''
    
    def make_ibd_clique_list(k, itr):
        start = itr
        edges = []
        for j in range(itr,itr+k):
            for i in range(itr+j-itr+1,itr+k):
                edge = (i,j)
                edges.append(edge)
        return edges

    def make_ibd_graph(ks, minlen=2):
        itr=0
        g=nx.Graph()
        for k in ks:
            edges=make_ibd_clique_list(k,itr)
            g.add_edges_from(edges)
            itr+=k
        return g
    
    color_mapk=[colors[0]]*sum(ks)
    color_mapj=[colors[1]]*sum(js)
    color_map=color_mapk+color_mapj
    g=make_ibd_graph(ks)
    h=make_ibd_graph(js)
    f=nx.disjoint_union(g,h)
    nx.draw(f,node_size=node_size, node_color=color_map)
    legend_elements = [
    Line2D([0], [0], marker='o', color='w', label='Adaptive', markerfacecolor='tab:red', markersize=node_size),
    Line2D([0], [0], marker='o', color='w', label='Ancestral', markerfacecolor='tab:blue', markersize=node_size),        
    ]
    plt.legend(handles=legend_elements)
    return None