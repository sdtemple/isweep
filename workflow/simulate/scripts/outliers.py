# find outlier haplotype groups

# importing
import sys
from isweep import *

#  i/o
short_ibd, folderout, K, scalar = sys.argv[1:]
K=int(float(K)/2)
scalar=int(float(scalar))

# forming graph
segs = read_ibd_file(short_ibd, header = 0, include_length = 0)
graph = make_ibd_graph(segs)

# detecting communities
communities = diameter_communities(graph, K=K, max_communities=np.inf)

# write clusters to file
def write_outliers(communities, folderout, scalar=5):
    '''Write haplotypes in outlier big communities

    Parameters
    ----------
    communities : list of sets
    scalar : float
        Multiplier of numpy .std()

    Returns
    -------
    None
    '''
    community_sizes = [len(community) for community in communities]
    community_sizes = np.array(community_sizes)
    cutoff = community_sizes.mean() + community_sizes.std() * scalar
    idx=1
    for community in communities:
        if len(community) > cutoff:
            f = open(folderout + '/outlier' + str(idx) + '.txt','w')
            for haplo in community:
                f.write(haplo)
                f.write('\n')
            f.close()
            idx += 1
    return None

# saving
write_outliers(communities, folderout, scalar)
