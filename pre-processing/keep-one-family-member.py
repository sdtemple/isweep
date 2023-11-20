# import
import sys
import gzip
import networkx as nx

# prepare
filein,fileout=sys.argv[1:]

# form graph; find connected components
G = nx.Graph()
with gzip.open(filein, 'rt') as f:
    f.readline()
    for line in f:
        cells = line.strip().split('\t')
        id1 = cells[0]
        id2 = cells[1]
        G.add_edge(id1,id2)
components = [list(comp) for comp in nx.connected_components(G)]

# ### probably slower
# import pandas as pd
# table=pd.read_csv(filein,sep='\t')
# M = table.shape[0]
# G = nx.Graph()
# for m in range(M):
#     row = table.iloc[m]
#     id1 = row[0]
#     id2 = row[1]
#     G.add_edge(id1,id2)
# components = [list(comp) for comp in nx.connected_components(G)]

# write to file
fileout1 = fileout + '.keep.txt'
fileout2 = fileout + '.excl.txt'
g = open(fileout1, 'w')
h = open(fileout2, 'w')
for comp in components:
    idkeep = comp[0]
    g.write(str(idkeep))
    g.write('\n')
    for idexcl in comp[1:]:
        h.write(str(idexcl))
        h.write('\n')
h.close()
g.close()
        
    