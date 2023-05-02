# importing
import sys
from isweep import *

#  i/o
short_ibd, short_vcf, fileout, K, F, Q1, Q2, scalar = sys.argv[1:]
K=int(float(K)/2)
F=float(F)
Q1=float(Q1)
Q2=float(Q2)
scalar=int(float(scalar))

# forming graph
segs = read_ibd_file(short_ibd, header = 0, include_length = 0)
graph = make_ibd_graph(segs)

# detecting communities
communities = diameter_communities(graph, K=K, max_communities=np.inf)
outliers = outlier_communities(communities, scalar=scalar)

# computing adaptive allele frequencies
tup=labeled_allele_frequencies(short_vcf, outliers)

# making table
pos=tup[0]
freq1, freq0, freqm = putative_allele_frequencies(tup[1], tup[2], tup[3])
table = format_allele_table(pos, freq1, freq0, freqm)
table.sort_values(['DELTA'],inplace=True,ascending=False)
table=table[table['AAF1']>=F]
table=table[table['AAF']>=Q1]
table=table[table['AAF']<=Q2]
table.reset_index(inplace=True,drop=True)
tab['DELTANORM']=tab['DELTA']/tab['DELTA'].sum()
table.to_csv(fileout,sep='\t',index=False)
