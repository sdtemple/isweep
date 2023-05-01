# setting up
from isweep import *
short_ibd=snakemake.input.ibdin
short_vcf=snakemake.input.vcfin
K=int(float(snakemake.config["ISWEEPPARAMS"]["IBDCOMMK"]))
F=float(snakemake.config["ISWEEPPARAMS"]["IBDCOMMF"])
Q1=float(snakemake.config["ISWEEPPARAMS"]["IBDCOMMQ1"])
Q2=float(snakemake.config["ISWEEPPARAMS"]["IBDCOMMQ2"])
header=int(float(snakemake.params.header))
scalar=int(float(snakemake.config["ISWEEPPARAMS"]["OUTLIERSCALAR"]))
fileout=snakemake.output.fileout

# forming graph
segs = read_ibd_file(short_ibd, header = header, include_length = 0)
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
table.to_csv(fileout,sep='\t',index=False)
