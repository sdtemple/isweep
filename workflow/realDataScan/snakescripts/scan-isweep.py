#!/bin/python

# Write locations where ibd exceeds a significance threshold
# Seth D. Temple, sdtemple@uw.edu
# April 26, 2023

# importing
import sys
import pandas as pd
import numpy as np

# inputtting
folder=snakemake.config['CHANGE']['FOLDERS']['STUDY']
folder+='/ibdsegs/modified/scan/'
chrlow=snakemake.config['CHANGE']['ISWEEP']['CHRLOW']
chrhigh=snakemake.config['CHANGE']['ISWEEP']['CHRHIGH']
cutoff1=snakemake.config['FIXED']['ISWEEP']['SCANCUTOFF']
cutoff2=snakemake.config['FIXED']['ISWEEP']['TELOCUTOFF']

# type casting
chrlow=int(float(chrlow))
chrhigh=int(float(chrhigh))
cutoff1=int(float(cutoff1))
cutoff2=int(float(cutoff2))

# reading in data
tab=pd.read_csv(folder+'chr'+chrlow+'.ibd.gz',sep='\t')
tab['CHROM']=chrlow
for i in range(chrlow+1,chrhigh+1):
    tabnow=pd.read_csv(folder+'chr'+i+'.ibd.gz',sep='\t')
    tabnow['CHROM']=i
    tab=pd.concat((tab,tabnow))

# calculating excess ibd
medi=np.quantile(tab['COUNT'],0.5)
stdv=tab['COUNT'].std()
a=medi-stdv*cutoff2
b=medi+stdv*cutoff2
sub=tab
sub=sub[sub['COUNT']>=a]
sub=sub[sub['COUNT']<=b]
medi=np.quantile(sub['COUNT'],0.5)
stdv=sub['COUNT'].std()
b=medi+stdv*cutoff1
out=sub[sub['COUNT']>=b]

# saving
out.to_csv(folder+'excess.ibd.tsv',sep'\t',index=False)
