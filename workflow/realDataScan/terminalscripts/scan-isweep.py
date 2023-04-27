#!/bin/python

# Write locations where ibd exceeds a significance threshold
# Seth D. Temple, sdtemple@uw.edu
# April 26, 2023

# importing
import sys
import pandas as pd
import numpy as np

# inputtting
folder,chrlow,chrhigh,cutoff1,cutoff2=sys.argv[1:]
firstfolder=folder
folder+='/ibdsegs/modified/scan/'

# type casting
chrlow=int(float(chrlow))
chrhigh=int(float(chrhigh))
cutoff1=int(float(cutoff1))
cutoff2=int(float(cutoff2))

# reading in data
tab=pd.read_csv(folder+'chr'+str(chrlow)+'.ibd.windowed.tsv.gz',sep='\t')
tab['CHROM']=chrlow
for i in range(chrlow+1,chrhigh+1):
    tabnow=pd.read_csv(folder+'chr'+str(i)+'.ibd.windowed.tsv.gz',sep='\t')
    tabnow['CHROM']=i
    tab=pd.concat((tab,tabnow))
# saving all
tab.to_csv(firstfolder+'/scan.ibd.tsv',sep='\t',index=False)

# calculating excess ibd
medi=np.quantile(tab['COUNT'],0.5)
stdv=tab['COUNT'].std()
a=medi-stdv*cutoff2
b=medi+stdv*cutoff2
sub=tab
sub=sub[(sub['COUNT']>=a)&(sub['COUNT']<=b)]
medi=np.quantile(sub['COUNT'],0.5)
stdv=sub['COUNT'].std()
b=medi+stdv*cutoff1
out=tab[tab['COUNT']>=b]

# saving excess
out.to_csv(firstfolder+'/excess.ibd.tsv',sep'\t',index=False)
