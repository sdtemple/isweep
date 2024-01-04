# Write locations where ibd exceeds a significance threshold

# importing
import sys
import pandas as pd
import numpy as np

# inputtting
folder,out1,out2,chrlow,chrhigh,cutoff1,cutoff2=sys.argv[1:]

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
tab.to_csv(out1,sep='\t',index=False)

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
out.to_csv(out2,sep='\t',index=False)
