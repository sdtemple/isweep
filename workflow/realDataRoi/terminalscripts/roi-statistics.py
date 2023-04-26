#!/bin/bash

# Compute mean, median, mode, min, max of ibd in region of interest
# Seth D. Temple, sdtemple@uw.edu
# April 25, 2023

# setup
import sys
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
folder=sys.argv[1]
roi=sys.argv[2]
left=sys.argv[3]
right=sys.argv[4]
filein=folder+'/excess.ibd.tsv'
fileout1=folder+roi+'excess.ibd.statistics.tsv'
fileout2=folder+roi+'excess.ibd.statistics.png'

# formatting
tab=pd.read_csv(filein, sep='\t')
tab['WEIGHT']=tab['COUNT']/tab['COUNT'].sum()
tab['CUMSUM']=np.cumsum(tab['WEIGHT'])

# central tendency
shape0=tab.shape[0]
mx=0
moCM=0
moBP=0
for j in range(shape0):
    row=tab.iloc[j]
    if row['COUNT'] > mx:
        moCM = row['CMWINDOW']
        moBP = row['BPWINDOW']
        mx = row['COUNT']
moBP=int(moBP)
meCM=(tab['WEIGHT']*tab['CMWINDOW']).sum()
meBP=(tab['WEIGHT']*tab['BPWINDOW']).sum()
meBP=int(meBP)
mdCM=tab[tab['CUMSUM']>=0.5]['CMWINDOW'].tolist()[0]
mdBP=tab[tab['CUMSUM']>=0.5]['BPWINDOW'].tolist()[0]
minCM=min(moCM,meCM,mdCM)
maxCM=max(moCM,meCM,mdCM)
minBP=min(moBP,meBP,mdBP)
maxBP=max(moBP,meBP,mdBP)

plt.plot(tab['CMWINDOW'],tab['COUNT'],lwd=3,color='k')
plt.xlabel('centiMorgan distance')
plt.ylabel('IBD segment count')
plt.title(roi)
plt.axvline(moCM,lty='dotted',lwd=2,color='tab:gray')
plt.axvline(mdCM,lty='dashed',lwd=2,color='tab:gray')
plt.axvline(meCM,lty='solid',lwd=2,color='tab:gray')
plt.savefig(fileout2,dpi=300)

# write output
f=open(fileout1,'w')
f.write('MODECM\tMEDIANCM\tMEANCM\tMINCM\tMAXCM\tMODEBP\tMEDIANBP\tMEANBP\tMINBP\tMAXBP\n')
f.write(str(moCM)); f.write('\t')
f.write(str(mdCM)); f.write('\t')
f.write(str(meCM)); f.write('\t')
f.write(str(minCM)); f.write('\t')
f.write(str(maxCM)); f.write('\t')
f.write(str(moBP)); f.write('\t')
f.write(str(mdBP)); f.write('\t')
f.write(str(meBP)); f.write('\t')
f.write(str(minBP)); f.write('\t')
f.write(str(maxBP)); f.write('\n')
f.close()
