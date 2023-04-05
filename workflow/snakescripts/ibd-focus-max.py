# get + some cM distance + max of mode, median, mean of ibd count windows

import pandas as pd
import numpy as np
import sys

filein, cmpm = sys.argv[1:]
cmpm=float(cmpm)

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


# cM leftstream
leftstream=tab[tab['CMWINDOW']>=(maxCM+cmpm)]['BPWINDOW'].tolist()[0]

# to stdout
print(leftstream) # actually rightstream
