# get mode, median, mean of ibd count windows

import pandas as pd
import numpy as np

filein=snakemake.input.filein
fileout=snakemake.output.fileout

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

# write output
f=open(fileout,'w')
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
