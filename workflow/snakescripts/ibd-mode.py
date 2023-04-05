# get mode of ibd count windows

import pandas as pd
import numpy as np
import sys

filein, = sys.argv[1:]

# formatting
tab=pd.read_csv(filein, sep='\t')

# central tendency: mode
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

# to stdout
print(moBP)
