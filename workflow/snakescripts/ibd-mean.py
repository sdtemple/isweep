# get mean of ibd count windows

import pandas as pd
import numpy as np
import sys

filein, = sys.argv[1:]

# formatting
tab=pd.read_csv(filein, sep='\t')
tab['WEIGHT']=tab['COUNT']/tab['COUNT'].sum()

# central tendency: mean
meCM=(tab['WEIGHT']*tab['CMWINDOW']).sum()
meBP=(tab['WEIGHT']*tab['BPWINDOW']).sum()
meBP=int(meBP)

# to stdout
print(meBP)
