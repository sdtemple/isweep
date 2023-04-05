# get median of ibd count windows

import pandas as pd
import numpy as np
import sys

filein, = sys.argv[1:]

# formatting
tab=pd.read_csv(filein, sep='\t')
tab['WEIGHT']=tab['COUNT']/tab['COUNT'].sum()
tab['CUMSUM']=np.cumsum(tab['WEIGHT'])

# central tendency: median
mdCM=tab[tab['CUMSUM']>=0.5]['CMWINDOW'].tolist()[0]
mdBP=tab[tab['CUMSUM']>=0.5]['BPWINDOW'].tolist()[0]

# to stdout
print(mdBP)
