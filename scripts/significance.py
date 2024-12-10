import sys
from scipy.stats import norm
import numpy as np
import pandas as pd
ibd,analytical,simulation,out,out2,=sys.argv[1:]

the_dictionary = dict()
with open(analytical,'r') as f:
    for line in f:
        key,val = line.strip().split('\t')
        the_dictionary[key] = float(val)
pvalue = the_dictionary['p-value']
lwrbnd = the_dictionary['initial-lower-bound']
uprbnd = the_dictionary['initial-upper-bound']
rmean = the_dictionary['revised-mean']
rstd = the_dictionary['revised-standard-deviation']
dz = the_dictionary['upper-discrete-z']
dp = norm.sf(dz)
cz = the_dictionary['upper-continuous-z']
cp = norm.sf(cz)
draw = the_dictionary['upper-discrete-raw']
craw = the_dictionary['upper-continuous-raw']

sims = []
with open(simulation,'r') as f:
    for line in f:
        sims.append(float(line.strip()))
simcut = np.quantile(sims,1-pvalue)
simraw = simcut * rstd + rmean
simp = norm.sf(simcut)

table = pd.read_csv(ibd,sep='\t')
table['UPPER_ANALYTICAL'] = draw
table['Z_UPPER_ANALYTICAL'] = dz
table['GW_LEVEL_ANALYTICAL'] = dp
table['UPPER_SIMULATE'] = simraw
table['Z_UPPER_SIMULATE'] = simcut
table['GW_LEVEL_SIMULATE'] = simp
table['ADJ_MEAN'] = rmean
table['ADJ_STDDEV'] = rstd
table['INIT_LOWER_BOUND'] = lwrbnd
table['INIT_UPPER_BOUND'] = uprbnd
table['UPPER_CONTINUOUS'] = craw
table['Z_UPPER_CONTINUOUS'] = cz
table['GW_LEVEL_CONTINUOUS'] = cp
table['PVALUE'] = pvalue
table.to_csv(out,sep='\t',header=True,index=False)

subtable = table[table['COUNT'] > draw]
subtable.to_csv(out2,sep='\t',header=True,index=False,)