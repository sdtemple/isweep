import sys
import numpy as np
filein,fileout=sys.argv[1:]
percentiles=[0.005,0.025,0.05,0.1,0.9,0.95,0.975,0.995,0.9986501,0.9999683]
rows = []
with open(filein,'r') as f:
    for line in f:
        row = [int(x) for x in line.strip().split('\t')]
        rows += row
rows = np.array(rows)
with open(fileout,'w') as h:
    h.write('percentile\tvalue\n')
    for p in percentiles:
        h.write(str(p)); h.write('\t');
        h.write(str(np.quantile(rows,p))); h.write('\n')
