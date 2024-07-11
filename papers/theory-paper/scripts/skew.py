import sys
import numpy as np
from scipy.stats import skew
filein,=sys.argv[1:]
rows = []
with open(filein,'r') as f:
    for line in f:
        row = [int(x) for x in line.strip().split('\t')]
        rows += row
rows = np.array(rows)
print('median')
print(np.quantile(rows,0.5))
print('mean')
print(rows.mean())
print('skew (bias)')
print(skew(rows,bias=True))
print('skew (adjusted)')
print(skew(rows,bias=False))
