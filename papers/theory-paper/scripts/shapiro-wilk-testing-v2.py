import numpy as np
from scipy.stats import shapiro
import sys

filein, fileout, cutidx, stopidx, validx = sys.argv[1:]

alphas = [0.10,0.05,0.01]
cutidx = int(cutidx)
stopidx = int(stopidx)
validx = int(validx)

g = open(fileout, 'w')
g.write('mean\tmedian\tstddev\tleftend\trightend\tshapirostat\tpvalue')
g.write('\talpha=' + str(alphas[0]))
g.write('\talpha=' + str(alphas[1]))
g.write('\talpha=' + str(alphas[2]))
g.write('\n')

ctr = 0
with open(filein, 'r') as f:
    f.readline()
    f.readline()
    f.readline()
    for line in f:
        vals = line.strip().split('\t')
        vals = [float(val.split(',')[validx]) for val in vals]
        vals = vals[:cutidx]
        vals = np.array(vals)
        av = vals.mean()
        md = np.quantile(vals,0.5)
        left = np.quantile(vals,0.025)
        right = np.quantile(vals,0.975) 
        dv = vals.std()
        sh = shapiro(vals)
        pv = sh.pvalue
        st = sh.statistic
        a1 = 1 if pv <= alphas[0] else 0
        a5 = 1 if pv <= alphas[1] else 0
        a0 = 1 if pv <= alphas[2] else 0
        g.write(str(av)); g.write('\t')
        g.write(str(md)); g.write('\t')
        g.write(str(dv)); g.write('\t')
        g.write(str(left)); g.write('\t')
        g.write(str(right)); g.write('\t')
        g.write(str(st)); g.write('\t')
        g.write(str(pv)); g.write('\t')
        g.write(str(a1)); g.write('\t')
        g.write(str(a5)); g.write('\t')
        g.write(str(a0)); g.write('\n')
        ctr += 1
        if ctr >= stopidx:
            break

g.close()