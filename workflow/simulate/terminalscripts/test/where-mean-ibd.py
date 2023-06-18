import sys
import gzip
import pandas as pd
import numpy as np
input_file, = sys.argv[1:]
f = gzip.open(input_file, 'rt')
f.readline()
cmwindows = []
bpwindows = []
counts = []
for line in f:
    bpwindow, cmwindow, count = line.strip().split('\t')
    bpwindows.append(int(float(bpwindow)))
    cmwindows.append(float(cmwindow))
    # use comments for cM based analysis
    counts.append(int(float(count)))
dat={'COUNT':counts,'CMWINDOW':cmwindows,'BPWINDOW':bpwindows}
datfra=pd.DataFrame(dat)
datfra['WEIGHT']=datfra['COUNT']/datfra['COUNT'].sum()
meBP=(datfra['WEIGHT']*datfra['BPWINDOW']).sum()
# caution: this will find the first window
# w/ the mean ibd segment count
print(int(meBP))
f.close()