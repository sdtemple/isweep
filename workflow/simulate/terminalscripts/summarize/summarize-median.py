
import os
import sys
import gzip
import logging
import pandas as pd
import numpy as np

# inputting
dirin, filein, fileout, left, right = sys.argv[1:]
left=float(left)
right=float(right)

g = open(dirin + '/' + fileout, 'w')

cmmedians=[]
# bpmeans=[]
for subdir in os.listdir(dirin):
    input_file = dirin + '/' + subdir + '/' + filein
    try:
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
        datfra=datfra[(datfra['BPWINDOW']>=left)&(datfra['BPWINDOW']<=right)]
        datfra['WEIGHT']=datfra['COUNT']/datfra['COUNT'].sum()
        datfra['CUMSUM']=np.cumsum(datfra['WEIGHT'])
        mdCM=datfra[datfra['CUMSUM']>=0.5]['CMWINDOW'].tolist()[0]
        # mdBP=datfra[datfra['CUMSUM']>=0.5]['BPWINDOW'].tolist()[0]
        cmmedians.append(mdCM)
        # bpmedians.append(mdBP)
        f.close()
    except NotADirectoryError:
        pass

# make a dataframe

# print a summary
avg=sum(cmmedians)/len(cmmedians)
print('Average of median: '+str(avg))
minmode=min(cmmedians)
maxmode=max(cmmedians)
print('Minimum of median: '+str(minmode))
print('Maximum of median: '+str(maxmode))

# writing mode file over replicate study
g.write('MEAN\n')
num=len(cmmedians)
for n in range(num):
    m=cmmedians[n]
    g.write(str(m)); g.write('\n')
g.close()
