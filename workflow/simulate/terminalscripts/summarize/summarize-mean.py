
import os
import sys
import gzip
import logging
import pandas as pd

# inputting
dirin, filein, fileout, left, right = sys.argv[1:]
left=float(left)
right=float(right)

g = open(dirin + '/' + fileout, 'w')

cmmeans=[]
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
        meCM=(datfra['WEIGHT']*datfra['CMWINDOW']).sum()
        cmmeans.append(meCM)
        # meBP=(tab['WEIGHT']*tab['BPWINDOW']).sum()
        # bpmeans.append(meBP)
        f.close()
    except NotADirectoryError:
        pass

# make a dataframe

# print a summary
avg=sum(cmmeans)/len(cmmeans)
print('Average of mean: '+str(avg))
minmode=min(cmmeans)
maxmode=max(cmmeans)
print('Minimum of mean: '+str(minmode))
print('Maximum of mean: '+str(maxmode))

# writing mode file over replicate study
g.write('MEAN\n')
num=len(cmmeans)
for n in range(num):
    m=cmmeans[n]
    g.write(str(m)); g.write('\n')
g.close()
