import sys
import numpy as np
from isweep import *

selcoefs = [0.01,0.02,0.03,0.04]

filepre = sys.argv[1]
fileout = sys.argv[2]
p = float(sys.argv[3])
ne = sys.argv[4]
neloc = sys.argv[5]
Ne = read_Ne(neloc)
size = int(sys.argv[6])
cM = float(sys.argv[7])
SIM = int(sys.argv[8])
REP = int(sys.argv[9])
pre = filepre

m = size * 2
M = m * (m-1) / 2

h = open(fileout,'w')

for selcoef in selcoefs:
    print(selcoef)
    filename = '-selcoef' + str(selcoef) + '-freq' + sys.argv[3] + '-' + ne + '-size' + str(size) + '-cM' + str(cM) + '-sim' + str(SIM) + '-rep' + str(REP)

    g=open(pre+filename+'.txt')
    g.readline()
    g.readline()
    g.readline()
    columns = g.readline().strip().split(',')
    table = dict()
    for col in columns:
        table[col] = []
    J = len(columns)
        
    g.readline()
    while True:
        row = g.readline().strip()
        reps = row.split('\t')
        if row == '':
            break
        for rep in reps:
            details = rep.split(',')
            for j in range(J):
                try:
                    table[columns[j]].append(int(details[j]))
                except:
                    table[columns[j]].append(np.NaN)
    
    g.close()
    
    pdtable = pd.DataFrame(table)

    h.write('selcoef: ' + str(selcoef) + '\n')

    # for i in range(2000):
    for i in range(pdtable.shape[0]-1):
        if (i % 1000) == 0:
            print(i)
        numTracts = list(pdtable['num_tracts'])[i]
        sestim = minimize_scalar(chi2_isweep, args=(p,Ne,M,(numTracts,),(cM,np.inf),'a'),
                     bounds=(0,0.5),
                     method='bounded'
                    ).x
        h.write(str(sestim))
        h.write('\t')
    numTracts = list(pdtable['num_tracts'])[-1]
    # numTracts = pdtable['num_tracts'][2000-1]
    sestim = minimize_scalar(chi2_isweep, args=(p,Ne,M,(numTracts,),(cM,np.inf),'a'),
                         bounds=(0,0.5),
                         method='bounded'
                        ).x
    h.write(str(sestim))
    h.write('\n')

h.close()
