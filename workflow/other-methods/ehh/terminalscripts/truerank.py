# importing
import pandas as pd
import sys

# inputting
filein, fileout, loc, ifheader, importantcolumn = sys.argv[1:]
loc = int(float(loc))
ifheader = int(float(ifheader))
importantcolumn = int(float(importantcolumn))

if ifheader:
     pass
else:
    ifheader = None


tab = pd.read_csv(filein, 
                  sep='\t', 
                  header=ifheader)
thecolumnname = list(tab.columns)[importantcolumn]
sub = tab[tab[thecolumnname]==loc]

# do something
try:
    idx = sub.index.tolist()[0] + 1
except:
    idx = 'NaN'

f=open(fileout,'w')
f.write(str(idx))
f.write('\n')
f.close()
