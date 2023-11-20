import sys
import pandas as pd
filein, ifheader, importantcol = sys.argv[1:]
ifheader = int(float(ifheader))
if ifheader:
        pass
else:
        ifheader = None
print(ifheader)
tablein = pd.read_csv(filein, header=ifheader, sep='\t')
#tablein = pd.read_csv(filein, sep='\t', header=0)
print(tablein.head())
tablecols = list(tablein.columns)
sortcol = tablecols[int(float(importantcol))]
tablein.sort_values(by=sortcol,ascending=False,inplace=True,ignore_index=True)
print(tablein.head())
tablein.to_csv(filein+'.sorted',sep='\t',index=False)