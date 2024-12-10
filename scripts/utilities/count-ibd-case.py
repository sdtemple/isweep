# Computing ibd counts over windows

# importing
import sys
import numpy as np
import pandas as pd

# inputting
input_file,map_file,output_file,start,end,ind1,ind2,casefile=sys.argv[1:]

# formatting
table = pd.read_csv(input_file, sep='\t')
map_file = pd.read_csv(map_file,sep='\t')
columns=list(table.columns)
encol=columns[int(float(end))]
stcol=columns[int(float(start))]
ind1col=columns[int(float(ind1))]
ind2col=columns[int(float(ind2))]

casedict=dict()
with open(casefile) as f:
    for line in f:
        ind,status=line.strip().split('\t')
        casedict[int(ind)]=int(float(status))   
table['case1']=table[ind1col].map(casedict)
table['case2']=table[ind2col].map(casedict)
table['match']=table['case1']==table['case2']
table['casemult']=table['case1']*table['case2']
table['case']=table.apply(lambda x: x['casemult'] if x['match'] else 2, axis=1)


# counting ibd segments
counts = [((table[stcol] <= i) & (table[encol] >= i)).sum() for i in map_file['bp']]
counts0 = [((table[stcol] <= i) & (table[encol] >= i) & (table['case']==0)).sum() for i in map_file['bp']]
counts1 = [((table[stcol] <= i) & (table[encol] >= i) & (table['case']==1)).sum() for i in map_file['bp']]
counter = {'BPWINDOW':map_file['bp'].to_list(),
           'CMWINDOW':map_file['cm'].to_list(),
           'COUNT':counts,
           'COUNT0':counts0,
           'COUNT1':counts1
          }
counter = pd.DataFrame(counter)

# count cutting
counter = counter[counter['COUNT']>0]

# saving
counter.to_csv(output_file, sep='\t', index=False)