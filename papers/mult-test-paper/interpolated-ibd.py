# Computing ibd counts over windows

# importing
import sys
import numpy as np
import pandas as pd

# inputting
input_file,map_file,output_file,start,end=sys.argv[1:]

# formatting
table = pd.read_csv(input_file, sep='\t')
map_file = pd.read_csv(map_file,sep='\t')
columns=list(table.columns)
encol=columns[int(float(end))]
stcol=columns[int(float(start))]

# counting ibd segments
counts = [((table[stcol] <= i) & (table[encol] >= i)).sum() for i in map_file['bp']]
counter = {'BPWINDOW':map_file['bp'].to_list(),
           'CMWINDOW':map_file['cm'].to_list(),
           'COUNT':counts
          }
counter = pd.DataFrame(counter)

# count cutting
counter = counter[counter['COUNT']>0]

# saving
counter.to_csv(output_file, sep='\t', index=False)