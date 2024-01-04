# find the rank of the causal allele

# importing
import pandas as pd
import sys

# inputting
filein, fileout, loc, typ = sys.argv[1:]
loc = int(float(loc))

tab = pd.read_csv(filein, sep='\t')
if typ == '1':
    try:
        sub = tab[tab['POS']==loc]
        idx = sub.index.tolist()[0] + 1
        pos = sub['POS'].tolist()[0]
        aaf = sub['AAF'].tolist()[0]
    except KeyError:
        # why does this happen?
        # is polymorphism not present?
        # i think i fixed this w/ .tolist()
        idx='NA'
        pos='NA'
        aaf='NA'
else: # for isafe
    try:
        sub = tab[tab['POS']==loc]
        idx = sub.index.tolist()[0] + 1
        pos = sub['POS'].tolist()[0]
        aaf = sub['DAF'].tolist()[0]
    except KeyError:
        # why does this happen?
        # is polymorphism not present?
        # i think i fixed this w/ .tolist()
        idx='NA'
        pos='NA'
        aaf='NA'

# writing
f=open(fileout,'w')
f.write('RANK\tPOS\tAAF\n')
f.write(str(idx)); f.write('\t')
f.write(str(pos)); f.write('\t')
f.write(str(aaf)); f.write('\n')
