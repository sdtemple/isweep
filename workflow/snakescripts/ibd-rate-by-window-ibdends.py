#!/bin/python

# importing
import sys
import numpy as np
import pandas as pd

# inputting
# end: basepair in Mb
# bp: basepair windows in kb
# min_span: cM threshold
# telocut: trimming telomere ends

map_file=snakemake.input.mapin
start=0
end=snakemake.config['FIRSTPASS']['GENOMEEND']
telocut=snakemake.config['FIRSTPASS']['TELOCUT']
by=snakemake.config['FIRSTPASS']['BY']
min_span=snakemake.config['FIRSTPASS']['IBDCUTOFF']
input_file=snakemake.input.filein
output_file=snakemake.output.fileout

start=0
end=int(float(end))*1000000
by=int(float(by))*1000
cmcut=telocut
cmcut=float(cmcut)
min_span=float(min_span)

# formatting
table = pd.read_csv(input_file, sep='\t', compression='gzip')
columns=list(table.columns)
cmcol=columns[11]
encol=columns[10]
stcol=columns[9]
table = table[table[cmcol] >= min_span]
table = table[table[encol] >= start]
table = table[table[stcol] <= end]
windows = [i for i in range(start, end+by, by)]

# setting up for genetic distance
f=open(map_file,'r')
currline=f.readline().strip()
currlist=currline.split('\t')
currcm=float(currlist[2])
currbp=float(currlist[3])
prevcm=currcm
prevbp=currbp
currline=f.readline().strip()
currlist=currline.split('\t')
currcm=float(currlist[2])
currbp=float(currlist[3])
cms=[]

# getting genetic distance
for w in windows:
    while w>=currbp:
        if currline=='':
            bpend=currbp
            cmend=currcm
            break
        prevline=currline
        prevcm=currcm
        prevbp=currbp
        prevline=currline
        currline=f.readline().strip()
        if currline!='':
            currlist=currline.split('\t')
            currcm=float(currlist[2])
            currbp=float(currlist[3])
        else:
            currcm=prevcm
            currbp=prevbp
            bpend=currbp
            cmend=currcm
    if currline=='':
        bpend=currbp
        cmend=currcm
    else:
        b=currbp-w
        a=w-prevbp
        d=a+b
        A=a/d
        B=b/d
        cm=B*prevcm+A*currcm
        cms.append(cm)
f.close()

# counting ibd segments
table=table[table[encol]<=bpend]
windows=[w for w in windows if w <= bpend]
counts = [((table[stcol] <= i) & (table[encol] >= i)).sum() for i in windows]
counter = {'BPWINDOW':windows,'CMWINDOW':cms,'COUNT':counts}
counter = pd.DataFrame(counter)
# cm cutting
counter = counter[counter['CMWINDOW']<=(cmend-cmcut)]
counter = counter[counter['CMWINDOW']>=cmcut]
# count cutting
counter = counter[counter['COUNT']>0]
# saving
counter.to_csv(output_file, sep='\t', index=False)
