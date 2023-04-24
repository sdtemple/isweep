#!/bin/python
# read and prepare IBDNe bootstraps for isweep
# seth d temple, sdtemple@uw.edu
# april 24, 2023
import sys
filein,folderout=sys.argv[1:]
fileoutpre='ibdneboot'
fileoutsuf='.ne'
fileitr=1
with open(filein,'r') as f:
    for line in f:
        g=open(folderout+'/'+fileoutpre+str(fileitr)+fileoutsuf,'w')
        g.write('GEN\tNE\n')
        dat=line.strip().split('\t')
        genitr=0
        for d in dat:
            g.write(str(genitr));g.write('\t')
            g.write(d);g.write('\n')
            genitr+=1
        g.close()
        fileitr+=1
