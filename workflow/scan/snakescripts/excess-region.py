#!/bin/python

# Combine windows in a region of excess IBD
# Seth D Temple, sdtemple@uw.edu
# April 27, 2023

import sys
filein=snakemake.input.filein
fileout=snakemake.output.fileout
cm=snakemake.fileout.cMgap
cm=float(cm)

g=open(fileout,'w')
g.write('CHROM\tBPLEFT\tBPRIGHT\tCMLEFT\tCMRIGHT\n')
with open(filein,'r') as f:
    f.readline()
    line=f.readline()
    row=line.strip().split('\t')
    prev=row
    towrite=prev
    for line in f:
        row=line.strip().split('\t')
        if int(row[3]) == int(prev[3]):
            if float(row[1]) >= (float(prev[1])+cm):
                # gap in chromosome
                g.write(towrite[3]); g.write('\t')
                g.write(towrite[0]); g.write('\t')
                g.write(prev[0]); g.write('\t')
                g.write(towrite[1]); g.write('\t')
                g.write(prev[1]); g.write('\n')
                towrite=row
        else:
            # chromosome changed
            g.write(towrite[3]); g.write('\t')
            g.write(towrite[0]); g.write('\t')
            g.write(prev[0]); g.write('\t')
            g.write(towrite[1]); g.write('\t')
            g.write(prev[1]); g.write('\n')
            towrite=row
        prev=row
g.write(towrite[3]); g.write('\t')
g.write(towrite[0]); g.write('\t')
g.write(prev[0]); g.write('\t')
g.write(towrite[1]); g.write('\t')
g.write(prev[1]); g.write('\n')
g.close()
