# importing
import os
import pandas as pd
import sys

# inputting
dirin, fileout, selcoef = sys.argv[1:]

# reading
f=open(dirin+'/'+fileout,'w')
f.write('SELCOEF\tISWEEPSELCOEF\tCLUESSELCOEF\n')

# computing
for subdir in os.listdir(dirin):

    if subdir!=fileout:

        isweepfile2=dirin+'/'+subdir+'/isweep.inference.tsv'
        cluesfile2=dirin+'/'+subdir+'/clues.v1.selcoef'
        # isweep inference
        try:
            h=open(cluesfile2,'r')
            # will error if clues file does not exist
            g=open(isweepfile2,'r')
            # will error if isweep file does not exist
            g.readline()
            line=g.readline().strip().split('\t')
            g.close()
            f.write(selcoef); f.write('\t')
            f.write(line[2]); f.write('\t')
            line=h.readline().strip()
            h.close()
            f.write(line)
            f.write('\n')
        except:
            pass

# saving
f.close()
