# importing
import os
import pandas as pd
import sys

# inputting
dirin, fileout, selcoef, loc, top = sys.argv[1:]
loc=int(float(loc))
top=int(float(top))

# reading
f=open(dirin+'/'+fileout,'w')
f.write('POS\tAAF\tISWEEPRANK\tISWEEPAAF\tISWEEPLOC\tISAFERANK\tISAFEAAF\tISAFELOC\tSELCOEF\tISWEEPSELCOEF\n')

# computing
for subdir in os.listdir(dirin):

    if subdir!=fileout:

        isafefile=dirin+'/'+subdir+'/isafe.ranks.tsv.gz'
        isweepfile=dirin+'/'+subdir+'/isweep.ranks.tsv.gz'
        isweepfile2=dirin+'/'+subdir+'/isweep.inference.tsv'
        # isweep
        try:
            tab=pd.read_csv(isweepfile,sep='\t')
            try:
                sub = tab[tab['POS']==loc]
                idx = sub.index.tolist()[0] + 1
                pos = sub['POS'].tolist()[0]
                aaf = sub['AAF'].tolist()[0]
            except:
                pos=loc
                idx=0 # signals not in file
                aaf=0 # signals not in file
            avgaaf=tab.iloc[:top]['AAF'].mean()
            avgloc=tab.iloc[:top]['POS'].mean()
            f.write(str(pos)); f.write('\t')
            f.write(str(aaf)); f.write('\t')
            f.write(str(idx)); f.write('\t')
            f.write(str(avgaaf)); f.write('\t')
            f.write(str(avgloc)); f.write('\t')
        except:
            pass
        # isafe
        try:
            tab=pd.read_csv(isafefile,sep='\t')
            try:
                sub = tab[tab['POS']==loc]
                idx = sub.index.tolist()[0] + 1
                pos = sub['POS'].tolist()[0]
                aaf = sub['DAF'].tolist()[0]
            except:
                pos=loc
                idx=0 # signals not in file
            avgaaf=tab.iloc[:top]['DAF'].mean()
            avgloc=tab.iloc[:top]['POS'].mean()
            f.write(str(idx)); f.write('\t')
            f.write(str(avgaaf)); f.write('\t')
            f.write(str(avgloc)); f.write('\t')
        except:
            pass
        # isweep inference
        try:
            g=open(isweepfile2,'r')
            g.readline()
            line=g.readline().strip().split('\t')
            g.close()
            f.write(selcoef); f.write('\t')
            f.write(line[2])
            f.write('\n')
        except:
            pass

# saving
f.close()
