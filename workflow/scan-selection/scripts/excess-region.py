# Combine windows in a region of excess IBD

import sys
filein,fileout,cm=sys.argv[1:]
cm=float(cm)

g=open(fileout,'w')
g.write('CHROM\tBPLEFT\tBPRIGHT\tCMLEFT\tCMRIGHT\tMINIBD\tMAXIBD\n')
try:
    with open(filein,'r') as f:
        f.readline()
        line=f.readline()
        row=line.strip().split('\t')
        prev=row
        towrite=prev
        maxct=int(float(prev[2]))
        minct=int(float(prev[2]))
        for line in f:
            row=line.strip().split('\t')
            if int(row[3]) == int(prev[3]):
                if float(row[1]) >= (float(prev[1])+cm):
                    # gap in chromosome
                    g.write(towrite[3]); g.write('\t')
                    g.write(towrite[0]); g.write('\t')
                    g.write(prev[0]); g.write('\t')
                    g.write(towrite[1]); g.write('\t')
                    g.write(prev[1]); g.write('\t')
                    g.write(str(minct)); g.write('\t')
                    g.write(str(maxct)); g.write('\n')
                    towrite=row
                    maxct=int(float(row[2]))
                    minct=int(float(row[2]))
            else:
                # chromosome changed
                g.write(towrite[3]); g.write('\t')
                g.write(towrite[0]); g.write('\t')
                g.write(prev[0]); g.write('\t')
                g.write(towrite[1]); g.write('\t')
                g.write(prev[1]); g.write('\t')
                g.write(str(minct)); g.write('\t')
                g.write(str(maxct)); g.write('\n')
                towrite=row
                maxct=int(float(row[2]))
                minct=int(float(row[2]))
            prev=row
            nextct= int(float(prev[2]))
            if nextct > maxct:
                maxct = nextct
            if nextct < minct:
                minc = nextct
    g.write(towrite[3]); g.write('\t')
    g.write(towrite[0]); g.write('\t')
    g.write(prev[0]); g.write('\t')
    g.write(towrite[1]); g.write('\t')
    g.write(prev[1]); g.write('\t')
    g.write(str(minct)); g.write('\t')
    g.write(str(maxct)); g.write('\n')
except IndexError:
    # this happens if there is no excess IBD
    # make a blank file
    pass
g.close()
