import sys
import gzip
filein,fileout,headersize,denom=sys.argv[1:]
headersize=int(float(headersize))
denom=float(denom)
h=open(fileout,'w')
ctr = 0
with gzip.open(filein, 'rt') as f:
    for i in range(headersize):
        f.readline()
    for line in f:
        ctr += 1
        lindat = line.strip().split('\t')
        rsid = 'rs' + str(ctr)
        bppos = int(float(lindat[1]))
        chrpos = lindat[0]
        cmpos = bppos / denom
        bppos = str(bppos)
        cmpos = str(cmpos)
        h.write(chrpos + '\t')
        h.write(rsid + '\t')
        h.write(bppos + '\t')
        h.write(cmpos + '\n')
h.close()