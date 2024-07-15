import sys
filein, fileout, validx = sys.argv[1:]
validx = int(validx)
g = open(fileout, 'w')
with open(filein, 'r') as f:
    f.readline()
    f.readline()
    f.readline()
    for line in f:
        vals = line.strip().split('\t')
        vals = [val.split(',')[validx] for val in vals]
        for v in vals[:-1]:
            g.write(v)
            g.write('\t')
        g.write(vals[-1])
        g.write('\n')
g.close()