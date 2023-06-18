import sys
filein,=sys.argv[1:]
with open(filein,'r') as f:
    for line in f:
        pass
vals=line.strip().split('\t')
print(vals[0])
