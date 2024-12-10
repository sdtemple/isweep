import sys
filein,row,col=sys.argv[1:]
row=int(float(row))
col=int(float(col))
with open(filein,'r') as f:
    for i in range(row):
        line=f.readline().strip().split('\t')
    val=line[col-1]
print(val)
