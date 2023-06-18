import sys
filein,=sys.argv[1:]
f=open(filein,'r')
line=f.readline().strip().split('\t')
print(line[0])
f.close()
