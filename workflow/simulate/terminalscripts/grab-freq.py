import sys
filein,=sys.argv[1:]
f=open(filein,'r')
f.readline().strip()
print(f.readline().strip())
f.close()
