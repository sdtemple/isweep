import sys
file1,file2,fileout=sys.argv[1:]
lst1 = []
with open(file1,'r') as f:
	for line in f:
		lst1.append(line.strip())
lst2 = []
with open(file2,'r') as f:
	for line in f:
		lst2.append(line.strip())
set1 = set(lst1)
set2 = set(lst2)
setint = set1.intersection(set2)
f=open(fileout,'w')
for i in setint:
	f.write(str(i))
	f.write('\n')
