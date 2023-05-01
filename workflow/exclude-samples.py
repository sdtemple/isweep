import sys
sampfile,subsampfile,exclsampfile=sys.argv[1:]
samplist=[]
with open(sampfile,'r') as f:
    for line in f:
        samplist.append(line.strip())
subsamplist=[]
with open(subsampfile,'r') as f:
    for line in f:
        subsamplist.append(line.strip())
exclsamplist=[]
for samp in samplist:
    if samp in subsamplist:
        pass
    else:
        exclsamplist.append(samp)
with open(exclsampfile,'w') as f:
    for exclsamp in exclsamplist:
        f.write(str(exclsamp))
        f.write('\n')
