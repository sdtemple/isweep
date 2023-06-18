import sys
nein,headin,fileout=sys.argv[1:]

# Ne file
def read_Ne(file):
    '''Read *.ne file
    Parameters
    ----------
    file : string
        Input file name
    Returns
    -------
    dict
        dict[generation] = size
    '''
    Ne = {}
    f = open(file, 'r')
    f.readline() # header
    for line in f:
        g, size = line.split('\t')[:2]
        Ne[int(g)] = int(float(size))
    f.close()
    return Ne
Ne=read_Ne(nein)

# file management
f=open(headin,'r')
line=f.readline()

g=open(fileout,'w')
g.write(line)

line=f.readline()
g.write(line)

# make the coalescent times integers
times=line.strip().split(' ')
times=[int(float(t)) for t in times[1:]]

# max time in Ne file
maxg=max(Ne.keys())
maxn=Ne[maxg]

# times but last time
g.write('0 0 ')
for t in times[:-1]:
    try:
        co=0.5/Ne[t]
    # if beyond the recent Ne file
    except KeyError:
        co=0.5/maxn
    g.write(str(co))
    g.write(' ')

# last time
try:
    co=0.5/Ne[times[-1]]
except KeyError:
    co=0.5/maxn

g.write(str(co))
g.write('\n')
g.close()
