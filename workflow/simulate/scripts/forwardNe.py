# Format Ne for forward simulations

import sys
from math import floor

filein, fileout = sys.argv[1:]

nes = [] # vector of effective population sizes

# read *.ne file
prior = 0
with open(filein, 'r') as f:
    line = f.readline()
    for line in f:
        g, ne = line.split('\t')[:2] # tab-separated
        try:
            ne1, ne2 = ne.split('E') # parse scientific notation
            ne1 = float(ne1)
            ne2 = int(ne2)
            ne = ne1 * (10 ** ne2)
        except ValueError:
            ne = int(float(ne))
        ne = floor(ne)
        while prior < int(g):
            prior += 1
            nes.append(ne)
        nes.append(ne)
        prior += 1

# for forward simulation
nes.reverse()

# write to file
f = open(fileout, 'w')
for ne in nes:
    f.write(str(ne))
    f.write('\n')
f.close()
