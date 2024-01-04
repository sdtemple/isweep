# compute ibd entropy for haplotype groups

# setup
import sys
import numpy as np
folder, fileout, samplesize, = sys.argv[1:]
samplesize=float(samplesize)

# read outlier files
ctr= 0 
numnodes = []
try:
    while True:
        ctr += 1
        numnode = 0
        with open(folder + '/outlier' + str(ctr) + '.txt', 'r') as f:
            for line in f:
                numnode += 1
        numnodes.append(numnode)
except:
    pass

# calculate entropy
numnodes = np.array(numnodes)
propnodes = numnodes / sum(numnodes)
entropy = - sum(propnodes * np.log10(propnodes))
print(entropy)

# write statistics
f = open(fileout, 'w')
f.write('entropy\tmax_prop\tmax_freq\tmean_freq\tnum_group\n')
statistics = [entropy,
              max(propnodes),
              max(numnodes)/samplesize,
              sum(numnodes)/samplesize,
              len(propnodes)
              ]
for statistic in statistics[:-1]:
    f.write(str(statistic) + '\t')
f.write(str(statistics[-1]) + '\n')
f.close()