import sys
import numpy as np
import pandas as pd

### script inputs in order
### see https://github.com/mfumagalli/ImaGene/blob/master/params.txt
# file name to write over
# location of msms jar file
# folder for sims
# N_e value for scaling
# *.ne file for demography
# mutation rate
# recombination rate
# number of basepairs in genome
# sample size for haploids
# lowest selection coefficient, e.g. s=0
# highest selection coefficient, e.g. s=0.05
# selection coefficient step, e.g. s=0.001
# frequency of selected allele at start, e.g. 0.001
# time of selection, e.g. 300 generations
# number of replicates per selection coefficient
# number of batches
# number of threads

fileout = sys.argv[1]
jar_loc = sys.argv[2]
sim_loc = sys.argv[3]
Ne_scale = int(float(sys.argv[4]))

# non-selection parameters
ne_loc = sys.argv[5]
mut_rate = float(sys.argv[6])
recomb_rate = float(sys.argv[7])
num_bp = int(float(sys.argv[8]))
sample_size = int(float(sys.argv[9]))

# selection parameters
sela = int(float(sys.argv[10]) * 2 * Ne_scale)
selb = int(float(sys.argv[11]) * 2 * Ne_scale)
selstep = int(float(sys.argv[12]) * 2 * Ne_scale)
selfreq=float(sys.argv[13])
seltime = int(float(sys.argv[14]))
seltimeden = int(float(Ne_scale * 4))

# simulation paramters
nrepl = int(float(sys.argv[15]))
nbatch = int(float(sys.argv[16]))
nthreads = int(float(sys.argv[17]))

# scaling
theta = Ne_scale * 4 * num_bp * mut_rate
rho = 4 * Ne_scale * recomb_rate * (num_bp - 1)
Ne_table = pd.read_csv(ne_loc, sep='\t')
Y = np.round(Ne_table['NE'] / Ne_scale, 4)
X = Ne_table['GEN'] / Ne_scale / 4
J = len(X)
demstring = ''
for j in range(J):
    x = X[j]
    y = Y[j]
    addstring = '-eN ' + str(x) + ' ' + str(y) + ' '
    demstring += addstring


### writing ###
f = open(fileout, 'w')

# directories
f.write("### 1) DIRECTORIES \n\n")
f.write("DIRMSMS="+str(jar_loc))
f.write('\n')
f.write("DIRDATA="+str(sim_loc))
f.write('\n\n')

# demographic model
f.write("### 2) DEMOGRAPHIC MODEL \n\n")
f.write("NREF="+str(int(Ne_scale)))
f.write('\n')
f.write("DEMO='" + demstring + "'")
f.write('\n\n')

# locus and sample size
f.write("### 3) LOCUS AND SAMPLE SIZE \n\n")
f.write("LEN="+str(num_bp))
f.write('\n')
f.write("THETA="+str(theta))
f.write('\n')
f.write("RHO="+str(rho))
f.write('\n')
f.write("NCHROMS="+str(int(sample_size)))
f.write('\n\n')

# selection
f.write("### 4) SELECTION\n\n")
f.write("SELPOS=`bc <<< 'scale=2; 1/2'`")
f.write('\n')
f.write("FREQ=`bc <<< 'scale=6; " + str(selfreq) + "'`")
f.write('\n')
f.write("SELRANGE=`seq " + str(sela) + ' ' + str(selstep) + ' ' + str(selb) + '`')
f.write('\n')
f.write("NREPL="+str(nrepl))
f.write('\n')
f.write("TIMERANGE=`bc <<< 'scale=4; " + str(seltime) + '/' + str(seltimeden) + "'`")
f.write('\n\n')

# sim details
f.write("### 5) Simulation study settings\n\n")
f.write("NBATCH="+str(nbatch))
f.write('\n')
f.write("NTHREADS="+str(nthreads))
f.close()