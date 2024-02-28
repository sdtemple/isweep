import sys
from isweep import *

fileout = sys.argv[1] # output file prefix
sample_size = sys.argv[2] # diploid sample size
Ne_file = sys.argv[3]  # file with Ne dictionary values
ibd_length = sys.argv[4] # IBD length cutoff
num_sim = sys.argv[5] # number of simulations (these are for histogram of each replicate)
num_rep = sys.argv[6] # number of replicates (each is an experiment)

num_sim = int(float(num_sim))
num_rep = int(float(num_rep))
sample_size = int(float(sample_size))
ibd_length = float(ibd_length)

print('number of replicates')
print(str(num_rep))
print('\n')
print('number of simulations')
print(str(num_sim))
print('\n')
print('sample size')
print(str(sample_size))
print('\n')
print('ibd length cutoff')
print(str(ibd_length))
print('\n')
print('simulation iterator')

Ne = read_Ne(Ne_file)

fileout=fileout+'-'+Ne_file+'-size'+str(sample_size)+'-cM'+str(ibd_length)+'-sim'+str(num_sim)+'-rep'+str(num_rep)+'.txt'
with open(fileout,'w') as f:
	for sim in range(num_sim):
		if sim % 10 == 0:
			print(sim)
		for rep in range(num_rep-1):
			ibd_out = simulate_ibd(sample_size, Ne, ibd_length, ibd_length, 2, False, False)
			num_tracts = ibd_out[0]
			f.write(str(num_tracts))
			f.write('\t')
		ibd_out = simulate_ibd(sample_size, Ne, ibd_length, ibd_length, 2, False, False)
		num_tracts = ibd_out[0]
		f.write(str(num_tracts))
		f.write('\n')
