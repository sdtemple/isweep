import sys
from isweep import *

fileout = sys.argv[1]
sample_size = sys.argv[2]
Ne_file = sys.argv[3]
ibd_length = sys.argv[4]
num_sim = sys.argv[5]
num_rep = sys.argv[6]

num_sim = int(float(num_sim))
num_rep = int(float(num_rep))
sample_size = int(float(sample_size))
ibd_length = float(ibd_length)

print('number of simulations')
print(str(num_sim))
print('\n')
print('number of replicate')
print(str(num_rep))
print('\n')
print('sample size')
print(str(sample_size))
print('\n')
print('ibd length cutoff')
print(str(ibd_length))
print('\n')
print('simulation iterator')

Ne = read_Ne(Ne_file)

with open(fileout,'w') as f:

	# number of simulations for calculating type 1 error
	for sim in range(num_sim):
	
		# print statement
		if sim % 10 == 0:
			print(sim)

		# number of replicates to do shapiro-wilks test
		for rep in range(num_rep-1):
			ibd_out = simulate_ibd(sample_size, Ne, ibd_length, ibd_length, 2, False, False)
			num_tracts = ibd_out[0]
			f.write(str(num_tracts))
			f.write('\t')
		
		# last replicate
		ibd_out = simulate_ibd(sample_size, Ne, ibd_length, ibd_length, 2, False, False)
		num_tracts = ibd_out[0]
		f.write(str(num_tracts))
		f.write('\n')

