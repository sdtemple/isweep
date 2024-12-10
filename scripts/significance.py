import argparse
from scipy.stats import norm
import numpy as np
import pandas as pd

# Set up the argument parser
parser = argparse.ArgumentParser(description='Process IBD, analytical, and simulation data.')

# Define required arguments with named options
parser.add_argument('--input_ibd_file', 
                    type=str, 
                    required=True, 
                    help='Input IBD scan file')
parser.add_argument('--input_analytical_file', 
                    type=str, 
                    required=True, 
                    help='Input analytical method file output')
parser.add_argument('--input_simulation_file', 
                    type=str, 
                    required=True, 
                    help='Input simulation method file output')
parser.add_argument('--output_modified_file', 
                    type=str, 
                    required=True, 
                    help='IBD scan file modified with summary statistics and significance thresholds')
parser.add_argument('--output_excess_file', 
                    type=str, 
                    required=True, 
                    help='Only those regions genome-wide significant')

# Parse the arguments
args = parser.parse_args()
ibd = args.input_ibd_file
analytical = args.input_analytical_file
simulation = args.input_simulation_file
out = args.output_modified_file
out2 = args.output_excess_file

# Process analytical file
the_dictionary = dict()
with open(analytical, 'r') as f:
    for line in f:
        key, val = line.strip().split('\t')
        the_dictionary[key] = float(val)

pvalue = the_dictionary['confidence-level:']
lwrbnd = the_dictionary['initial-lower-bound:']
uprbnd = the_dictionary['initial-upper-bound:']
rmean = the_dictionary['revised-mean:']
rstd = the_dictionary['revised-standard-deviation:']
dz = the_dictionary['upper-discrete-z:']
dp = norm.sf(dz)
cz = the_dictionary['upper-continuous-z:']
cp = norm.sf(cz)
draw = the_dictionary['upper-discrete-raw:']
craw = the_dictionary['upper-continuous-raw:']

# Process simulation file
sims = []
with open(simulation, 'r') as f:
    for line in f:
        sims.append(float(line.strip()))
simcut = np.quantile(sims, 1 - pvalue)
simraw = simcut * rstd + rmean
simp = norm.sf(simcut)

# Process IBD file
table = pd.read_csv(ibd, sep='\t')
table['PVALUE'] = norm.sf(table['Z'])
table['UPPER_ANALYTICAL'] = draw
table['Z_UPPER_ANALYTICAL'] = dz
table['GW_LEVEL_ANALYTICAL'] = dp
table['UPPER_SIMULATE'] = simraw
table['Z_UPPER_SIMULATE'] = simcut
table['GW_LEVEL_SIMULATE'] = simp
table['ADJ_MEAN'] = rmean
table['ADJ_STDDEV'] = rstd
table['INIT_LOWER_BOUND'] = lwrbnd
table['INIT_UPPER_BOUND'] = uprbnd
table['UPPER_CONTINUOUS'] = craw
table['Z_UPPER_CONTINUOUS'] = cz
table['GW_LEVEL_CONTINUOUS'] = cp
table['CONFLEVEL'] = pvalue
table.to_csv(out, sep='\t', header=True, index=False)

# Filtered subtable
subtable = table[table['COUNT'] > draw]
subtable.to_csv(out2, sep='\t', header=True, index=False)