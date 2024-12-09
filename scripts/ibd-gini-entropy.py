# compute ibd entropy

# setup
import argparse
import numpy as np

# Set up the argument parser
parser = argparse.ArgumentParser(description='Compute IBD entropy from outlier files.')

parser.add_argument('--input_folder', 
                    type=str,
                    required=True, 
                    help='Path to the folder containing outlier files')
parser.add_argument('--output_file', 
                    type=str,
                    required=True, 
                    help='Path to the output file')
parser.add_argument('--sample_size', 
                    type=float,
                    required=True, 
                    help='Sample size')

# Parse the arguments
args = parser.parse_args()
folder = args.input_folder
fileout = args.output_file
samplesize = args.sample_size

# read outlier files
ctr = 0  
numnodes = []
while True:
    try:
        ctr += 1
        numnode = 0
        with open(f"{folder}/outlier{ctr}.txt", 'r') as f:
            for line in f:
                numnode += 1
        numnodes.append(numnode)
    except FileNotFoundError:
        break

# calculate entropy
numnodes = np.array(numnodes)
propnodes = numnodes / sum(numnodes)
entropy = -sum(propnodes * np.log10(propnodes))
gini = 1 - sum(propnodes * propnodes)
print(entropy)

# write statistics
with open(fileout, 'w') as f:
    f.write('entropy\tmax_prop\tmax_freq\tmean_freq\tnum_group\tgini\n')
    statistics = [
        entropy,
        max(propnodes),
        max(numnodes) / samplesize,
        sum(numnodes) / samplesize,
        len(propnodes),
        gini,
    ]
    for statistic in statistics[:-1]:
        f.write(f"{statistic}\t")
    f.write(f"{statistics[-1]}\n")