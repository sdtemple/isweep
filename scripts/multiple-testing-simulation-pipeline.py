# Simulating Ornstein-Uhlenbeck process
# For multiple testing correction
# Seth Temple, GitHub: sdtemple
# Last modified: 2024-07-18
# Description: This script is used to estimate a multiple testing correction,
# using simulation of an Ornstein-Uhlenbeck process.

# Example usage:
# python ibd-selscan-simtest.py 0.05 35.50 22 1.5 0.0005 1000

import argparse
import numpy as np

def main(input_file, output_file, numsims):

    the_dictionary = dict()
    with open(input_file,'r') as f:
        for line in f:
            key,val = line.strip().split('\t')
            the_dictionary[key] = float(val)
    stepsize = the_dictionary['step-size-morgan:']
    chrnum = the_dictionary['chromosome-number:']
    chrlen = the_dictionary['average-chromosome-length-morgan:']
    theta = the_dictionary['estimated-theta:']
    pvalue = the_dictionary['confidence-level:']

    # Calculate K
    K = int(np.floor(chrlen / stepsize) * chrnum)
    
    # Simulate Ornstein-Uhlenbeck process
    maxs = []
    for j in range(1, numsims + 1):
        # if j % 100 == 0:
        #     print(j)
        x = np.random.normal()
        xs = [x]
        for i in range(2, K + 1):
            newx = np.random.normal(loc=x * np.exp(-stepsize * theta), scale=np.sqrt(2 - 2 * np.exp(-stepsize * theta)))
            xs.append(newx)
            x = newx
        xs = np.array(xs)
        zs = (xs - np.mean(xs)) / np.std(xs)
        maxs.append(np.max(zs))
    
    # Compute multiple testing correction
    # I.e., find the 1 - pvalue quantile of the maximums
    maxs = np.array(maxs)
    if len(maxs) > 0:
        out = np.quantile(maxs, 1 - pvalue)
    else:
        out = -100
        maxs = [-100]
    print(out)

    with open(output_file,'w') as f:
        for m in maxs:
            f.write(str(m))
            f.write('\n')

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Estimate a multiple testing correction using simulation of an Ornstein-Uhlenbeck process.')
    
    parser.add_argument('--input_file',
                        type=str,
                        required=True,
                        help='File output from multiple-testing-analytical.py'
                        )
    parser.add_argument('--output_file',
                        type=str,
                        required=True,
                        help='Output file with maxima of simulations'
                        )
    parser.add_argument('--num_sims', 
                        type=int, 
                        default=1000, 
                        help='(default: 1000) Number of simulations')
    
    args = parser.parse_args()
    
    main(args.input_file,
         args.output_file, 
         args.num_sims)
