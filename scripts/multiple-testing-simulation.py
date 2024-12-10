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

def main(pvalue, theta, chrnum, chrlen, stepsize, numsims):
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
    out = np.quantile(maxs, 1 - pvalue)
    print(out)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Estimate a multiple testing correction using simulation of an Ornstein-Uhlenbeck process.')

    parser.add_argument('--output_file',
                        type=str,
                        help='Name of file with simulated Zs'
                        )
    
    parser.add_argument('--confidence_level', 
                        type=float, 
                        default=0.05, 
                        help='(default: 0.05) p-value for the multiple testing correction')
    parser.add_argument('--theta', 
                        type=float, 
                        default=35, 
                        help='(default: 35) Theta parameter for Ornstein-Uhlenbeck process')
    parser.add_argument('--chromosome_number', 
                        type=int, 
                        default=22, 
                        help='(default: 22) Chromosome number')
    parser.add_argument('--chr_average_size', 
                        type=float, 
                        default=1.5, 
                        help='(default: 1.5) Average chromosome length (in Morgans)')
    parser.add_argument('--cM_step_size', 
                        type=float, 
                        default=0.0005, 
                        help='(default: 0.0005) Step size for each test (in Morgans)')
    parser.add_argument('--num_sims', 
                        type=int, 
                        default=1000, 
                        help='(default: 1000) Number of simulations')
    
    args = parser.parse_args()
    
    main(args.output_file,
         args.confidence_level, 
         args.theta, 
         args.chromosome_number, 
         args.chr_average_size, 
         args.cM_step_size, 
         args.num_sims)
