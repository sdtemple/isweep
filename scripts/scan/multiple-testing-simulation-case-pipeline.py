#!/bin/python
import argparse
import numpy as np
from scipy.stats import pearsonr
from math import floor

def main():

    parser = argparse.ArgumentParser(
        description='Simulation-based threshold for 2-dim Ornstein-Uhlenbeck process (cases and controls)'
    )

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

    the_dictionary = dict()
    with open(args.input_file,'r') as f:
        for line in f:
            key,val = line.strip().split('\t')
            the_dictionary[key] = float(val)
    stepsize = the_dictionary['step-size-morgan:']
    chrnum = the_dictionary['chromosome-number:']
    chrlen = the_dictionary['average-chromosome-length-morgan:']
    theta1 = the_dictionary['estimated-theta1:']
    theta0 = the_dictionary['estimated-theta0:']
    rho = the_dictionary['estimated-rho:']
    pvalue = the_dictionary['confidence-level:']


    # Parameters
    mu_X, mu_Y = 0.0, 0.0
    sigma_X, sigma_Y = 1., 1.
    
    J = args.num_sims
    
    dt = stepsize
    N = int(floor(chrlen/stepsize)*chrnum)
    
    # Initialize arrays
    X = np.zeros(N)
    Y = np.zeros(N)
    X[0], Y[0] = 0.0, 0.0  # Initial conditions
    
    # Cholesky decomposition for correlated Brownian motions
    L = np.linalg.cholesky([[1, rho], [rho, 1]])
    
    
    f = open(args.output_file,'w')
    
    # Simulate the process many times
    for j in range(J):
        # Simulate the process once
        for t in range(1, N):
            dW = np.sqrt(dt) * np.random.normal(size=2)
            dW_X, dW_Y = L @ dW
            X[t] = X[t-1] + theta1 * (mu_X - X[t-1]) * dt + sigma_X * dW_X
            Y[t] = Y[t-1] + theta0 * (mu_Y - Y[t-1]) * dt + sigma_Y * dW_Y
        Yz = (Y - np.mean(Y)) / np.std(Y)
        Xz = (X - np.mean(X)) / np.std(X)
        Z = Yz - Xz
        Zz = ( Z - np.mean(Z) ) / np.std(Z)
        f.write(str(max(Zz))); f.write('\n')
    
    f.close()

if __name__ == '__main__':
    main()
    
