#!/bin/python
import argparse
import numpy as np
from scipy.stats import pearsonr
from math import floor

def main():

    parser = argparse.ArgumentParser(
        description='Simulation-based threshold for 2-dim Ornstein-Uhlenbeck process (cases and controls)'
    )

    parser.add_argument('--output_file',
                        type=str,
                        help='Name of file with simulated Zs'
                        )
    parser.add_argument('--confidence_level', 
                        type=float,
                        default=0.05, 
                        help='(default: 0.05) p-value for the multiple testing correction')
    parser.add_argument('--theta1', 
                        type=float, 
                        default=35, 
                        help='(default: 35) Theta parameter for OU process in sample cases')
    parser.add_argument('--theta0', 
                        type=float, 
                        default=35, 
                        help='(default: 35) Theta parameter for OU process in sample controls')
    parser.add_argument('--rho', 
                        type=float, 
                        default=0., 
                        help='(default: 0.) Rho parameter for cross-correlated OU process')
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


    # Parameters
    theta1 = args.theta1
    theta0 = args.theta0
    mu_X, mu_Y = 0.0, 0.0
    sigma_X, sigma_Y = 1., 1.
    rho = args.rho  # Correlation coefficient
    
    J = args.num_sims
    
    stepsize = args.cM_step_size
    dt = stepsize
    chrnum = args.chromosome_number
    chrlen = args.chr_average_size
    N = floor(chrlen/stepsize)*chrnum
    
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
    
