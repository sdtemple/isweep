# Selection coefficient estimation using IBD data
# Using normality-based confidence intervals
# Seth Temple, GitHub: sdtemple
# Converted to argparse on 2024-07-18

import argparse
from isweep import *
import pandas as pd
import numpy as np
import gzip

def main():
    # Set up argparse
    parser = argparse.ArgumentParser(
        description="Estimate the selection coefficient using IBD data (Gaussian CIs)."
    )
    
    # Define arguments
    parser.add_argument(
        '--output_file', 
        type=str,
        required=True, 
        help="Output file to save the results."
    )

    parser.add_argument(
        '--sample_size', 
        type=int,
        required=True, 
        help="Number of samples."
    )
    
    parser.add_argument(
        '--Ne_est', 
        type=str,
        required=True, 
        help="Estimated effective population size file."
    )

    parser.add_argument(
        '--p_est', 
        type=float,
        required=True, 
        help="Estimated allele frequency."
    )
    
    parser.add_argument(
        '--ibd_count', 
        type=int,
        required=True, 
        help="The count of detected IBD segments."
    )
    
    parser.add_argument(
        '--ibd_cutoff', 
        type=float,
        default=3., 
        help="(default: 3.0) The cutoff value for IBD sgements."
    )
    
    parser.add_argument(
        '--model', 
        type=str,
        default='a', 
        help="(default: a for additive) Genetic model for positive selection."
    )
    
    parser.add_argument(
        '--alpha', 
        type=float,
        default=0.05, 
        help="(default: 0.05) Significance level for confidence intervals."
    )
    
    parser.add_argument(
        '--ploidy', 
        type=int,
        default=2, 
        help="(default: 2) Ploidy level."
    )

    parser.add_argument(
        '--num_bootstraps', 
        type=int,
        default=100, 
        help="(default: 100) Number of bootstrap samples."
    )
    
    # Parse the arguments
    args = parser.parse_args()
    
    alpha = args.alpha
    alpha1 = (1 - alpha) / 2
    alpha2 = 1 - alpha1
    inhs = [args.model]

    # Setting up
    B = args.num_bootstraps
    CUTOFF = args.ibd_cutoff
    long_ibd = CUTOFF
    Ne = read_Ne(args.Ne_est)
    ab = [CUTOFF, np.inf]
    n = args.sample_size
    ploidy = args.ploidy
    m = ploidy * n
    N = m * (m - 1) / 2 - m
    p_est = args.p_est
    numTracts = int(float(args.ibd_count))

    sinhs = []
    for inh in inhs:
        s_est = minimize_scalar(
            chi2_isweep,
            args=(p_est, Ne, N, (numTracts,), ab, inh),
            bounds=(0, 0.5),
            method='bounded'
        ).x
        sinhs.append(s_est)

    # Bootstrap
    sbsinhs = [[] for _ in range(len(sinhs))]

    for j in range(len(inhs)):
        sinh = sinhs[j]
        inh = inhs[j]
        for b in range(B):
            simdata = simulate_ibd_isweep(
                n, sinh, p_est, Ne,
                long_ibd=long_ibd,
                short_ibd=long_ibd,
                ploidy=ploidy,
                one_step_model=inh
            )
            simibd = simdata[0][0]
            sb = minimize_scalar(
                chi2_isweep,
                args=(p_est, Ne, N, (simibd,), ab, inh),
                bounds=(0, 0.5),
                method='bounded'
            ).x
            sbsinhs[j].append(sb)

    # Writing results
    with open(args.output_file, 'w') as f:
        f.write('VarFreqEst\t')
        f.write('SelCoefEst\t')
        f.write('SelCoefLow\t')
        f.write('SelCoefUpp\t')
        f.write('Model\n')

        for j in range(len(inhs)):
            sinh = sinhs[j]
            sbsinh = sbsinhs[j]
            inh = inhs[j]
            sl, sm, su = bootstrap_standard(sinh, sbsinh, alpha1, alpha2)
            f.write(f"{p_est}\t{sm}\t{sl}\t{su}\t{inh}\n")

if __name__ == "__main__":
    main()
