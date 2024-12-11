# Make a plot of estimated autocovariances.
# Seth Temple, GitHub: sdtemple
# Last modified: 2024-12-10
# Current bug: some unnamed space at last column

import argparse
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
pd.options.display.float_format = '{:.8f}'.format

def main():

    # Set up argparse
    parser = argparse.ArgumentParser(
        description="Make a histogram of some standardized statistic."
    )

    # Define arguments
    parser.add_argument(
        '--input_autocov_file', 
        type=str,
        required=True, 
        help="Input file with autocovariance data."
    )

    parser.add_argument(
        '--input_analytical_file', 
        type=str,
        required=True, 
        help="Input file with estimated theta."
    )

    parser.add_argument(
        '--output_prefix', 
        type=str,
        required=True, 
        help="Output file prefix for saving the plot."
    )

    parser.add_argument(
        '--theta_type', 
        type=str,
        default='estimated-theta:', 
        help="Estimator defined in output FWER analysis file."
    )
    
    parser.add_argument(
        '--title', 
        type=str,
        default=None, 
        help="(default: None) Title of the plot."
    )

    parser.add_argument(
        '--xlabel',
        type=str,
        default='cM distance',
        help="(default: cM distance) X-axis label."
    )

    parser.add_argument(
        '--ylabel',
        type=str,
        default='Autocovariance',
        help="(default: Autocovariance) Y-axis label."
    )

    parser.add_argument(
        '--yupp', 
        type=float,
        default=1.5, 
        help="(default: 1.5) Upper bound on y-scale (for density=True)."
    )

    parser.add_argument(
        '--color',
        type=str,
        default='tab:blue',
        help="(default: tab:blue) Choose bar colors."
    )

    parser.add_argument(
        '--fontsize',
        type=int,
        default=12,
        help="(default: 12) Choose font size."
    )
    
    parser.add_argument(
        '--width', 
        type=float,
        default=6.4, 
        help="(default: 6.4) Width of the plot."
    )
    
    parser.add_argument(
        '--height', 
        type=float,
        default=4.8, 
        help="(default: 4.8) Height of the plot."
    )

    # Parse the arguments
    args = parser.parse_args()

    plt.figure(figsize=(args.width, args.height))
    plt.tight_layout()
    plt.rc('font',size=args.fontsize)

    # plotting
    table = pd.read_csv(args.input_autocov_file,sep='\t')
    the_dictionary = dict()
    with open(args.input_analytical_file,'r') as f:
        for line in f:
            key,val = line.strip().split('\t')
            the_dictionary[key] = float(val)
    stepsize = the_dictionary['step-size-morgan:']
    theta = the_dictionary[args.theta_type]
    
    xs = list(table.columns)
    xs = xs[:-1]
    xs = [float(x)*100 for x in xs]
    for i in range(table.shape[0]):
        ys = table.iloc[i]
        ys = ys[:-1]
        plt.plot(xs,
                 ys,
                 color=args.color)
    ys = table.iloc[0]
    ys = ys[:-1]
    plt.plot(xs,ys,
             color=args.color,
             label='Specific chromosomes')

    def func(x, a):
        return np.exp(-a * x)

    # Generate x values
    ws = np.array([w * stepsize for w in range(len(xs))])
    zs = func(ws, theta)
    ws = ws * 100
    plt.plot(ws,
             zs,
             color='k',
             linewidth=3,
             label='Log linear model fit')
    plt.legend(loc='upper right',title=None)

    plt.xlim(0, max(ws))
    plt.title(args.title, loc='center')
    plt.ticklabel_format(axis='y', style='sci', scilimits=(0,0))


    plt.ylim(-0.1,args.yupp)
    plt.grid(alpha=0.25)
    plt.xlabel(args.xlabel)
    plt.ylabel(args.ylabel)

    # Saving plots
    for pic in ['jpeg', 'png', 'tiff']:
        plt.savefig(f"{args.output_prefix}.{pic}")

if __name__ == '__main__':
    main()