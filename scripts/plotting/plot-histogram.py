# Make a histogram of IBD rates.
# Seth Temple, GitHub: sdtemple
# Last modified: 2024-12-10

import pandas as pd
import matplotlib.pyplot as plt
import argparse
pd.options.display.float_format = '{:.8f}'.format

def main():

    # Set up argparse
    parser = argparse.ArgumentParser(
        description="Make a histogram of some standardized statistic."
    )

    # Define arguments
    parser.add_argument(
        '--input_file', 
        type=str,
        required=True, 
        help="Input file with IBD data."
    )

    parser.add_argument(
        '--output_file', 
        type=str,
        required=True, 
        help="Output file prefix for saving the plot."
    )

    parser.add_argument(
        '--chr_low', 
        type=int,
        default=1, 
        help="(default: 1) Lowest chromosome number."
    )
    
    parser.add_argument(
        '--chr_high', 
        type=int,
        default=22, 
        help="(default: 22) Highest chromosome number."
    )

    parser.add_argument(
        '--chrom',
        type=str,
        default="CHROM",
        help="(default: CHROM) Name of chromosome numbers to use."
    )

    parser.add_argument(
        '--statistic',
        type=str,
        default="COUNT",
        help="(default: COUNT) Statistic to plot."
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
        default='Standardized statistic',
        help="(default: Standardized statistic) X-axis label."
    )

    parser.add_argument(
        '--ylabel',
        type=str,
        default='Density',
        help="(default: Histogram density) Y-axis label."
    )

    parser.add_argument(
        '--xupp', 
        type=float,
        default=4., 
        help="(default: 4.) Upper (and symmetric lower) bound on x-scale."
    )

    parser.add_argument(
        '--yupp', 
        type=float,
        default=None, 
        help="(default: None) Upper bound on y-scale (for density=True)."
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

    # read and process input
    statistic = args.statistic
    chrom = args.chrom
    table = pd.read_csv(args.input_file,sep='\t')
    subtable = table[(table[args.chrom] >= args.chr_low) & 
                     (table[args.chrom] <= args.chr_high)] # control which chromosomes
    y = subtable[statistic]
    z = (y - y.mean()) / y.std() # normalize

    plt.figure(figsize=(args.width, args.height))
    plt.tight_layout()
    plt.rc('font',size=args.fontsize)

    # plotting
    plt.hist(z,
             50,
             edgecolor='k',
             density=True,
             color=args.color,
             range=(-args.xupp,args.xupp)
             )
    
    # formatting
    if args.yupp is not None:
        plt.ylim(0,args.yupp)
    plt.grid(alpha=0.25)
    plt.xlabel(args.xlabel)
    plt.ylabel(args.ylabel)
    plt.title(args.title)
    plt.xlim(-args.xupp,args.xupp)

    # saving
    plt.savefig(args.output_file)

if __name__ == "__main__":
    main()