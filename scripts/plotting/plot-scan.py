# Make a Manhattan plot of IBD rates.
# Seth Temple, GitHub: sdtemple
# Last modified: 2024-07-18
# Current bug: Must have a chromosome labeled 1 contiguous up to last chromosome.

import argparse
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

def main():
    # Set up argparse
    parser = argparse.ArgumentParser(
        description="Make a Manhattan plot of IBD rates."
    )
    
    # Define arguments
    parser.add_argument(
        '--input_file', 
        type=str,
        required=True, 
        help="Input file with IBD data."
    )
    
    parser.add_argument(
        '--output_prefix', 
        type=str,
        required=True, 
        help="Output file prefix for saving the plot."
    )

    parser.add_argument(
        '--sample_size', 
        type=int,
        required=True, 
        help="Sample size. Try 1 if not dividing by a count variable."
    )

    parser.add_argument(
        '--ploidy', 
        type=int,
        default=2, 
        help="(default: 2) Ploidy number"
    )

    parser.add_argument(
        '--sim_cutoff', 
        type=float,
        default=None, 
        help="(default: None) Simulation cutoff value for extreme IBD rates."
    )

    parser.add_argument(
        '--analytical_cutoff', 
        type=float,
        default=None, 
        help="(default: None) Analytical cutoff value for extreme IBD rates."
    )
    
    parser.add_argument(
        '--outlier_cutoff', 
        type=float,
        default=3.0, 
        help="(default: 3.0) Heuristic cutoff value for outlier IBD rates. First pass."
    )
    
    parser.add_argument(
        '--heuristic_cutoff', 
        type=float,
        default=4.0, 
        help="(default: 4.0) Heuristic cutoff value for extreme IBD rates. Second pass."
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
        default='Base pair along genome',
        help="(default: Base pair) X-axis label."
    )

    parser.add_argument(
        '--ylabel',
        type=str,
        default='IBD rate',
        help="(default: IBD rate) Y-axis label."
    )

    parser.add_argument(
        '--rotation', 
        type=int,
        default=90, 
        help="(default: 90) Rotation parameter in matplotlib."
    )

    parser.add_argument(
        '--yupp', 
        type=float,
        default=None, 
        help="(default: None) Upper bound on y-scale."
    )

    parser.add_argument(
        '--alpha',
        type=float,
        default=0.25,
        help="(0. for no grid lines) Intensity of grid lines."
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
        default=12., 
        help="(default: 12.0) Width of the plot."
    )
    
    parser.add_argument(
        '--height', 
        type=float,
        default=4.0, 
        help="(default: 4.0) Height of the plot."
    )
    
    # Parse the arguments
    args = parser.parse_args()

    plt.rc('font',size=args.fontsize)
    
    # Reading in data
    ibd = pd.read_table(args.input_file, sep='\t')

    numsamples = args.sample_size
    m = numsamples * args.ploidy
    if args.ploidy == 2:
        M =  m * (m-1) / 2 - numsamples
    elif args.ploidy == 1:
        M =  m * (m-1) / 2
    else:
        M = 1

    statistic = args.statistic


    # Maths
    medi = ibd[statistic].median()
    stdv = ibd[statistic].std()
    a = medi - stdv * args.outlier_cutoff
    b = medi + stdv * args.outlier_cutoff
    low = a
    newibd = ibd[(ibd[statistic] >= a) & (ibd[statistic] <= b)].copy()

    # Compute browning and browning 2020 cutoff
    medi = newibd[statistic].median(numeric_only=True)
    stdv = newibd[statistic].std(numeric_only=True)
    upp = medi + stdv * args.heuristic_cutoff
    md = medi
    md = md / M
    upp = upp / M

    # Chromosome maths
    ibd["CUMPOS"] = None
    s = 0
    mbp = [0] * (args.chr_high + 1)
    chr = {i: i for i in range(args.chr_low, args.chr_high + 1)}
    nchr = len(chr)

    for i in range(args.chr_low, args.chr_high + 1):
        try:
            mbp[i] = ibd.loc[ibd["CHROM"] == chr[i], "CMWINDOW"].max()
            ibd.loc[ibd["CHROM"] == chr[i], "CUMPOS"] = ibd.loc[ibd["CHROM"] == chr[i], "CMWINDOW"] + s
            if not np.isnan(mbp[i]):
                s = s + mbp[i]
        except:
            pass
    ibd["CHROM"] = ibd["CHROM"].astype(int)

    # Plotting
    plt.figure(figsize=(args.width, args.height))
    mxy = ibd[statistic].max() * 1.5
    mxy = mxy / M
    for chrom, group in ibd.groupby("CHROM"):
        clr = "black" if chrom % 2 == 0 else "gray"
        plt.plot(group["CUMPOS"], group[statistic]/M, c=clr)
        
    plt.axhline(y=md, color="tab:blue", linestyle="solid", label="Median", linewidth=2)
    plt.axhline(y=upp, color="tab:orange", linestyle='dashed', label='Heuristic',alpha=0.75)

    if args.analytical_cutoff is not None:
        ac = args.analytical_cutoff / M
        plt.axhline(y=ac, color="tab:green", linestyle=(0,(5,10)), label="Analytical",linewidth=1,alpha=0.75)
    if args.sim_cutoff is not None:
        sc = args.sim_cutoff / M
        plt.axhline(y=sc, color="tab:red", linestyle='dotted', label='Simulation',linewidth=2,alpha=0.75)

    plt.xlabel(args.xlabel,loc='left')
    plt.ylabel(args.ylabel)

    if args.yupp is not None:
        plt.ylim(0, args.yupp)
    else:
        plt.ylim(0, mxy)

    plt.xlim(ibd['CUMPOS'].min() - 25, ibd['CUMPOS'].max() + 25)
    plt.title(args.title, loc='center')
    plt.tight_layout()
    plt.legend(loc="upper right", title=None)
    plt.ticklabel_format(axis='y', style='sci', scilimits=(0,0))

    # Customize x-axis ticks
    chromosome_values = sorted(ibd["CHROM"].unique())
    x_ticks = [np.mean(ibd[ibd["CHROM"] == chrom]["CUMPOS"]) for chrom in chromosome_values]
    x_labels = [f"{chrom}" for chrom in chromosome_values]
    plt.xticks(x_ticks, x_labels, rotation=args.rotation)

    plt.grid(alpha=args.alpha)

    # Saving plots
    for pic in ['jpeg', 'png', 'tiff']:
        plt.savefig(f"{args.output_prefix}.{pic}")

if __name__ == "__main__":
    main()
