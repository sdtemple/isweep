# Scans IBD data for regions with excess IBD.
# Seth Temple, GitHub: sdtemple
# Last modified: 2024-07-18

import argparse
import pandas as pd
import numpy as np

def main():
    # Set up argparse
    parser = argparse.ArgumentParser(
        description="Scan IBD data for regions with excess IBD."
    )
    
    # Define arguments

    parser.add_argument(
        'study',
        type=str,
        help="Study name. The main folder."
    )

    parser.add_argument(
        '--folder', 
        type=str,
        default='ibdsegs/ibdends/scan', 
        help="Folder path to IBD data."
    )
    
    parser.add_argument(
        '--all_genome', 
        type=str,
        default='scan.ibd.tsv', 
        help="Output file for saving all IBD data."
    )
    
    parser.add_argument(
        '--excess_genome', 
        type=str,
        default='excess.ibd.tsv', 
        help="Output file for saving excess IBD data."
    )
    
    parser.add_argument(
        '--chrlow', 
        type=int,
        default=1, 
        help="(default: 1) Lowest chromosome number."
    )
    
    parser.add_argument(
        '--chrhigh', 
        type=int,
        default=22, 
        help="(default: 22) Highest chromosome number."
    )

    parser.add_argument(
        '--prefix',
        type=str,
        default='chr',
        help="(default: chr) Prefix of chromosome files."
    )

    parser.add_argument(
        '--suffix',
        type=str,
        default='.ibd.windowed.tsv.gz',
        help="(default: .ibd.windowed.tsv.gz) Suffix of chromosome files."
    )
    
    parser.add_argument(
        '--heuristic_cutoff', 
        type=float,
        default=4.0, 
        help="(default: 4.0) Heuristic cutoff value for extreme IBD rates. First pass."
    )
    
    parser.add_argument(
        '--outlier_cutoff', 
        type=float,
        default=3.0, 
        help="(default: 3.0) Heuristic cutoff value for outlier IBD rates. Second pass."
    )
    
    # Parse the arguments
    args = parser.parse_args()

    all_genome = f"{args.study}/{args.all_genome}" 
    excess_genome = f"{args.study}/{args.excess_genome}" 

    # Reading in data
    tab = pd.read_csv(f"{args.study}/{args.folder}/{args.prefix}{args.chrlow}{args.suffix}", sep='\t')
    tab['CHROM'] = args.chrlow
    for i in range(args.chrlow + 1, args.chrhigh + 1):
        tabnow = pd.read_csv(f"{args.study}/{args.folder}/{args.prefix}{i}{args.suffix}", sep='\t')
        tabnow['CHROM'] = i
        tab = pd.concat((tab, tabnow))
    
    # Saving all data
    tab.to_csv(all_genome, sep='\t', index=False)

    # Calculating excess IBD
    medi = np.quantile(tab['COUNT'], 0.5)
    stdv = tab['COUNT'].std()
    a = medi - stdv * args.outlier_cutoff
    b = medi + stdv * args.outlier_cutoff
    sub = tab[(tab['COUNT'] >= a) & (tab['COUNT'] <= b)]
    medi = np.quantile(sub['COUNT'], 0.5)
    stdv = sub['COUNT'].std()
    b = medi + stdv * args.heuristic_cutoff
    out = tab[tab['COUNT'] >= b]

    # Saving excess IBD data
    out.to_csv(excess_genome, sep='\t', index=False)

if __name__ == "__main__":
    main()
