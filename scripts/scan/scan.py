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
        '--input_study',
        type=str,
        help="Study name. The main folder."
    )

    parser.add_argument(
        '--input_folder', 
        type=str,
        default='ibdsegs/ibdends/scan', 
        help="Folder path to IBD data."
    )
    
    parser.add_argument(
        '--output_all_genome_file', 
        type=str,
        default='scan.ibd.tsv', 
        help="Output file for saving all IBD data."
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
        '--input_prefix',
        type=str,
        default='chr',
        help="(default: chr) Prefix of chromosome files."
    )

    parser.add_argument(
        '--input_suffix',
        type=str,
        default='.ibd.windowed.tsv.gz',
        help="(default: .ibd.windowed.tsv.gz) Suffix of chromosome files."
    )
    
    parser.add_argument(
        '--outlier_cutoff', 
        type=float,
        default=4.0, 
        help="(default: 4.0) Heuristic cutoff value for outlier IBD rates. Second pass."
    )
    
    # Parse the arguments
    args = parser.parse_args()

    all_genome = f"{args.input_study}/{args.output_all_genome_file}"

    # Reading in data
    tab = pd.read_csv(f"{args.input_study}/{args.input_folder}/{args.input_prefix}{args.chr_low}{args.input_suffix}", sep='\t')
    tab['CHROM'] = args.chr_low
    for i in range(args.chr_low + 1, args.chr_high + 1):
        try:
            tabnow = pd.read_csv(f"{args.input_study}/{args.input_folder}/{args.input_prefix}{i}{args.input_suffix}", sep='\t')
            tabnow['CHROM'] = i
            tab = pd.concat((tab, tabnow))
        except:
            pass

    # Calculating excess IBD
    medi = np.median(tab['COUNT'])
    # medi = np.quantile(tab['COUNT'], 0.5)
    stdv = tab['COUNT'].std()
    a = medi - stdv * args.outlier_cutoff
    b = medi + stdv * args.outlier_cutoff
    sub = tab[(tab['COUNT'] >= a) & (tab['COUNT'] <= b)]
    # medi = np.quantile(sub['COUNT'], 0.5)
    avg = np.mean(sub['COUNT'])
    stdv = sub['COUNT'].std()

    # Saving all data
    tab['Z'] = (tab['COUNT'] - avg) / stdv
    tab.to_csv(all_genome, sep='\t', index=False)

if __name__ == "__main__":
    main()
