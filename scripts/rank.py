# Score variants based on differences between IBD clusters.
# Seth Temple, GitHub: sdtemple
# Converted to argparse on 2024-07-18

import argparse
import numpy as np
from isweep import *

def main():
    # Set up argparse
    parser = argparse.ArgumentParser(
        description="Score variants based on differences between IBD clusters."
    )
    
    # Define arguments
    parser.add_argument(
        '--ibd_file', 
        type=str, 
        help="Path to the input IBD file."
    )
    
    parser.add_argument(
        '--vcf', 
        type=str, 
        help="Path to the input VCF file."
    )
    
    parser.add_argument(
        '--file_out', 
        type=str, 
        help="Output file to save the scored variants."
    )
    
    parser.add_argument(
        '--graph_diameter', 
        type=int,
        default=6, 
        help="(default: 6) The graph diameter."
    )
    
    parser.add_argument(
        '--group_cutoff', 
        type=float,
        default=3.0, 
        help="(default: 3.0) Scalar for community size cutoff."
    )

    parser.add_argument(
        '--lowest_freq', 
        type=float,
        default=0.1, 
        help="(default: 0.10) Lower bound for allele frequency."
    )
    
    # Parse the arguments
    args = parser.parse_args()
    
    # Adjust values as needed
    K = int(args.graph_diameter / 2)
    Q1 = args.lowest_freq
    Q2 = 1 - Q1
    scalar = args.group_cutoff

    # Form graph
    segs = read_ibd_file(args.ibd_file, header=0, include_length=0)
    graph = make_ibd_graph(segs)

    # Detect communities
    communities = diameter_communities(graph, K=K, max_communities=np.inf)
    outliers = outlier_communities(communities, scalar=scalar)

    # Compute adaptive allele frequencies
    tup = labeled_allele_frequencies(args.vcf, outliers)

    # Make table
    pos = tup[0]
    freq1, freq0, freqm = putative_allele_frequencies(tup[1], tup[2], tup[3])
    table = format_allele_table(pos, freq1, freq0, freqm)
    table['ZDELTA'] = table['DELTA'] / np.sqrt(table['AAF'] * (1 - table['AAF']))
    table.sort_values(['ZDELTA'], inplace=True, ascending=False)
    table = table[table['AAF'] >= Q1]
    table = table[table['AAF'] <= Q2]
    table.reset_index(inplace=True, drop=True)
    table.to_csv(args.file_out, sep='\t', index=False)

if __name__ == "__main__":
    main()
