# Find haplotype groups of excess IBD.
# Seth Temple, GitHub: sdtemple
# Converted to argparse on 2024-07-18

import argparse
import numpy as np
from isweep import *

def main():
    # Set up argparse
    parser = argparse.ArgumentParser(
        description="Find haplotype groups of excess IBD."
    )
    
    # Define arguments with defaults
    parser.add_argument(
        '--input_ibd_file', 
        type=str,
        required=True, 
        help="Path to the input IBD file."
    )

    # Define arguments with defaults
    parser.add_argument(
        '--input_phenotype_file', 
        type=str,
        required=True, 
        help="Path to the input phenotype file."
    )
    
    parser.add_argument(
        '--output_folder', 
        type=str,
        required=True, 
        help="Output folder for writing results."
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
    
    # Parse the arguments
    args = parser.parse_args()

    phen=dict()
    with open(args.input_phenotype_file,'r') as f:
        for line in f:
            hap,sta=line.strip().split('\t')
            phen[hap]=float(sta)
    
    # Adjust values as needed
    K = int(args.graph_diameter / 2)
    scalar = args.group_cutoff

    # Form graph
    segs = read_ibd_file(args.input_ibd_file, header=0, include_length=0)
    graph = make_ibd_graph(segs)

    # Detect communities
    communities = diameter_communities(graph, K=K, max_communities=np.inf)

    # Write clusters to file
    def write_outliers(communities, phen, folderout, scalar=3):
        '''Write haplotypes in outlier big communities.

        Parameters
        ----------
        communities : list of sets
            List of communities, each represented as a set of haplotypes.
        phen: dict
            Mapping between individual and categorical phenotype
        folderout : str
            Directory to write outlier files.
        scalar : float
            Multiplier of numpy .std().

        Returns
        -------
        None
        '''
        community_sizes = [len(community) for community in communities]
        community_sizes = np.array(community_sizes)
        cutoff = community_sizes.mean() + community_sizes.std() * scalar
        idx = 1
        for community in communities:
            if len(community) > cutoff:
                with open(f"{folderout}/outlier{idx}.phenotype.tsv", 'w') as f:
                    for haplo in community:
                        f.write(haplo)
                        f.write('\t')
                        f.write(str(phen[haplo[:-2]]))
                        f.write('\n')
                idx += 1

    # Saving
    write_outliers(communities, phen, args.output_folder, scalar)

if __name__ == "__main__":
    main()
