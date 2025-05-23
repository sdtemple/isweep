# Find haplotype groups of excess IBD.
# Seth Temple, GitHub: sdtemple
# Converted to argparse on 2024-07-18

import argparse
import numpy as np
from copy import deepcopy
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
        communities2 = deepcopy(communities)
        communities2 = sorted(communities2, key=len, reverse=True) 
        communities2 = [list(c) for c in communities2]
        community_sizes = [len(community) for community in communities]
        community_sizes = np.array(community_sizes)
        try:
            cutoff = community_sizes.mean() + community_sizes.std() * scalar
            idx = 1
            for community in communities:
                if len(community) >= cutoff:
                    with open(f"{folderout}/outlier{idx}.phenotype.tsv", 'w') as f:
                        for haplo in community:
                            f.write(haplo)
                            f.write('\t')
                            f.write(str(phen[haplo[:-2]]))
                            f.write('\n')
                    idx += 1
            if idx == 1:
                with open(f"{folderout}/outlier{idx}.phenotype.tsv", 'w') as f:
                    for haplo in communities2[0]:
                        f.write(haplo)
                        f.write('\t')
                        f.write(str(phen[haplo[:-2]]))
                        f.write('\n')
            with open(f"{folderout}/communities.phenotype.csv",'w') as f:
                for community in communities2:
                    for haplo in community[:-1]:
                        f.write(haplo); f.write(':')
                        f.write(str(phen[haplo[:-2]])) ;f.write(',')
                    f.write(haplo); f.write(':')
                    f.write(str(phen[haplo[:-2]])); f.write('\n')
        except IndexError:
            # scenario in which there is no IBD
            print(len(communities))
            f=open(f"{folderout}/outlier1.phenotype.tsv", 'w')
            f.close()

    # Saving
    write_outliers(communities, phen, args.output_folder, scalar)

if __name__ == "__main__":
    main()
