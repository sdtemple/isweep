from isweep import *
import networkx as nx
from math import floor
import argparse

def main():

    # Set up argparse
    parser = argparse.ArgumentParser(
        description="Detect IBD graph communities"
    )
    
    # Define arguments

    parser.add_argument(
        '--input_ibd_file',
        type=str,
        default='second.filt.ibd.gz',
        help="File name of IBD segments at a locus"
    )

    parser.add_argument(
        '--output_connected_file', 
        type=str,
        default='connected.communities.csv', 
        help="File name of connected components"
    )
    
    parser.add_argument(
        '--output_diameter_file', 
        type=str,
        default='diameter.communities.csv', 
        help="File name of diameter communities"
    )
    
    parser.add_argument(
        '--diameter', 
        type=int,
        default=6, 
        help="(default: 6) Graph diameter."
    )

    # Parse the arguments
    args = parser.parse_args()

    # Implementing the code for outlier detection in rank.py
    segs = read_ibd_file(args.input_ibd_file, header=0, include_length=0)
    connected_file = args.output_connected_file
    diameter_file = args.output_diameter_file

    if len(segs) > 0:

        graph = make_ibd_graph(segs)
        connected = [list(community) for community in nx.connected_components(graph)]
        connected = sorted(connected,key=len,reverse=True)

        K = floor(args.diameter / 2)
        dia = sorted(diameter_communities(graph,K=3),key=len,reverse=True)
        dia = [list(community) for community in dia]

        with open(connected_file,'w') as f:
            for community in connected:
                for c in community[:-1]:
                    f.write(str(c)); f.write(str(','))
                f.write(str(community[-1])); f.write('\n')

        with open(diameter_file,'w') as f:
            for community in dia:
                for c in community[:-1]:
                    f.write(str(c)); f.write(str(','))
                f.write(str(community[-1])); f.write('\n')

    else:
        # when there is no IBD
        f = open(connected_file,'w'); f.close()
        g = open(diameter_file,'w'); g.close()

    return None

if __name__ == "__main__":
    main()