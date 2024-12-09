import pandas as pd
import argparse

# Set up argument parser
parser = argparse.ArgumentParser(description='Trim telomeres from genetic map.')
parser.add_argument('--input_map_file', 
                    type=str,
                    required=True, 
                    help='Input genetic map')
parser.add_argument('--output_map_file', 
                    type=str,
                    required=True, 
                    help='Output genetic map')
parser.add_argument('--cM_telo_trim', 
                    type=float,
                    default=2., 
                    help='Cut this many cM from map ends')

# Parse the arguments
args = parser.parse_args()
map_in = args.input_map_file
map_out = args.output_map_file
telocut = args.cM_telo_trim

# Load the table
table = pd.read_csv(map_in, sep='\t', header=None)
cmcol = table.columns.to_list()[2]
cut = float(telocut)
mn = table[cmcol].min() + cut
mx = table[cmcol].max() - cut
subtable = table[(table[cmcol] >= mn) & (table[cmcol] <= mx)]
subtable.to_csv(map_out, sep='\t', index=False, header=False)