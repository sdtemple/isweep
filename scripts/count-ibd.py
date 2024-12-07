# Computing IBD counts over windows

# importing
import argparse
import numpy as np
import pandas as pd

# Set up the argument parser
parser = argparse.ArgumentParser(description='Compute IBD counts over windows.')
parser.add_argument('input_file', 
                    type=str, 
                    help='Input file path')
parser.add_argument('map_file', 
                    type=str, 
                    help='Map file path')
parser.add_argument('output_file', 
                    type=str, 
                    help='Output file path')
parser.add_argument('--start', 
                    type=int,
                    default=5, 
                    help='Start bp column index')
parser.add_argument('--end', 
                    type=int,
                    default=6, 
                    help='End bp column index')

# Parse the arguments
args = parser.parse_args()
input_file = args.input_file
map_file = args.map_file
output_file = args.output_file
start = args.start
end = args.end

# formatting
table = pd.read_csv(input_file, sep='\t')
map_file = pd.read_csv(map_file, sep='\t')
columns = list(table.columns)
encol = columns[end]
stcol = columns[start]
map_file.columns = ['chrom', 'rsid', 'cm', 'bp']

# counting IBD segments
counts = [((table[stcol] <= i) & (table[encol] >= i)).sum() for i in map_file['bp']]
counter = {
    'BPWINDOW': map_file['bp'].to_list(),
    'CMWINDOW': map_file['cm'].to_list(),
    'COUNT': counts
}
counter = pd.DataFrame(counter)

# count cutting
counter = counter[counter['COUNT'] > 0]

# saving
counter.to_csv(output_file, sep='\t', index=False)