import argparse
import pandas as pd

# Set up the argument parser
parser = argparse.ArgumentParser(description='Find the best scored SNP position and frequency.')

# Define required arguments with named options
parser.add_argument('--input_snp_file', 
                    type=str, 
                    required=True, 
                    help='Input table file')
parser.add_argument('--output_file', 
                    type=str, 
                    required=True, 
                    help='Output file')
parser.add_argument('--lowest_freq', 
                    type=float, 
                    default=0.1, 
                    help='(default: 0.1) Lower frequency bound of AAF to filter rows')

# Parse the arguments
args = parser.parse_args()
snp_file_input = args.input_snp_file
file_output = args.output_file
low_freq = args.lowest_freq

# Read the input table
table = pd.read_csv(snp_file_input, sep='\t')

# Filter the table rows where 'AAF' is greater than or equal to 'low_freq'
table = table[table['AAF'] >= low_freq]

# Extract the best position and frequency from the filtered table
bestbp = table['POS'].tolist()[0]
bestaf = table['AAF'].tolist()[0]

# Write the results to the output file
with open(file_output, 'w') as f:
    f.write('bp\t')
    f.write(str(int(bestbp)))
    f.write('\n')
    f.write('frequency\t')
    f.write(str(round(bestaf, 4)))
    f.write('\n')