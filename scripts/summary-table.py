## File to summarize all of the region-specific analyses
## July 8, 2024

import argparse
import pandas as pd

# Set up the argument parser
parser = argparse.ArgumentParser(description='Summarize all region-specific analyses.')

# Define required arguments with named options
parser.add_argument('--output_file', 
                    type=str, 
                    required=True, 
                    help='Path to the output file')
parser.add_argument('--input_folder', 
                    type=str, 
                    required=True, 
                    help='Path to the folder containing region data')
parser.add_argument('--input_roi_file', 
                    type=str, 
                    required=True, 
                    help='Region of interest file name')
parser.add_argument('--file_type', 
                    type=str, 
                    required=True, 
                    help='Haplotype or SNP-based analysis')
parser.add_argument('--uncertainty_type', 
                    type=str, 
                    required=True, 
                    help='Normal or percentile bootstrap intervals')

# Parse the arguments
args = parser.parse_args()
file_output = args.output_file
folder = args.input_folder
roi = args.input_roi_file
file_type = args.file_type
uncertainty_type = args.uncertainty_type

assert file_type in ['hap','snp']
assert uncertainty_type in ['norm','perc']

# Read the region of interest file
tablein = pd.read_csv(f"{folder}/{roi}", sep='\t')
tablein = tablein[['NAME', 'CHROM', 'MAXIBD', 'ALPHA']]
nms = list(tablein['NAME'])

p0 = []
se = []
sl = []
su = []
ms = []
bs = []

# Process each analysis file
for nm in nms:
    with open(f"{folder}/{nm}/third.best.{file_type}.txt", 'r') as f:
        line = f.readline().strip().split('\t')
        bp = int(float(line[1]))
    bs.append(bp)

    restab = pd.read_csv(f"{folder}/{nm}/results.{file_type}.{uncertainty_type}.tsv", sep='\t')
    p0.append(restab['VarFreqEst'][0])
    se.append(restab['SelCoefEst'][0])
    sl.append(restab['SelCoefLow'][0])
    su.append(restab['SelCoefUpp'][0])
    ms.append(restab['Model'][0])

# Add the results to the table
tablein['LOCHAT'] = bs
tablein['PHAT'] = p0
tablein['SHAT'] = se
tablein['CONF_INT_LOW'] = sl
tablein['CONF_INT_UPP'] = su
tablein['MODEL'] = ms

# Write the summarized table
tablein.to_csv(file_output, sep='\t', index=False, header=True)