import argparse

# Set up the argument parser
parser = argparse.ArgumentParser(description='Exclude samples not in sub-sample list.')
parser.add_argument('--sample_file', 
                    type=str,
                    required=True, 
                    help='Path to the sample file')
parser.add_argument('--subsample_file', 
                    type=str,
                    required=True, 
                    help='Path to the sub-sample file')
parser.add_argument('--exclude_sample_file', 
                    type=str,
                    required=True, 
                    help='Path to the output file for excluded samples')

# Parse the arguments
args = parser.parse_args()
sampfile = args.sample_file
subsampfile = args.subsample_file
exclsampfile = args.exclude_sample_file

# Read sample list
samplist = []
with open(sampfile, 'r') as f:
    for line in f:
        samplist.append(line.strip())

# Read sub-sample list
subsamplist = []
with open(subsampfile, 'r') as f:
    for line in f:
        subsamplist.append(line.strip())

# set minus operation
sampset = set(samplist)
subsampset = set(subsamplist)
exclsampset = sampset - subsampset
exclsamplist = list(exclsampset)

# Write excluded samples to output file
with open(exclsampfile, 'w') as f:
    for exclsamp in exclsamplist:
        f.write(str(exclsamp) + '\n')