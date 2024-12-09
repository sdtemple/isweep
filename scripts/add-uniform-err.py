import gzip
import sys
import argparse
import numpy as np

parser = argparse.ArgumentParser(description="Add genotype errors by switching alleles. No dosages. Biallelic only.")

parser.add_argument('--input_file', 
                    type=str,
                    required=True, 
                    help='Input gzipped VCF file')
parser.add_argument('--output_file', 
                    type=str,
                    required=True, 
                    help='Output gzipped VCF file with genotyped errors')
parser.add_argument('--unif_err_rate',
                    type=float,
                    default=0.,
                    help='Rate at with genotype error happens'
                    )
parser.add_argument('--phase_status',
                    type=int,
                    default=1,
                    help='If 0, it is not phased'
                    )

args = parser.parse_args()
err_rate = args.unif_err_rate
if args.phase_status:
    splitter = '|'
else:
    splitter = '/'

g = gzip.open(args.output_file,'wt')

header = True

with gzip.open(args.input_file,'rt') as f:

    # process the header
    while header:
        line = f.readline()
        if line[0] == '#':
            g.write(line)
        else:
            header = False

    # do introduce errors to first line not header
    full_row = line.split()
    info_data = full_row[:9]
    geno_data = full_row[9:]
    geno_data = [splitter.join([
        str(( int(j) + (np.random.uniform() <= err_rate) ) % 2)
        for j in geno.split(splitter)]) 
                    for geno in geno_data]
    g.write('\t'.join(info_data))
    g.write('\t')
    g.write('\t'.join(geno_data))
    g.write('\n')

    for line in f:
        full_row = line.split()
        info_data = full_row[:9]
        geno_data = full_row[9:]
        geno_data = [splitter.join([
            str(( int(j) + (np.random.uniform() <= err_rate) ) % 2)
            for j in geno.split(splitter)]) 
                        for geno in geno_data]
        g.write('\t'.join(info_data))
        g.write('\t')
        g.write('\t'.join(geno_data))
        g.write('\n')
            
g.close()