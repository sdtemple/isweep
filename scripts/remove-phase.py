import gzip
import sys
import argparse

parser = argparse.ArgumentParser(description="Unphase a phased VCF file.")

parser.add_argument('phased_file', 
                    type=str, 
                    help='Input VCF file with phased genotypes')
parser.add_argument('unphased_file', 
                    type=str, 
                    help='Output VCF file with unphased genotypes')

args = parser.parse_args()

g = gzip.open(args.unphased_file,'wb')

with gzip.open(args.phased_file,'rb') as f:
    for line in f:
        if line[:2] == b'##':
            g.write(line)
        elif line[:2] == b'#C':
            g.write(line)
        else:
            full_row = line.split()
            info_data = full_row[:9]
            geno_data = full_row[9:]
            geno_data = [geno.replace(b'|',b'/') for geno in geno_data]
            g.write(b'\t'.join(info_data))
            g.write(b'\t')
            g.write(b'\t'.join(geno_data))
            g.write(b'\n')
            
g.close()