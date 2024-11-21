import gzip
import sys
import argparse

parser = argparse.ArgumentParser(description="Unphase a gzipped phased VCF file.")

parser.add_argument('phased_file', 
                    type=str, 
                    help='Input gzipped VCF file with phased genotypes')
parser.add_argument('unphased_file', 
                    type=str, 
                    help='Output gzipped VCF file with unphased genotypes')

args = parser.parse_args()

g = gzip.open(args.unphased_file,'wb')

header = True

with gzip.open(args.phased_file,'rb') as f:

    # process the header
    while header:
        line = f.readline()
        if line[0] == b'#':
            g.write(line)
        else:
            header = False

    # do unphase first line not header
    full_row = line.split()
    info_data = full_row[:9]
    geno_data = full_row[9:]
    geno_data = [geno.replace(b'|',b'/') for geno in geno_data]
    g.write(b'\t'.join(info_data))
    g.write(b'\t')
    g.write(b'\t'.join(geno_data))
    g.write(b'\n')

    for line in f:
        full_row = line.split()
        info_data = full_row[:9]
        geno_data = full_row[9:]
        geno_data = [geno.replace(b'|',b'/') for geno in geno_data]
        g.write(b'\t'.join(info_data))
        g.write(b'\t')
        g.write(b'\t'.join(geno_data))
        g.write(b'\n')
            
g.close()