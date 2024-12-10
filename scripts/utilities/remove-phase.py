import gzip
import sys
import argparse

parser = argparse.ArgumentParser(description="Unphase a gzipped phased VCF file.")

parser.add_argument('--phased_vcf_file', 
                    type=str,
                    required=True, 
                    help='Input gzipped VCF file with phased genotypes')
parser.add_argument('--unphased_vcf_file', 
                    type=str,
                    required=True, 
                    help='Output gzipped VCF file with unphased genotypes')

args = parser.parse_args()

g = gzip.open(args.unphased_vcf_file,'wt')

header = True

with gzip.open(args.phased_vcf_file,'rt') as f:

    # process the header
    while header:
        line = f.readline()
        if line[0] == '#':
            g.write(line)
        else:
            header = False

    # do unphase first line not header
    full_row = line.split()
    info_data = full_row[:9]
    geno_data = full_row[9:]
    geno_data = [geno.replace('|','/') for geno in geno_data]
    g.write('\t'.join(info_data))
    g.write('\t')
    g.write('\t'.join(geno_data))
    g.write('\n')

    for line in f:
        full_row = line.split()
        info_data = full_row[:9]
        geno_data = full_row[9:]
        geno_data = [geno.replace('|','/') for geno in geno_data]
        g.write('\t'.join(info_data))
        g.write('\t')
        g.write('\t'.join(geno_data))
        g.write('\n')
            
g.close()