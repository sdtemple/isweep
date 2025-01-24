import argparse
import gzip

parser = argparse.ArgumentParser(description='Convert a VCF of any ploidy to a haploid VCF.')
parser.add_argument('--input_vcf_file', 
                    type=str, 
                    help='Path to the input VCF file (gzipped).')
parser.add_argument('--output_vcf_file', 
                    type=str, 
                    help='Path to the output VCF file (gzipped).')
parser.add_argument('--input_samples_file', 
                    type=str, 
                    help='Path to file with sample IDs.')
parser.add_argument('--ploidy', 
                    type=int, 
                    help='Ploidy of the input VCF.')

args = parser.parse_args()
vcf_input = args.input_vcf_file
vcf_output = args.output_vcf_file
sample_text_file = args.input_samples_file
ploidy = args.ploidy

samples = []
with open(sample_text_file,'r') as f:
    for line in f:
        identity = line.strip()
        for p in range(1,ploidy+1):
            samples.append(identity+'_'+str(p))
samples_header = '\t'.join(samples)
non_samples_header = '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT'
header = non_samples_header + '\t' + samples_header

g = gzip.open(vcf_output,'wt')
with gzip.open(vcf_input,'rt') as f:
    for line in f:
        if line[:2] == '##':
            g.write(line.strip())
            g.write('\n')
        elif line[:2] == '#C':
            g.write(header)
            g.write('\n')
        else:
            columns = line.strip().split('\t')
            metadata = columns[:9]
            g.write('\t'.join(metadata))
            nonmetadata = columns[9:]
            g.write('\t')
            for individual in nonmetadata[:-1]:
                alleles = individual.split('|')
                alleles = [allele[0] for allele in alleles]
                g.write('\t'.join(alleles))
                g.write('\t')
            individual = nonmetadata[-1]
            alleles = individual.split('|')
            alleles = [allele[0] for allele in alleles]
            g.write('\t'.join(alleles))
            g.write('\n')
g.close()