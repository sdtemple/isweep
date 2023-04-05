### setup ###

import gzip
import pandas as pd

vcf_file=snakemake.input.vcf
labels1=snakemake.output.labels1
labels0=snakemake.output.labels0

def get_vcf_names(vcf_path):
    with gzip.open(vcf_path, "rt") as ifile:
        for line in ifile:
            if line.startswith("#CHROM"):
                vcf_names = [x for x in line.split('\t')]
                break
    ifile.close()
    vcf_names[-1] = vcf_names[-1].strip()
    return vcf_names

### true labels for haploids ###

names = get_vcf_names(vcf_file)
vcf = pd.read_csv(vcf_file,
                  compression='gzip',
                  comment='#',
                  delim_whitespace=True,
                  header=None,
                  names=names)

vcf.drop(['#CHROM','POS','ID','REF','ALT','QUAL','FILTER','INFO','FORMAT'], inplace = True, axis=1)
vcf = vcf.transpose()

hap1 = vcf[0].apply(lambda x: x[0])
hap2 = vcf[0].apply(lambda x: x[2])

ids = list(hap1.index)

ct = 0
rows = {}
for i in range(len(ids)):
    row = []
    row.append(ids[i] + '_1')
    row.append(int(hap1[i]))
    rows[ct] = row
    ct += 1
    row = []
    row.append(ids[i] + '_2')
    row.append(int(hap2[i]))
    rows[ct] = row
    ct += 1

lbl = pd.DataFrame(rows)
lbl = lbl.transpose()
lbl.columns = ['ID', 'LABEL']

lbl0 = lbl[lbl['LABEL']==0]
lbl1 = lbl[lbl['LABEL']==1]
lbl0.to_csv(labels0, header=True, index=False, compression='gzip', sep='\t')
lbl1.to_csv(labels1, header=True, index=False, compression='gzip', sep='\t')
