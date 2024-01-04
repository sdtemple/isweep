# isweep real data analysis part 2 (roi)
# analyses of regions of interest

import os
import pandas as pd
macro=str(config['CHANGE']['FOLDERS']['STUDY'])
micro=str(config['CHANGE']["ISWEEP"]["ROI"])
sims = pd.read_csv(macro+'/'+micro, sep='\t', header=0)
J = sims.shape[0]
for j in range(J):
	row = sims.iloc[j]
	if not os.path.exists(macro+'/'+str(row.NAME)):
		os.mkdir(macro+'/'+str(row.NAME))
	f=open(macro+'/'+str(row.NAME)+'/locus.txt','w')
	f.write('NAME\t')
	f.write(str(row.NAME))
	f.write('\n')
	f.write('CHROM\t')
	f.write(str(int(row.CHROM)))
	f.write('\n')
	f.write('CENTER\t')
	f.write(str(int(row.BPCENTER)))
	f.write('\n')
	f.write('LEFT\t')
	f.write(str(int(row.BPLEFTCENTER)))
	f.write('\n')
	f.write('RIGHT\t')
	f.write(str(int(row.BPRIGHTCENTER)))
	f.write('\n')
	f.close()
sims['FOLDER'] = [(macro +'/'+str(sims.iloc[j].NAME)).strip() for j in range(J)]

# include .smk files with rules
include: 'rules/first.smk'
include: 'rules/second.smk'
include: 'rules/third.smk'

# snakemake all -c1 -n
rule all:
	input:
		# for best snp analysis
		[(macro +'/'+str(sims.iloc[j].NAME)).strip()+'/results.snp.tsv' for j in range(J)],
		# for best haplotype analysis
		[(macro +'/'+str(sims.iloc[j].NAME)).strip()+'/results.hap.tsv' for j in range(J)],
		# for entropy
		[(macro +'/'+str(sims.iloc[j].NAME)).strip()+'/ibd.entropy.tsv' for j in range(J)],

