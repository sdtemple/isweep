# The main organizing file for the sweep modeling,
# where the rule all is and many of the yaml parameters.
# You run this file with the -s parameter in snakemake.

import os
import pandas as pd
macro=str(config['change']['files']['study'])
cohort=macro
micro=str(config['change']["files"]["regions_of_interest"])
roi=micro
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
	f.write('MODEL\t')
	f.write(str(row.MODEL))
	f.write('\n')
	f.write('ALPHA\t')
	f.write(str(float(row.ALPHA)))
	f.close()
sims['FOLDER'] = [(macro +'/'+str(sims.iloc[j].NAME)).strip() for j in range(J)]

samplesize=0
with open(cohort+'/subsample.txt','r') as f:
    for line in f:
        samplesize+=1
ploidy=int(float(config['change']['files']['ploidy']))
samplesize_ploidy=samplesize*ploidy

maf=str(config['fixed']['isweep']['min_adaptive_allele_frequency'])
diameter=str(config['fixed']['isweep']['diameter'])
group_cutoff=str(config['fixed']['isweep']['group_cutoff'])

# include .smk files with rules
include: 'rules/first.smk'
include: 'rules/second.smk'
include: 'rules/third.smk'

# snakemake all -c1 -n
rule all:
	input:
		# for best snp analysis
		#[(macro +'/'+str(sims.iloc[j].NAME)).strip()+'/results.snp.tsv' for j in range(J)],
		# for best haplotype analysis
		#[(macro +'/'+str(sims.iloc[j].NAME)).strip()+'/results.hap.tsv' for j in range(J)],
		# for entropy
		# [(macro +'/'+str(sims.iloc[j].NAME)).strip()+'/ibd.gini.tsv' for j in range(J)],
		# for types of summary analyses
		[macro+'/summary.hap.norm.tsv'], # best haploptye, normal intervals
		# [macro+'/summary.hap.perc.tsv'], # best haplotype, percentile intervals
		[macro+'/summary.snp.norm.tsv'], # best snp, normal intervals
		# [macro+'/summary.snp.perc.tsv'], # best snp, percentile intervals
	output:
		yaml=macro+'/arguments.roi.yaml',
	params:
		yaml=str(config['change']['files']['yaml']),
	shell:
		'cp {params.yaml} {output.yaml}'

