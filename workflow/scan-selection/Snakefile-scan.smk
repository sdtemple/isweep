# The main organizing file for the selection scan,
# where the rule all is and many of the yaml parameters.
# You run this file with the -s parameter in snakemake.

import shutil
import os

# setup macro folder
macro=str(config['change']['files']['study'])
chrlow=int(str(config['change']['files']['chromosome_low']))
chrhigh=int(str(config['change']['files']['chromosome_high']))
low=chrlow
high=chrhigh

# vcf data
vcffolder=str(config['change']['files']['vcfs'])
vcfpre=str(config['change']['files']['vcf_prefix'])
vcfsuf=str(config['change']['files']['vcf_suffix'])

# calculate the sample size
subsample=str(config['change']['files']['subsample'])
samplesize=0
with open(subsample,'r') as f:
    for line in f:
        samplesize+=1

ploidy = int(float(config['change']['files']['ploidy']))

# number of sims for the simulation method
numsims = int(str(config['change']['isweep']['num_sims'])) 

# copy and paste maps
folder_name = macro
if not os.path.exists(folder_name):
    os.makedirs(folder_name)
folder_name = macro + '/maps'
if not os.path.exists(folder_name):
    os.makedirs(folder_name)
folder_name = macro + '/ibdsegs'
if not os.path.exists(folder_name):
    os.makedirs(folder_name)
folder_name = macro + '/ibdsegs/hapibd'
if not os.path.exists(folder_name):
    os.makedirs(folder_name)
folder_name = macro + '/ibdsegs/ibdends'
if not os.path.exists(folder_name):
    os.makedirs(folder_name)
folder_name = macro + '/ibdsegs/ibdends/mle'
if not os.path.exists(folder_name):
    os.makedirs(folder_name)
folder_name = macro + '/ibdsegs/ibdends/scan'
if not os.path.exists(folder_name):
    os.makedirs(folder_name)

mapfol=str(config['change']['files']['genetic_maps'])
mappre=str(config['change']['files']['map_prefix'])
mapsuf=str(config['change']['files']['map_suffix'])

mapexcl=str(config['change']['files']['chromosome_exclude'])
exclude=[]
if not os.path.exists(mapexcl):
    pass
else:
    with open(mapexcl,'r') as f:
        for line in f:
            exclude.append(int(float(line.strip())))
chroms = [i for i in range(chrlow,chrhigh+1) if i not in exclude]

import pandas as pd
sizes_file = macro+'/chromosome-sizes.tsv'
sizes_out = open(sizes_file,'w')
sizes_out.write('CHROM\tCMSIZE\n')
for i in range(chrlow,chrhigh+1):
    source_file = mapfol+'/'+mappre+str(i)+mapsuf
    sizes_table = pd.read_csv(source_file,sep='\t',header=None)
    start_cm = float(sizes_table[2].tolist()[0])
    end_cm = float(sizes_table[2].tolist()[-1])
    size_cm = end_cm - start_cm
    sizes_out.write(str(i)); sizes_out.write('\t')
    sizes_out.write(str(size_cm)); sizes_out.write('\n')
sizes_out.close()

sizes_kept = macro+'/chromosome-sizes-kept.tsv'
sizes_kept_out = open(sizes_kept,'w')
sizes_kept_out.write('CHROM\tCMSIZE\n')
min_chr_size = float(str(config['fixed']['isweep']['chromosome_size_cutoff']))
exclude2 = []
total_size = 0
total_chr = 0
for i in chroms:
    source_file = mapfol+'/'+mappre+str(i)+mapsuf
    sizes_table = pd.read_csv(source_file,sep='\t',header=None)
    start_cm = float(sizes_table[2].tolist()[0])
    end_cm = float(sizes_table[2].tolist()[-1])
    size_cm = end_cm - start_cm
    if size_cm >= min_chr_size:
        total_size += size_cm
        total_chr += 1
        sizes_kept_out.write(str(i)); sizes_kept_out.write('\t')
        sizes_kept_out.write(str(size_cm)); sizes_kept_out.write('\n')
    else:
        exclude2.append(i)
sizes_kept_out.close()
chroms2 = [i for i in chroms if i not in exclude2]

for i in chroms2:
    source_file = mapfol+'/'+mappre+str(i)+mapsuf
    destination_file = macro+'/maps/chr'+str(i)+'.map'
    if not os.path.exists(destination_file):
        shutil.copy(source_file, destination_file)

# include .smk files with rules
include: 'rules/scan.smk'
include: 'rules/interesting.smk'
include: 'rules/samples.smk'
include: 'rules/maps.smk'
include: 'rules/fwer.smk'

# record the configuration settings used
rule yaml:
    input:
        macro+'/excess.ibd.tsv',
        macro+'/excess.region.ibd.tsv',
        macro+'/roi.tsv',
        macro+'/excludesamples.txt',
        macro+'/zhistogram.png',
        macro+'/autocovariance.png',
        macro+'/scan.png'
    output:
        yaml=macro+'/arguments.scan.yaml',
    params:
        yaml=str(config['change']['files']['yaml']),
    shell:
        'cp {params.yaml} {output.yaml}'

# snakemake all -c1 -n
rule all:
    input:
        macro+'/arguments.scan.yaml',
