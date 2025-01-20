##### An analysis pipeline for flare, hap-ibd, and beagle

import os

### setup macro folder
macro=str(config['change']['want-data']['your-analysis-folder'])
low=int(str(config['change']['want-data']['chr-low']))
high=int(str(config['change']['want-data']['chr-high']))

### setup the folder structure

# make the main folder
folder_name = macro
if not os.path.exists(folder_name):
    os.makedirs(folder_name)

# make the genotypes data folder
folder_name = macro + '/gtdata'
if not os.path.exists(folder_name):
    os.makedirs(folder_name)

# make the genotypes data for references
folder_name = macro + '/gtdata/refpop'
if not os.path.exists(folder_name):
    os.makedirs(folder_name)

# make the genotypes data for admixed samples to study
folder_name = macro + '/gtdata/adxpop'
if not os.path.exists(folder_name):
    os.makedirs(folder_name)

# make the genotypes folder for refs + admixeds
folder_name = macro + '/gtdata/all'
if not os.path.exists(folder_name):
    os.makedirs(folder_name)

# make the ibd segments subfolder
folder_name = macro + '/ibdsegs'
if not os.path.exists(folder_name):
    os.makedirs(folder_name)

# make the lai subfolder
folder_name = macro + '/lai'
if not os.path.exists(folder_name):
    os.makedirs(folder_name)

# make the maps subfolder
folder_name = macro + '/maps'
if not os.path.exists(folder_name):
    os.makedirs(folder_name)

# prime the chromosome wildcards
import shutil
mapfol=str(config['change']['existing-data']['maps-folder'])
mappre=str(config['change']['existing-data']['maps-prefix'])
mapsuf=str(config['change']['existing-data']['maps-suffix'])
for i in range(low,high+1):
    source_file = mapfol+'/'+mappre+str(i)+mapsuf
    destination_file = macro+'/maps/chr'+str(i)+'.map'
    if not os.path.exists(destination_file):
        shutil.copy(source_file, destination_file)

### what the pipeline demands

# include .smk files with rules
include: 'rules/ibd.smk'
include: 'rules/lai.smk'
include: 'rules/prepping.smk'
include: 'rules/phasing.smk'

# do the dry run `snakemake all -c1 -n`
# this tells you what the pipeline will do
# this rule will also copy the arguments you used

rule record_yaml:
    input:
        [macro+'/lai/chr'+str(i)+'.referencephased.flare.anc.vcf.gz' for i in range(low,high+1)],
        # [macro+'/lai/chr'+str(i)+'.rephased.flare.anc.vcf.gz' for i in range(low,high+1)],
        # local ancestry inference
        [macro+'/gtdata/adxpop/chr'+str(i)+'.referencephased.vcf.gz' for i in range(low,high+1)],
        # [macro+'/gtdata/adxpop/chr'+str(i)+'.rephased.vcf.gz' for i in range(low,high+1)],
        # [macro+'/gtdata/refpop/chr'+str(i)+'.rephased.vcf.gz' for i in range(low,high+1)],
        # phasing
        # [macro+'/ibdsegs/chr'+str(i)+'.referencephased.ref.hapibd.ibd.gz' for i in range(low,high+1)],
        [macro+'/ibdsegs/chr'+str(i)+'.referencephased.adx.hapibd.ibd.gz' for i in range(low,high+1)],
        # [macro+'/ibdsegs/chr'+str(i)+'.rephased.ref.hapibd.ibd.gz' for i in range(low,high+1)],
        # [macro+'/ibdsegs/chr'+str(i)+'.rephased.adx.hapibd.ibd.gz' for i in range(low,high+1)],
        # ibd segment detection
    output:
        yaml=macro+'/arguments.yaml',
    params:
        yaml=str(config['change']['pipe']['yaml']),
    shell:
        'cp {params.yaml} {output.yaml}'

rule all:
    input:
        macro+'/arguments.yaml'
