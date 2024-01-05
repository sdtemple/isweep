# scanning for excess ibd

# setup macro folder
macro=str(config['CHANGE']['FOLDERS']['STUDY'])
low=int(str(config['CHANGE']['ISWEEP']['CHRLOW']))
high=int(str(config['CHANGE']['ISWEEP']['CHRHIGH']))

# copy and paste maps
import os
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
# if not os.path.exists(folder_name):
#     os.makedirs(folder_name)
# folder_name = macro + '/ibdsegs/ibdends/modified'
if not os.path.exists(folder_name):
    os.makedirs(folder_name)
# folder_name = macro + '/ibdsegs/ibdends/modified/mle'
folder_name = macro + '/ibdsegs/ibdends/mle'
if not os.path.exists(folder_name):
    os.makedirs(folder_name)
# folder_name = macro + '/ibdsegs/ibdends/modified/scan'
folder_name = macro + '/ibdsegs/ibdends/scan'
if not os.path.exists(folder_name):
    os.makedirs(folder_name)

mapfol=str(config['CHANGE']['EXISTING']['MAPS'])
mappre=str(config['CHANGE']['EXISTING']['MAPPRE'])
mapsuf=str(config['CHANGE']['EXISTING']['MAPSUF'])
import shutil
source_file = "path/to/source/file"
destination_file = "path/to/destination/file"
for i in range(low,high+1):
    source_file = mapfol+'/'+mappre+str(i)+mapsuf
    destination_file = macro+'/maps/chr'+str(i)+'.map'
    if not os.path.exists(destination_file):
        shutil.copy(source_file, destination_file)
    # shutil.copy(source_file, destination_file)

# include .smk files with rules
include: 'rules/scan.smk'
include: 'rules/interesting.smk'
include: 'rules/samples.smk'

rule yaml:
    input:
        macro+'/excess.ibd.tsv',
        macro+'/excess.region.ibd.tsv',
        macro+'/roi.tsv',
        macro+'/excludesamples.txt',
    output:
        yaml=macro+'/arguments.scan.yaml',
    params:
        yaml=str(config['CHANGE']['FOLDERS']['YAML']),
    shell:
        'cp {params.yaml} {output.yaml}'

# snakemake all -c1 -n
rule all:
    input:
        macro+'/arguments.scan.yaml',
