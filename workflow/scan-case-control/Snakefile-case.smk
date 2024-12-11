# scanning for excess ibd

import shutil

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
    raise Exception("Files don't exist. Run workflow/scan-selection beforehand.")
folder_name = macro + '/ibdsegs/ibdends'
if not os.path.exists(folder_name):
    raise Exception("Files don't exist. Run workflow/scan-selection beforehand.")
folder_name = macro + '/ibdsegs/ibdends/scan'
if not os.path.exists(folder_name):
    raise Exception("Files don't exist. Run workflow/scan-selection beforehand.")

include: 'rules/interesting.smk'
include: 'rules/fwer.smk'

rule yaml:
    input:
        macro+'/roi.case.tsv',
        macro+'/zhistogram.case.png',
        macro+'/zhistogram.control.png',
        macro+'/zhistogram.diff.png',
        macro+'/autocovariance.diff.png',
        macro+'/autocovariance.case.png',
        macro+'/autocovariance.control.png',
        macro+'/fwer.crosscovariance.diff.tsv',
        macro+'/scan.case.control.png',
    output:
        yaml=macro+'/arguments.case.yaml',
    params:
        yaml=str(config['CHANGE']['FOLDERS']['YAML']),
    shell:
        'cp {params.yaml} {output.yaml}'

# snakemake all -c1 -n
rule all:
    input:
        macro+'/arguments.case.yaml',
