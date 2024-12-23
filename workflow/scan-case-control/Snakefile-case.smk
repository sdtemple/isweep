# The main organizing file for the case-control scan,
# where the rule all is and many of the yaml parameters.
# You run this file with the -s parameter in snakemake.

import shutil
import os

# setup macro folder
macro=str(config['change']['files']['study'])
low=int(str(config['change']['files']['chromosome_low']))
high=int(str(config['change']['files']['chromosome_high']))
chrlow=low
chrhigh=high
cases=str(config['change']['files']['cases'])

# copy and paste maps
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
file_name = macro + '/phenotypes.txt'
if not os.path.exists(file_name):
    shutil.copy(cases, file_name)


include: 'rules/interesting.smk'
include: 'rules/fwer.smk'

# record the configuration settings used
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
        macro+'/fwer.simulation.case.txt',
    output:
        yaml=macro+'/arguments.case.yaml',
    params:
        yaml=str(config['change']['files']['yaml']),
    shell:
        'cp {params.yaml} {output.yaml}'

# snakemake all -c1 -n
rule all:
    input:
        macro+'/arguments.case.yaml',
