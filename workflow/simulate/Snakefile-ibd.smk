# iSWEEP simulation studies
# Seth D. Temple, sdtemple@uw.edu
# November 17, 2023
# editing to do ibd analysis after simulating data
# you have to have first run snakemake -s Snakefile-simulate.smk

# setup macro folder
import os
macro=str(config['CHANGE']['FOLDERS']['MACRO'])
if not os.path.exists(macro):
	os.mkdir(macro)

# load in experiments, set up micro, simname folders
import pandas as pd
micro=str(config['CHANGE']["FOLDERS"]["MICRO"])
sims = pd.read_csv(micro, sep='\t', header=0)
J = sims.shape[0]
for j in range(J):
	row = sims.loc[j,]
	if not os.path.exists(macro+'/'+str(row.MICROEXP)):
		os.mkdir(macro+'/'+str(row.MICROEXP))
	if not os.path.exists(macro+'/'+str(row.MICROEXP)+'/'+str(row.SIMNAME)):
		os.mkdir(macro+'/'+str(row.MICROEXP)+'/'+str(row.SIMNAME))
sims['FOLDER'] = [macro + '/' + sims.loc[j].MICROEXP + '/' + str(sims.loc[j].SIMNAME) for j in range(J)]
sims = sims.set_index("SIMNAME", drop=False)

include: 'rules/scan.smk'
include: 'rules/first.smk'
include: 'rules/second.smk'
include: 'rules/third.smk'

rule yaml:
    input:
        [f"{sim.FOLDER}/results.snp.tsv" for sim in sims.itertuples()],
        [f"{sim.FOLDER}/results.hap.tsv" for sim in sims.itertuples()],
        [f"{sim.FOLDER}/ibd.entropy.tsv" for sim in sims.itertuples()],
    output:
        yaml=macro+'/arguments.simulate.yaml',
    params:
        yaml=str(config['CHANGE']['FOLDERS']['YAML']),
    shell:
        'cp {params.yaml} {output.yaml}'

rule all:
    input:
        yaml=macro+'/arguments.simulate.yaml',