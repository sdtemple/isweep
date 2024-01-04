# simulating data
# run snakemake -s Snakefile-simulate.smk --configfile *.yaml

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

include: 'rules/simulate.smk'

rule all:
    input:
        [f"{sim.FOLDER}/slimulation.trees.tsz".replace(" ","") for sim in sims.itertuples()],
        [f"{sim.FOLDER}/large.chr1.vcf.gz".replace(" ","") for sim in sims.itertuples()],
