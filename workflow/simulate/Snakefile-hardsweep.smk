wildcard_constraints:
	SIMNAME = '\w+',

n=int(float(config['CHANGE']['SIMULATE']['SAMPSIZE']))
ploidy=int(float(config['CHANGE']['SIMULATE']['PLOIDY']))
sample_size = n * ploidy

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

rule all:
    input:
        [f"{sim.FOLDER}/if.hard.sweep.txt" for sim in sims.itertuples()],

rule if_hard_sweep:
    input:
        filedone='{macro}/{micro}/{seed}/results.hap.tsv',
    output:
        fileout='{macro}/{micro}/{seed}/if.hard.sweep.txt',
    params:
        folder='{macro}/{micro}/{seed}',
        scripts=str(config['CHANGE']['FOLDERS']['TERMINALSCRIPTS']),
        sample_size=str(sample_size),
        percent_rule=str(config['FIXED']['HARDSWEEP']['PERCENTRULE']),
        plurality_rule=str(config['FIXED']['HARDSWEEP']['PLURALITYRULE']),
    shell:
        """
        python {params.scripts}/if-hard-sweep.py \
            {params.folder} \
            {output.fileout} \
            {params.percent_rule} \
            {params.plurality_rule} \
            {params.sample_size}
        """
