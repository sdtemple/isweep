# setup macro folder
import os
macro=str(config['MACRO'])
if not os.path.exists(macro):
	os.mkdir(macro)

# load in experiments, set up micro, simname folders
import pandas as pd
micro=str(config["MICRO"])
sims = pd.read_csv(micro, sep='\t', header=0)
J = sims.shape[0]
for j in range(J):
	row = sims.loc[j,]
	if not os.path.exists(macro+'/'+str(row.MICROEXP)):
		os.mkdir(macro+'/'+str(row.MICROEXP))
	if not os.path.exists(macro+'/'+str(row.MICROEXP)+'/'+str(row.SIMNAME)):
		os.mkdir(macro+'/'+str(row.MICROEXP)+'/'+str(row.SIMNAME))
sims['FOLDER'] = [macro + '/' + sims.loc[j].MICROEXP + '/' + str(sims.loc[j].SIMNAME) for j in range(J)]
sims['FILE'] = sims['FOLDER'] + '/second.ranks.tsv.gz'
sims['EXISTS'] = sims['FILE'].apply(os.path.isfile)
sims=sims[sims['EXISTS']]
sims = sims.set_index("SIMNAME", drop=False)

rule all:
	input:
		[f"{sim.FOLDER}/imagene.results.tsv" for sim in sims.itertuples()],

rule unzip:
    input:
        vcf='{macro}/{micro}/{seed}/relate.chr1.vcf.gz',
    output:
        out='{macro}/{micro}/{seed}/relate.chr1.vcf',
    shell:
        """
        gunzip -c {input.vcf} > {output.out}
        """

rule imagene:
    input:
        vcf='{macro}/{micro}/{seed}/relate.chr1.vcf',
    output:
        out='{macro}/{micro}/{seed}/imagene.results.tsv',
    params:
       sampsize=str(config['SAMPSIZE']),
       effsize=str(config['EFFSIZE']),
       script=str(config['SCRIPT']),
       model=str(config['MODEL']),
    shell:
        """
        python {params.script} \
            {params.model} \
            {output.out} \
            {input.vcf} \
            {params.sampsize} \
            {params.effsize}
        """