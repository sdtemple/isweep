wildcard_constraints:
	SIMNAME = '\w+',

n=int(float(config['CHANGE']['SIMULATE']['SAMPSIZE']))
ploidy=int(float(config['FIXED']['HAPIBD']['PLOIDY']))
maf3=float(config['FIXED']['HAPIBD']['MINMAF'])
mac3=int(ploidy*n*maf3)

# iSWEEP simulation studies
# Seth D. Temple, sdtemple@uw.edu
# July 8, 2023

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
        [f"{sim.FOLDER}/isweep.ranks.tsv.gz" for sim in sims.itertuples()],


rule fair_region:
    input:
        vcfin='{macro}/{micro}/{seed}/large.chr1.vcf.gz',
    output:
        vcfout='{macro}/{micro}/{seed}/fair.chr1.vcf.gz',
    params:
        folder='{macro}/{micro}/{seed}',
        pm=str(config['FIXED']['SIMULATE']['BUFFER']),
        script=str(config['CHANGE']['FOLDERS']['TERMINALSCRIPTS'])+'/lines.py',
    shell: # to bgz and back is being consertative
        """
        gunzip -c {input.vcfin} | bgzip  > {params.folder}/chrtemp.vcf.bgz
        tabix -fp vcf {params.folder}/chrtemp.vcf.bgz
        left=2800000
        right=5200000
        bcftools view {params.folder}/chrtemp.vcf.bgz \
            -r 1:${{left}}-${{right}} \
            -Oz -o {output.vcfout}
        rm {params.folder}/chrtemp.vcf.bgz
        """

### call hap-ibd ###

rule fair_hapibd:
    input:
        vcf='{macro}/{micro}/{seed}/fair.chr1.vcf.gz',
        map=str(config["CHANGE"]["FOLDERS"]["MACRO"]) + '/uniform.map',
    params:
        minmac=str(mac3),
        out='{macro}/{micro}/{seed}/fair.chr1',
        soft=str(config['CHANGE']['FOLDERS']['SOFTWARE']),
        prog=str(config['CHANGE']['PROGRAMS']['HAPIBD']),
        minsee=str(config['FIXED']['HAPIBD']['MINSEED']),
        minext=str(config['FIXED']['HAPIBD']['MINEXT']),
        minout=str(config['FIXED']['HAPIBD']['MINOUT']),
    output:
        ibd='{macro}/{micro}/{seed}/fair.chr1.ibd.gz',
        hbd='{macro}/{micro}/{seed}/fair.chr1.hbd.gz',
        log='{macro}/{micro}/{seed}/fair.chr1.log',
    resources:
        mem_gb='{config[CHANGE][CLUSTER][LARGEMEM]}',
    shell:
        """
        java -Xmx{config[CHANGE][CLUSTER][LARGEMEM]}g -jar {params.soft}/{params.prog} \
            gt={input.vcf} \
            map={input.map} \
            out={params.out} \
            min-seed={params.minsee} \
            min-extend={params.minext} \
            min-output={params.minout} \
            min-mac={params.minmac}
        """

### filter ibd file ###

rule fair_filt:
    input:
        ibd='{macro}/{micro}/{seed}/fair.chr1.ibd.gz',
    output:
        ibd='{macro}/{micro}/{seed}/fair.filt.chr1.ibd.gz',
    params:
        soft=str(config['CHANGE']['FOLDERS']['SOFTWARE']),
        prog=str(config['CHANGE']['PROGRAMS']['FILTER']),
    resources:
        mem_gb='{config[CHANGE][CLUSTER][LARGEMEM]}'
    shell:
        """
        thecenter=4000000
        zcat {input.ibd} | \
            java -Xmx{config[CHANGE][CLUSTER][LARGEMEM]}g -jar {params.soft}/{params.prog} \
            "I" 6 0.00 $thecenter | \
            java -Xmx{config[CHANGE][CLUSTER][LARGEMEM]}g -jar {params.soft}/{params.prog} \
            "I" 7 $thecenter 10000000000 | \
            gzip > {output.ibd}
        """

### rank snps ###

rule fair_rank:
    input:
        short='{macro}/{micro}/{seed}/fair.filt.chr1.ibd.gz',
        vcf='{macro}/{micro}/{seed}/fair.chr1.vcf.gz',
    output:
        fileout='{macro}/{micro}/{seed}/isweep.ranks.tsv.gz',
    params:
        scripts=str(config['CHANGE']['FOLDERS']['TERMINALSCRIPTS']),
        diameter=str(config['FIXED']['ISWEEP']['DIAMETER']),
        q1=str(config['FIXED']['ISWEEP']['MINAAF1']),
        rulesigma=str(config['FIXED']['ISWEEP']['RULESIGMA']),
    shell:
        """
        python {params.scripts}/rank.py \
            {input.short} \
            {input.vcf} \
            {output.fileout} \
            {params.diameter} \
            {params.q1} \
            {params.rulesigma}
        """
