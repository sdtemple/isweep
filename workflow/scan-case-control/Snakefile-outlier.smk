# The main organizing file for the outlier detection + phenotypes,
# where the rule all is and many of the yaml parameters.
# You run this file with the -s parameter in snakemake.

import os
import shutil
import pandas as pd
macro=str(config['change']['files']['study'])
micro=str(config['change']["files"]["regions_of_interest"])

folder_name = macro
if not os.path.exists(folder_name):
    os.makedirs(folder_name)
if not os.path.exists(macro+'/'+micro):
    shutil.copy(micro,macro+'/'+micro)

sims = pd.read_csv(macro+'/'+micro, sep='\t', header=0)
J = sims.shape[0]
for j in range(J):
	row = sims.iloc[j]
	if not os.path.exists(macro+'/'+str(row.NAME)):
		os.mkdir(macro+'/'+str(row.NAME))
	f=open(macro+'/'+str(row.NAME)+'/locus.case.txt','w')
	f.write('NAME\t')
	f.write(str(row.NAME))
	f.write('\n')
	f.write('CHROM\t')
	f.write(str(int(row.CHROM)))
	f.write('\n')
	f.write('CENTER\t')
	f.write(str(int(row.BPCENTER)))
	f.write('\n')
	f.write('LEFT\t')
	f.write(str(int(row.BPLEFTCENTER)))
	f.write('\n')
	f.write('RIGHT\t')
	f.write(str(int(row.BPRIGHTCENTER)))
	f.close()
sims['FOLDER'] = [(macro +'/'+str(sims.iloc[j].NAME)).strip() for j in range(J)]

# copy and paste maps
folder_name = macro + '/maps'
if not os.path.exists(folder_name):
    os.makedirs(folder_name)

mapfol=str(config['change']['files']['genetic_maps'])
mappre=str(config['change']['files']['map_prefix'])
mapsuf=str(config['change']['files']['map_suffix'])

low=str(config['change']['files']['chromosome_low'])
high=str(config['change']['files']['chromosome_high'])
low=int(low)
high=int(high)

for i in range(low,high+1):
    source_file = mapfol+'/'+mappre+str(i)+mapsuf
    destination_file = macro+'/maps/chr'+str(i)+'.map'
    if not os.path.exists(destination_file):
        shutil.copy(source_file, destination_file)

cases=str(config['change']['files']['cases'])
file_name = macro + '/phenotypes.txt'
if not os.path.exists(file_name):
    shutil.copy(cases, file_name)


# snakemake all -c1 -n
rule all:
	input:
		[(macro +'/'+str(sims.iloc[j].NAME)).strip()+'/matrix.outlier.phenotypes.tsv' for j in range(J)],
	output:
		yaml=macro+'/arguments.outliers.yaml',
	params:
		yaml=str(config['change']['files']['yaml']),
	shell:
		'cp {params.yaml} {output.yaml}'


# zooming in on a region of interest

# some inputs, string managements, count sample size
subsamplefile=cases
cohort=str(config['change']['files']['study'])
samplesize=0
with open(subsamplefile,'r') as f:
    for line in f:
        samplesize+=1
ploidy=int(float(config['change']['files']['ploidy']))
maf3=float(config['fixed']['hap_ibd']['min_minor_allele_frequency'])
mac3=int(ploidy*samplesize*maf3)

# subset vcf to region of interest
rule first_region: # focus vcf on region of interest
    input:
        locus='{cohort}/{hit}/locus.case.txt',
        subsample="{cohort}/subsample.txt",
    output:
        subvcf='{cohort}/{hit}/narrowing.case.vcf.gz',
    params:
        qmaf=maf3,
        chrpre=str(config['change']['files']['chromosome_prefix']),
        vcfs=str(config['change']['files']['vcfs']),
        vcfprefix=str(config['change']['files']['vcf_prefix']),
        vcfsuffix=str(config['change']['files']['vcf_suffix']),
    resources:
        mem_gb='{config[change][isweep][xmx_mem]}',
    shell: # if chromosome is huge (greater than 10000 Mb), may need to modify the third pipe
        """
        chr=$(python ../../scripts/utilities/lines.py {input.locus} 2 2)
        left=$(python ../../scripts/utilities/lines.py {input.locus} 4 2)
        right=$(python ../../scripts/utilities/lines.py {input.locus} 5 2)
        vcf={params.vcfs}/{params.vcfprefix}${{chr}}{params.vcfsuffix}
        tabix -fp vcf $vcf
        bcftools view ${{vcf}} -r {params.chrpre}${{chr}}:${{left}}-${{right}} -Ob | \
            bcftools view -S {input.subsample} -Ob | \
            bcftools view -q {params.qmaf}:nonmajor -Oz -o {output.subvcf}
        """

### call hap-ibd ###
rule first_hapibd:
    input:
        vcf='{cohort}/{hit}/narrowing.case.vcf.gz',
        locus='{cohort}/{hit}/locus.case.txt',
    params:
        minmac=str(mac3),
        out='{cohort}/{hit}/narrowing.case',
        minsee=str(config['fixed']['hap_ibd']['min_seed']),
        minext=str(config['fixed']['hap_ibd']['min_extend']),
        minout=str(config['fixed']['hap_ibd']['min_output']),
    output:
        ibd='{cohort}/{hit}/narrowing.case.ibd.gz',
        hbd='{cohort}/{hit}/narrowing.case.hbd.gz',
        log='{cohort}/{hit}/narrowing.case.log',
    resources:
        mem_gb='{config[change][isweep][xmx_mem]}',
    shell:
        """
        chr=$(python ../../scripts/utilities/lines.py {input.locus} 2 2)
        java -Xmx{config[change][isweep][xmx_mem]}g -jar ../../software/hap-ibd.jar \
            gt={input.vcf} \
            map={wildcards.cohort}/maps/chr${{chr}}.map \
            out={params.out} \
            min-seed={params.minsee} \
            min-extend={params.minext} \
            min-output={params.minout} \
            min-mac={params.minmac}
        """

### filter ibd file ###
rule first_filt:
    input:
        ibd='{cohort}/{hit}/narrowing.case.ibd.gz',
        locus='{cohort}/{hit}/locus.case.txt',
    output:
        ibd='{cohort}/{hit}/narrowing.filt.case.ibd.gz',
    shell:
        """
        thecenter=$(python ../../scripts/utilities/lines.py {input.locus} 3 2)
        python ../../scripts/utilities/filter-lines.py \
            --input_file {input.ibd} \
            --output_file {wildcards.cohort}/{wildcards.hit}/intermediate.ibd.gz \
            --column_index 6 \
            --upper_bound $thecenter \
            --complement 0
        python ../../scripts/utilities/filter-lines.py \
            --input_file {wildcards.cohort}/{wildcards.hit}/intermediate.ibd.gz \
            --output_file {output.ibd} \
            --column_index 7 \
            --lower_bound $thecenter \
            --upper_bound 10000000000 \
            --complement 0
        rm {wildcards.cohort}/{wildcards.hit}/intermediate.ibd.gz
        """


### write outliers ###
rule outlier:
    input:
        short='{cohort}/{hit}/narrowing.filt.case.ibd.gz',
    output:
        out1='{cohort}/{hit}/outlier1.phenotype.tsv',
    params:
        diameter=str(config['fixed']['isweep']['diameter']),
        rulesigma=str(config['fixed']['isweep']['group_cutoff']),
        cases=cases,
    shell:
        """
        python ../../scripts/model/outliers-phenotype.py \
            --input_ibd_file {input.short} \
            --input_phenotype_file {params.cases} \
            --output_folder {wildcards.cohort}/{wildcards.hit} \
            --graph_diameter {params.diameter} \
            --group_cutoff {params.rulesigma}
        """

rule design:
    input:
        short='{cohort}/{hit}/outlier1.phenotype.tsv',
    output:
        out1='{cohort}/{hit}/matrix.outlier.phenotypes.tsv',
        out2='{cohort}/{hit}/matrix.outlier.phenotypes.sorted.tsv',
    params:
        diameter=str(config['fixed']['isweep']['diameter']),
        rulesigma=str(config['fixed']['isweep']['group_cutoff']),
        cases=cases,
        ploidy=str(ploidy),
    shell:
        """
        python ../../scripts/utilities/make-design-matrix.py \
            --output_design_matrix_prefix {wildcards.cohort}/{wildcards.hit}/matrix.outlier.phenotypes \
            --input_phenotype_file {params.cases} \
            --input_outlier_prefix {wildcards.cohort}/{wildcards.hit}/outlier \
            --input_outlier_suffix .phenotype.tsv \
            --first_index 1 \
            --ploidy {params.ploidy}
        """

