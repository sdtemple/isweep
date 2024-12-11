import os
import pandas as pd
macro=str(config['CHANGE']['FOLDERS']['STUDY'])
micro=str(config['CHANGE']["ISWEEP"]["ROI"])
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

# snakemake all -c1 -n
rule all:
	input:
		[(macro +'/'+str(sims.iloc[j].NAME)).strip()+'/outlier1.phenotype.txt' for j in range(J)],
	output:
		yaml=macro+'/arguments.outliers.yaml',
	params:
		yaml=str(config['CHANGE']['FOLDERS']['YAML']),
	shell:
		'cp {params.yaml} {output.yaml}'


# zooming in on a region of interest

# some inputs, string managements, count sample size
subsamplefile=str(config['CHANGE']['ISWEEP']['CASES'])
cohort=str(config['CHANGE']['FOLDERS']['STUDY'])
samplesize=0
with open(subsamplefile,'r') as f:
    for line in f:
        samplesize+=1
ploidy=2
# ploidy=int(float(config['FIXED']['HAPIBD']['PLOIDY']))
maf3=float(config['FIXED']['HAPIBD']['MINMAF'])
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
        chrpre=str(config['CHANGE']['ISWEEP']['CHRPRE']),
        vcfs=str(config['CHANGE']['EXISTING']['VCFS']),
        vcfprefix=str(config['CHANGE']['EXISTING']['VCFPRE']),
        vcfsuffix=str(config['CHANGE']['EXISTING']['VCFSUF']),
    resources:
        mem_gb='{config[CHANGE][ISWEEP][XMXMEM]}',
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
        minsee=str(config['FIXED']['HAPIBD']['MINSEED']),
        minext=str(config['FIXED']['HAPIBD']['MINEXT']),
        minout=str(config['FIXED']['HAPIBD']['MINOUT']),
    output:
        ibd='{cohort}/{hit}/narrowing.case.ibd.gz',
        hbd='{cohort}/{hit}/narrowing.case.hbd.gz',
        log='{cohort}/{hit}/narrowing.case.log',
    resources:
        mem_gb='{config[CHANGE][ISWEEP][XMXMEM]}',
    shell:
        """
        chr=$(python ../../scripts/utilities/lines.py {input.locus} 2 2)
        java -Xmx{config[CHANGE][ISWEEP][XMXMEM]}g -jar ../../software/hap-ibd.jar \
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
        out1='{cohort}/{hit}/outlier1.phenotype.txt',
    params:
        diameter=str(config['FIXED']['ISWEEP']['DIAMETER']),
        rulesigma=str(config['FIXED']['ISWEEP']['GROUPCUTOFF']),
        cases=str(config['CHANGE']['ISWEEP']['CASES']),
    shell:
        """
        python ../../scripts/model/outliers-phenotype.py \
            --input_ibd_file {input.short} \
            --input_phenotype_file {params.cases} \
            --output_folder {wildcards.cohort}/{wildcards.hit} \
            --graph_diameter {params.diameter} \
            --group_cutoff {params.rulesigma}
        touch {output.out1}
        """

# ### map the phenotypes to the outlier haplotypes
# rule phenotypes:
#     input:
#         inn='{cohort}/{hit}/outlier{num}.case.txt',
#     output:
#         out='{cohort}/{hit}/outlier{num}.case.phenotype.txt',
#     params:
#         cases=str(config['CHANGE']['ISWEEP']['CASES']),
#     shell:
#         """
#         python ../../scripts/model/map-phenotype.py \
#             {input.inn} \
#             {params.cases} \
#             {output.out} \
#         """

