wildcard_constraints:
	SIMNAME = '\w+',

rule uniform_map:
	input:
	output:
		mapout=str(config["FOLDERS"]["MACRO"]) + '/uniform.map',
	script:
		'../../../snakescripts/uniformMap.py'

rule forward_Ne:
    input: str(config["FOLDERS"]["MACRO"]) + '/uniform.map',
	output: str(config["FOLDERS"]["MACRO"]) + '/forward.ne',
	script:
		'../../../snakescripts/forwardNe.py'

rule slim_script:
	input:
		trueNe = str(config["FOLDERS"]["MACRO"]) + '/forward.ne',
		map = str(config["FOLDERS"]["MACRO"]) + '/uniform.map',
	output:
		[f"{sim.FOLDER}/slimulation.slim" for sim in sims.itertuples()],
	script:
		'../../../snakescripts/writeSlimDemography.py'

### generate tree, vcf data ###

# slim forward

rule slim:
	input:
		"{macro}/{micro}/{seed}/slimulation.slim",
	output:
		trees="{macro}/{micro}/{seed}/slimulation.trees",
		freq="{macro}/{micro}/{seed}/slimulation.freq",
	shell:
		'{config[PATHS][SLiM]} {wildcards.macro}/{wildcards.micro}/{wildcards.seed}/slimulation.slim'

# msprime backward

rule msprime:
	input:
		trees = "{macro}/{micro}/{seed}/slimulation.trees",
	output:
		bcf="{macro}/{micro}/{seed}/slimulation.bcf.gz",
	script:
		'../../../snakescripts/treeVCF.py'

# bcf, vcf magic

rule tabix_bcf:
	input:
		bcf='{macro}/{micro}/slimulation.bcf.gz',
	output:
		csi='{macro}/{micro}/slimulation.bcf.gz.csi',
	shell:
		'tabix {input.bcf}'


rule vcf:
	input:
		bcf='{macro}/{micro}/slimulation.bcf.gz',
		csi='{macro}/{micro}/slimulation.bcf.gz.csi',
	output:
		vcf='{macro}/{micro}/slimulation.vcf.gz',
	shell:
		'bcftools view {input.bcf} -Oz -m2 -M2 -o {output.vcf}'

rule genotyping_error:
	input:
		vcf='{macro}/{micro}/slimulation.vcf.gz',
	output:
		out='{macro}/{micro}/large.vcf.gz',
	shell:
		'zcat {input.vcf} | java -jar {config[PATHS][GTERR]} {config[FIXED][GTERR]} | gzip > {output.out}'

rule remove_vcf_after_phasing:
	input:
		vcf='{macro}/{micro}/slimulation.vcf.gz',
		pvcf='{macro}/{micro}/large.vcf.gz',
	shell:
		'rm {input.vcf}'

### true labels ###

# work on location
rule causal_vcf:
	input:
		bcf='{macro}/{micro}/slimulation.bcf.gz',
		csi='{macro}/{micro}/slimulation.bcf.gz.csi'
	output:
		vcf='{macro}/{micro}/causal.vcf.gz'
	shell:
		'bcftools view {input.bcf} -Oz -m2 -M2 -t 1:{config[FIXED][LOC]} -o {output.vcf}'

rule true_labels:
	input:
		vcf='{macro}/{micro}/causal.vcf.gz'
	output:
		labels1='{macro}/{micro}/slimulation.true1.labels.gz',
		labels0='{macro}/{micro}/slimulation.true0.labels.gz'
	script:
		'../../../snakescripts/true-labels.py'
