wildcard_constraints:
	SIMNAME = '\w+',

rule uniform_map:
	input:
	output:
		mapout=str(config["CHANGE"]["FOLDERS"]["MACRO"]) + '/uniform.map',
	script:
		'{config[CHANGE][FOLDERS][SNAKESCRIPTS]}/uniformMap.py'

rule forward_Ne:
    input: str(config["CHANGE"]["FOLDERS"]["MACRO"]) + '/uniform.map',
	output: str(config["CHANGE"]["FOLDERS"]["MACRO"]) + '/forward.ne',
	script:
		'{config[CHANGE][FOLDERS][SNAKESCRIPTS]}/forwardNe.py'

rule slim_script:
	input:
		trueNe = str(config["CHANGE"]["FOLDERS"]["MACRO"]) + '/forward.ne',
		map = str(config["CHANGE"]["FOLDERS"]["MACRO"]) + '/uniform.map',
	output:
		[f"{sim.FOLDER}/slimulation.slim" for sim in sims.itertuples()],
	script:
		'{config[CHANGE][FOLDERS][SNAKESCRIPTS]}/writeSlimDemography.py'

### generate tree, vcf data ###

# slim forward
rule slim:
	input:
		"{macro}/{micro}/{seed}/slimulation.slim",
	output:
		trees="{macro}/{micro}/{seed}/slimulation.trees",
		freq="{macro}/{micro}/{seed}/slimulation.freq",
	shell:
		'{config[CHANGE][FOLDERS][SOFTWARE]}/{config[CHANGE][PROGRAMS][SLiM]} {wildcards.macro}/{wildcards.micro}/{wildcards.seed}/slimulation.slim'

# msprime backward
rule msprime:
	input:
		trees = "{macro}/{micro}/{seed}/slimulation.trees",
	output:
		bcf="{macro}/{micro}/{seed}/slimulation.bcf.gz",
	script:
		'{config[CHANGE][FOLDERS][SNAKESCRIPTS]}/treeVCF.py'

# bcf, vcf magic
rule tabix_bcf:
	input:
		bcf='{macro}/{micro}/{seed}/slimulation.bcf.gz',
	output:
		csi='{macro}/{micro}/{seed}/slimulation.bcf.gz.csi',
	shell:
		'tabix {input.bcf}'

rule vcf:
	input:
		bcf='{macro}/{micro}/{seed}/slimulation.bcf.gz',
		csi='{macro}/{micro}/{seed}/slimulation.bcf.gz.csi',
	output:
		vcf='{macro}/{micro}/{seed}/slimulation.vcf.gz',
	shell:
		'bcftools view {input.bcf} -Oz -m2 -M2 -o {output.vcf}'

rule genotyping_error:
	input:
		vcf='{macro}/{micro}/{seed}/slimulation.vcf.gz',
	output:
		out='{macro}/{micro}/{seed}/large.chr1.vcf.gz',
	shell:
		'zcat {input.vcf} | java -Xmx{config[CHANGE][PROGRAMS][XMXMEM]}g -jar {config[CHANGE][FOLDERS][SOFTWARE]}/{config[CHANGE][PROGRAMS][GTERR]} {config[FIXED][SIMULATE][GTERR]} | gzip > {output.out}'

### true labels ###

# work on location
rule causal_vcf:
	input:
		bcf='{macro}/{micro}/{seed}/slimulation.bcf.gz',
		csi='{macro}/{micro}/{seed}/slimulation.bcf.gz.csi'
	output:
		vcf='{macro}/{micro}/{seed}/causal.chr1.vcf.gz'
	shell:
		'bcftools view {input.bcf} -Oz -m2 -M2 -t 1:{config[FIXED][SIMULATE][LOC]} -o {output.vcf}'

rule true_labels:
	input:
		vcf='{macro}/{micro}/{seed}/causal.chr1.vcf.gz'
	output:
		labels1='{macro}/{micro}/{seed}/slimulation.true1.labels.gz',
		labels0='{macro}/{micro}/{seed}/slimulation.true0.labels.gz'
	script:
		'{config[CHANGE][FOLDERS][SOFTWARE]}/{config[CHANGE][PROGRAMS][SNAKESCRIPTS]}/true-labels.py'

### zip tree sequence ###

rule tszip:
	input:
		trees='{macro}/{micro}/{seed}/slimulation.trees',
		yaml='{macro}/{micro}/{seed}/slimulation.bcf.gz',
	output:
		tsz='{macro}/{micro}/{seed}/slimulation.trees.tsz'
	shell:
		'tszip {input.trees}'

### write yaml as log ###

rule record:
	input:
		'{macro}/{micro}/{seed}/slimulation.slim'
	output:
		yaml='{macro}/{micro}/{seed}/arguments.yaml'
	shell:
		'cp ./simulatey.yaml {output.yaml}'
