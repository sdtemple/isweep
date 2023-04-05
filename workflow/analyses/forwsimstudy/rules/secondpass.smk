wildcard_constraints:
	SIMNAME = '\w+',

### focus VCF file for short IBD ###

rule focus_region:
	input:
		windowin='{macro}/{micro}/firstpass.{process}.window.tsv.gz',
	output:
		vcfout='{macro}/{micro}/small.{process}.vcf.gz',
		btbi='{macro}/{micro}/large.{process}.vcf.bgz.tbi',
		bvcf='{macro}/{micro}/large.{process}.vcf.bgz',
	params:
		vcfin='{macro}/{micro}/large',
		process='{process}',
	shell:
		'bash ../../snakescripts/peak-vcf.sh {params.vcfin} {input.windowin} ../../snakescripts {config[FIXED][PM]} {output.vcfout} {config[ISWEEPPARAMS][IBDCOMMA]} {params.process}'

### short IBD calls for focus region ###

n=int(float(config['FIXED']['SAMPSIZE']))
ploidy=int(float(config['FIXED']['PLOIDY']))
maf3=float(config['THIRDHAP']['MINMAF'])
mac3=int(ploidy*n*maf3)

rule small_hapibd:
    input:
        vcf='{macro}/{micro}/small.{process}.vcf.gz',
        map=str(config["FOLDERS"]["MACRO"]) + '/uniform.map',
    params:
        minmac=str(mac3),
        out='{macro}/{micro}/small.{process}',
    output:
        ibd='{macro}/{micro}/small.{process}.ibd.gz',
        hbd='{macro}/{micro}/small.{process}.hbd.gz',
        log='{macro}/{micro}/small.{process}.log',
    shell:
        'java -jar {config[PATHS][HAPIBD]} gt={input.vcf} map={input.map} out={params.out} min-seed={config[THIRDHAP][MINSEED]} min-extend={config[THIRDHAP][MINEXT]} min-output={config[THIRDHAP][MINOUT]} min-markers={config[THIRDHAP][MINMARK]} min-mac={params.minmac}'

### long IBD ###

rule long_ibd:
    input:
        ibd='{macro}/{micro}/firstpass.{process}.ibd.gz',
    output:
        ibd='{macro}/{micro}/long.{process}.ibd.gz',
    shell:
        'zcat {input.ibd} | java -jar {config[PATHS][FILTER]} "I" -8 0.00 {config[ISWEEPPARAMS][IBDCUTOFF]} | java -jar {config[PATHS][FILTER]} "I" 6 0 {config[FIXED][LOC]} | java -jar {config[PATHS][FILTER]} "I" 7 {config[FIXED][LOC]} 100000000 | gzip > {output.ibd}'

### small IBD ###

rule small_ibd:
    input:
        ibd='{macro}/{micro}/small.{process}.ibd.gz',
    output:
        ibd='{macro}/{micro}/short.{process}.ibd.gz',
    shell:
        'zcat {input.ibd} | java -jar {config[PATHS][FILTER]} "I" -8 0.00 {config[THIRDHAP][MINOUT]} | java -jar {config[PATHS][FILTER]} "I" 6 0 {config[FIXED][LOC]} | java -jar {config[PATHS][FILTER]} "I" 7 {config[FIXED][LOC]} 100000000 | gzip > {output.ibd}'
