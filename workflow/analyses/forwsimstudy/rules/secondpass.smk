wildcard_constraints:
	SIMNAME = '\w+',

### focus VCF file for short IBD ###

rule focus_region:
	input:
		windowin='{macro}/{micro}/firstpass.{process}.window.tsv.gz',
	output:
		vcfout='{macro}/{micro}/short.{process}.vcf.gz',
		btbi='{macro}/{micro}/large.{process}.vcf.bgz.tbi',
		bvcf='{macro}/{micro}/large.{process}.vcf.bgz',
	params:
		vcfin='{macro}/{micro}/large',
		process='{process}',
	shell:
		'bash ../../snakescripts/peak-vcf.sh {params.vcfin} {input.windowin} ../snakescripts/babysnakescripts {config[FIXED][PM]} {output.vcfout} {config[ISWEEPPARAMS][IBDCOMMA]} {params.process}'

### short IBD calls for focus region ###

n=int(float(config['FIXED']['SAMPSIZE']))
ploidy=int(float(config['FIXED']['PLOIDY']))
maf3=float(config['THIRDHAP']['MINMAF'])
mac3=int(ploidy*n*maf3)

rule short_hapibd:
    input:
        vcf='{macro}/{micro}/short.{process}.vcf.gz',
        map=str(config["FOLDERS"]["MACRO"]) + '/uniform.map',
    params:
        minmac=str(mac3),
        out='{macro}/{micro}/short.{process}',
    output:
        ibd='{macro}/{micro}/short.{process}.ibd.gz',
        hbd='{macro}/{micro}/short.{process}.hbd.gz',
        log='{macro}/{micro}/short.{process}.log',
    shell:
        'java -jar {config[PATHS][HAPIBD]} gt={input.vcf} map={input.map} out={params.out} min-seed={config[THIRDHAP][MINSEED]} min-extend={config[THIRDHAP][MINEXT]} min-output={config[THIRDHAP][MINOUT]} min-markers={config[THIRDHAP][MINMARK]} min-mac={params.minmac}'

### long IBD ###

rule long_ibd_true:
    input:
        ibd='{macro}/{micro}/firstpass.{process}.ibd.gz',
    output:
        ibd='{macro}/{micro}/secondpass.{process}.ibd.gz',
    shell:
        'zcat {input.ibd} | java -jar {config[PATHS][FILTER]} "I" -8 0.00 {config[ISWEEPPARAMS][IBDCUTOFF]} | java -jar {config[PATHS][FILTER]} "I" 6 0 {config[FIXED][LOC]} | java -jar {config[PATHS][FILTER]} "I" 7 {config[FIXED][LOC]} 100000000 | gzip > {output.ibd}'

### short IBD ###

rule short_ibd_true:
    input:
        ibd='{macro}/{micro}/short.{process}.ibd.gz',
    output:
        ibd='{macro}/{micro}/secondpass.{process}.ibd.gz',
    shell:
        'zcat {input.ibd} | java -jar {config[PATHS][FILTER]} "I" -8 0.00 {config[THIRDHAP][MINOUT]} | java -jar {config[PATHS][FILTER]} "I" 6 0 {config[FIXED][LOC]} | java -jar {config[PATHS][FILTER]} "I" 7 {config[FIXED][LOC]} 100000000 | gzip > {output.ibd}'
