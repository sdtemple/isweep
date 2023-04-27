wildcard_constraints:
	SIMNAME = '\w+',

### focus VCF file for short IBD ###

rule focus_region:
	input:
		vcfin='{macro}/{micro}/{seed}/large.chr1.vcf.gz',
	output:
		subvcf='{macro}/{micro}/{seed}/small.chr1.vcf.gz',
	shell:
		'bash {config[CHANGE][FOLDERS][TERMINALSCRIPTS]}/roi-vcf.sh {macro}/{micro}/{seed} {input.vcfin} {output.vcfout} {config[FIXED][SIMULATE][LEFT]} {config[FIXED][SIMULATE][RIGHT]} 1 {config[CHANGE][SIMULATE][MSPMAF]}'

### short IBD calls for focus region ###

n=int(float(config['CHANGE']['SIMULATE']['SAMPSIZE']))
ploidy=int(float(config['FIXED']['HAPIBD']['PLOIDY']))
maf3=float(config['FIXED']['HAPIBD']['MINMAF'])
mac3=int(ploidy*n*maf3)
largethreads=int(float(str(config['CHANGE']['CLUSTER']['LARGETHREAD'])))

rule small_hapibd:
    input:
        vcf='{macro}/{micro}/{seed}/small.chr1.vcf.gz',
        map=str(config["FOLDERS"]["MACRO"]) + '/uniform.map',
    params:
        minmac=str(mac3),
        out='{macro}/{micro}/{seed}/small.chr1',
    output:
        ibd='{macro}/{micro}/{seed}/small.chr1.ibd.gz',
        hbd='{macro}/{micro}/{seed}/small.chr1.hbd.gz',
        log='{macro}/{micro}/{seed}/small.chr1.log',
	threads:
		largethreads,
	resources:
		mem_gb='{config[CHANGE][CLUSTER][LARGEMEM]}',
    shell:
        'java -Xmx{config[CHANGE][CLUSTER][LARGEMEM]}g -jar {config[CHANGE][FOLDERS][SOFTWARE]}/{config[CHANGE][PROGRAMS][HAPIBD]} gt={input.vcf} map={input.map} out={params.out} min-seed={config[FIXED][HAPIBD][MINSEED]} min-extend={config[FIXED][HAPIBD][MINEXT]} min-output={config[FIXED][HAPIBD][MINOUT]} min-markers={config[FIXED][HAPIBD][MINMARK]} min-mac={params.minmac}'

### long IBD ###

rule long_ibd:
    input:
        ibd='{macro}/{micro}/{seed}/scan.chr1.ibd.gz',
    output:
        ibd='{macro}/{micro}/{seed}/long.chr1.ibd.gz',
    shell:
        'zcat {input.ibd} | java -jar {config[CHANGE][FOLDERS][SOFTWARE]}/{config[CHANGE][PROGRAMS][FILTER]} "I" -8 0.00 {config[FIXED][ISWEEP][MOMCUTOFF]} | java -jar {config[CHANGE][FOLDERS][SOFTWARE]}/{config[CHANGE][PROGRAMS][FILTER]} "I" 6 0 {config[FIXED][SIMULATE][LOC]} | java -jar {config[CHANGE][FOLDERS][SOFTWARE]}/{config[CHANGE][PROGRAMS][FILTER]} "I" 7 {config[FIXED][SIMULATE][LOC]} 100000000 | gzip > {output.ibd}'

### small IBD ###

rule small_ibd:
    input:
        ibd='{macro}/{micro}/{seed}/small.chr1.ibd.gz',
    output:
        ibd='{macro}/{micro}/{seed}/short.chr1.ibd.gz',
    shell:
        'zcat {input.ibd} | java -jar {config[CHANGE][FOLDERS][SOFTWARE]}/{config[CHANGE][PROGRAMS][FILTER]} "I" -8 0.00 {config[FIXED][HAPIBD][MINOUT]} | java -jar {config[CHANGE][FOLDERS][SOFTWARE]}/{config[CHANGE][PROGRAMS][FILTER]} "I" 6 0 {config[FIXED][SIMULATE][LOC]} | java -jar {config[CHANGE][FOLDERS][SOFTWARE]}/{config[CHANGE][PROGRAMS][FILTER]} "I" 7 {config[FIXED][SIMULATE][LOC]} 100000000 | gzip > {output.ibd}'
