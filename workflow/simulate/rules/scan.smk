wildcard_constraints:
	SIMNAME = '\w+',

n=int(float(config['CHANGE']['ISWEEP']['SAMPSIZE']))
ploidy=int(float(config['FIXED']['CANDHAPIBD']['PLOIDY']))
maf1=float(config['FIXED']['CANDHAPIBD']['MINMAF'])
mac1=int(ploidy*n*maf1)

### ibd-ends from candidate segments ###

rule hapibd:
    input:
        vcf='{macro}/{micro}/{seed}/large.chr1.vcf.gz',
        map=str(config["FOLDERS"]["MACRO"]) + '/uniform.map',
    params:
        minmac=str(mac1),
        out='{macro}/{micro}/{seed}/large.chr1.hapibd.candidate',
    output:
        ibd='{macro}/{micro}/{seed}/large.chr1.hapibd.candidate.ibd.gz',
        hbd='{macro}/{micro}/{seed}/large.chr1.hapibd.candidate.hbd.gz',
        log='{macro}/{micro}/{seed}/large.chr1.hapibd.candidate.log',
	threads:
		'{config[CHANGE][CLUSTER][SMALLTHREAD]}',
	resources:
		mem_gb='{config[CHANGE][CLUSTER][SMALLMEM]}'
    shell:
        'java -Xmx{config[CHANGE][CLUSTER][SMALLMEM]}g -jar {config[CHANGE][FOLDERS][SOFTWARE]}/{config[CHANGE][PROGRAMS][HAPIBD]} gt={input.vcf} map={input.map} out={params.out} min-seed={config[FIXED][CANDHAPIBD][MINSEED]} min-extend={config[FIXED][CANDHAPIBD][MINEXT]} min-output={config[FIXED][CANDHAPIBD][MINOUT]} min-mac={params.minmac}'

rule ibdends:
    input:
        vcf='{macro}/{micro}/{seed}/large.chr1.vcf.gz',
        map=str(config["FOLDERS"]["MACRO"]) + '/uniform.map',
        ibd='{macro}/{micro}/{seed}/large.chr1.hapibd.candidate.ibd.gz'
    params:
        out='{macro}/{micro}/{seed}/large.chr1.ibdends',
    output:
        ibd='{macro}/{micro}/{seed}/large.chr1.ibdends.ibd.gz',
        log='{macro}/{micro}/{seed}/large.chr1.ibdends.log',
	threads:
		'{config[CHANGE][CLUSTER][LARGETHREAD]}',
	resources:
		mem_gb='{config[CHANGE][CLUSTER][LARGEMEM]}'
    shell:
        'java -Xmx{config[CHANGE][CLUSTER][LARGEMEM]}g -jar {config[CHANGE][FOLDERS][SOFTWARE]}/{config[CHANGE][PROGRAMS][IBDENDS]} gt={input.vcf} ibd={input.ibd} map={input.map} out={params.out} min-maf={config[FIXED][IBDENDS][MINMAF]} quantiles={config[FIXED][IBDENDS][QUANTILES]} nsamples={config[FIXED][IBDENDS][NSAMPLES]} err={config[FIXED][IBDENDS][ERRRATE]}'

rule format_ibdends:
	input:
		bigibd='{macro}/{micro}/{seed}/large.chr1.ibdends.ibd.gz',
	output:
		cutibd='{macro}/{micro}/{seed}/large.chr1.ibdends.modified.ibd.gz',
	shell:
		'zcat {input.bigibd} | tail -n +2 | cut -f 1-5,10,11,12 | gzip > {output.cutibd}'

rule filter_ibdends:
    input:
        ibd='{macro}/{micro}/{seed}/large.chr1.ibdends.modified.ibd.gz',
    output:
        fipass='{macro}/{micro}/{seed}/scan.chr1.ibd.gz',
    shell:
        'zcat {input.ibd} | java -jar {config[CHANGE][FOLDERS][SOFTWARE]}/{config[CHANGE][PROGRAMS][FILTER]} "I" -8 0.00 {config[FIXED][ISWEEP][SCANCUTOFF]} | gzip > {output.fipass}'

rule count_ibdends:
	input:
		filein='{macro}/{micro}/{seed}/scan.chr1.ibd.gz',
		mapin=str(config["FOLDERS"]["MACRO"]) + '/uniform.map',
	output:
		fileout='{macro}/{micro}/{seed}/scan.chr1.windowed.tsv.gz',
	script:
		'{config[CHANGE][FOLDERS][SNAKESCRIPTS]}/ibd-windowed.py'
