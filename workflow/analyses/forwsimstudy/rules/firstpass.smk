wildcard_constraints:
	SIMNAME = '\w+',

n=int(float(config['FIXED']['SAMPSIZE']))
ploidy=int(float(config['FIXED']['PLOIDY']))
maf1=float(config['FIRSTHAP']['MINMAF'])
mac1=int(ploidy*n*maf1)
maf2=float(config['SECONDHAP']['MINMAF'])
mac2=int(ploidy*n*maf2)

### ibd-ends from candidate segments ###

rule hapibd:
    input:
        vcf='{macro}/{micro}/large.vcf.gz',
        map=str(config["FOLDERS"]["MACRO"]) + '/uniform.map',
    params:
        minmac=str(mac1),
        out='{macro}/{micro}/large.hapibd.cand',
    output:
        ibd='{macro}/{micro}/large.hapibd.cand.ibd.gz',
        hbd='{macro}/{micro}/large.hapibd.cand.hbd.gz',
        log='{macro}/{micro}/large.hapibd.cand.log',
    shell:
        'java -jar {config[PATHS][HAPIBD]} gt={input.vcf} map={input.map} out={params.out} min-seed={config[FIRSTHAP][MINSEED]} min-extend={config[FIRSTHAP][MINEXT]} min-output={config[FIRSTHAP][MINOUT]} min-mac={params.minmac}'

rule ibdends:
    input:
        vcf='{macro}/{micro}/large.vcf.gz',
        map=str(config["FOLDERS"]["MACRO"]) + '/uniform.map',
        ibd='{macro}/{micro}/large.hapibd.cand.ibd.gz'
    params:
        out='{macro}/{micro}/large.ibdends.p1',
    output:
        ibd='{macro}/{micro}/large.ibdends.p1.ibd.gz',
        log='{macro}/{micro}/large.ibdends.p1.log',
    shell:
        'java -jar {config[PATHS][IBDENDS]} gt={input.vcf} ibd={input.ibd} map={input.map} out={params.out} min-maf={config[FIRSTEND][MINMAF]} quantiles={config[FIRSTEND][QUANTILES]} nsamples={config[FIRSTEND][NSAMPLES]} err={config[FIRSTEND][ERRRATE]}'

rule format_ibdends:
	input:
		bigibd='{macro}/{micro}/large.ibdends.p1.ibd.gz',
	output:
		cutibd='{macro}/{micro}/large.ibdends.ibd.gz',
	shell:
		'zcat {input.bigibd} | tail -n +2 | cut -f 1-5,10,11,12 | gzip > {output.cutibd}'

rule filter_ibdends:
    input:
        ibd='{macro}/{micro}/large.ibdends.ibd.gz',
    output:
        fipass='{macro}/{micro}/firstpass.ibdends.ibd.gz',
    shell:
        'zcat {input.ibd} | java -jar {config[PATHS][FILTER]} "I" -8 0.00 {config[ISWEEPPARAMS][IBDCUTOFF]} | gzip > {output.fipass}'

rule count_ibdends:
	input:
		filein='{macro}/{micro}/firstpass.ibdends.ibd.gz',
		mapin=str(config["FOLDERS"]["MACRO"]) + '/uniform.map',
	output:
		fileout='{macro}/{micro}/firstpass.ibdends.window.tsv.gz',
	script:
		'../../../snakescripts/ibd-rate-by-window-ibdends.py'

rule central_ibdends:
	input:
		filein='{macro}/{micro}/firstpass.ibdends.window.tsv.gz',
	output:
		fileout='{macro}/{micro}/firstpass.ibdends.central.tsv',
	script:
		'../../../snakescripts/ibd-central.py'

### hap-ibd with sequence parameters ###

rule hapibd_seq:
    input:
        vcf='{macro}/{micro}/large.vcf.gz',
        map=str(config["FOLDERS"]["MACRO"]) + '/uniform.map',
    params:
        minmac=str(mac2),
        out='{macro}/{micro}/large.hapibd',
    output:
        ibd='{macro}/{micro}/large.hapibd.ibd.gz',
        hbd='{macro}/{micro}/large.hapibd.hbd.gz',
        log='{macro}/{micro}/large.hapibd.log',
    shell:
        'java -jar {config[PATHS][HAPIBD]} gt={input.vcf} map={input.map} out={params.out} min-seed={config[SECONDHAP][MINSEED]} min-extend={config[SECONDHAP][MINEXT]} min-output={config[ISWEEPPARAMS][IBDCUTOFF]} min-mac={params.minmac}'

rule filter_hapibd:
	input:
		ibd='{macro}/{micro}/large.hapibd.ibd.gz',
	output:
		fhpass='{macro}/{micro}/firstpass.hapibd.ibd.gz',
	shell:
		'zcat {input.ibd} | java -jar {config[PATHS][FILTER]} "I" -8 0.00 {config[ISWEEPPARAMS][IBDCUTOFF]} | gzip > {output.fhpass}'

rule count_hapibd:
	input:
		filein='{macro}/{micro}/firstpass.hapibd.ibd.gz',
		mapin=str(config["FOLDERS"]["MACRO"]) + '/uniform.map',
	output:
		fileout='{macro}/{micro}/firstpass.hapibd.window.tsv.gz',
	script:
		'../../../snakescripts/ibd-rate-by-window-hapibd.py'

rule central_hapibd:
	input:
		filein='{macro}/{micro}/firstpass.hapibd.window.tsv.gz',
	output:
		fileout='{macro}/{micro}/firstpass.hapibd.central.tsv',
	script:
		'../../../snakescripts/ibd-central.py'
