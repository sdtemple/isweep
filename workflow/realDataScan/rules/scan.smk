n=int(float(config['CHANGE']['ISWEEP']['SAMPSIZE']))
ploidy=int(float(config['FIXED']['CANDHAPIBD']['PLOIDY']))
maf1=float(config['FIXED']['CANDHAPIBD']['MINMAF'])
mac1=int(ploidy*n*maf1)
macro=str(config['CHANGE']['FOLDERS']['STUDY'])
low=int(float(str(config['CHANGE']['ISWEEP']['CHRLOW'])))
high=int(float(str(config['CHANGE']['ISWEEP']['CHRHIGH'])))

rule hapibd: # candidate segments from hap-ibd.jar
    input:
        vcf='{cohort}/vcfs/chr{num}.vcf.gz',
        map='{cohort}/maps/chr{num}.map',
    params:
        minmac=str(mac1),
        out='{cohort}/ibdsegs/hapibd/chr{num}',
    output:
        ibd='{cohort}/ibdsegs/hapibd/chr{num}.ibd.gz',
        hbd='{cohort}/ibdsegs/hapibd/chr{num}.hbd.gz',
        log='{cohort}/ibdsegs/hapibd/chr{num}.log',
    threads:
        int(float('{config[CHANGE][CLUSTER][SMALLTHREAD]}')),
    resources:
        mem_gb='{config[CHANGE][PROGRAMS][SMALLMEM]}',
    shell:
        'java -Xmx{config[CHANGE][PROGRAMS][SMALLMEM]}g -jar {config[CHANGE][FOLDERS][SOFTWARE]}/{config[CHANGE][PROGRAMS][HAPIBD]} gt={input.vcf} map={input.map} out={params.out} min-seed={config[FIXED][CANDHAPIBD][MINSEED]} min-extend={config[FIXED][CANDHAPIBD][MINEXT]} min-output={config[FIXED][CANDHAPIBD][MINOUT]} min-mac={params.minmac}'

rule ibdends: # ibd-ends.jar
    input:
        vcf='{cohort}/vcfs/chr{num}.vcf.gz',
        map='{cohort}/maps/chr{num}.map',
        ibd='{cohort}/ibdsegs/hapibd/chr{num}.ibd.gz'
    params:
        out='{cohort}/ibdsegs/ibdends/chr{num}',
    output:
        ibd='{cohort}/ibdsegs/ibdends/chr{num}.ibd.gz',
        log='{cohort}/ibdsegs/ibdends/chr{num}.log',
    threads:
        int(float('{config[CHANGE][CLUSTER][LARGETHREAD]}')),
    resources:
        mem_gb='{config[CHANGE][PROGRAMS][LARGEMEM]}',
    shell:
        'java -Xmx{config[CHANGE][PROGRAMS][LARGEMEM]}g -jar {config[CHANGE][FOLDERS][SOFTWARE]}/{config[CHANGE][PROGRAMS][IBDENDS]} gt={input.vcf} ibd={input.ibd} map={input.map} out={params.out} min-maf={config[FIXED][IBDENDS][MINMAF]} quantiles={config[FIXED][IBDENDS][QUANTILES]} nsamples={config[FIXED][IBDENDS][NSAMPLES]} err={config[FIXED][IBDENDS][ERRRATE]}'

rule format_ibdends: # reformatting
	input:
		bigibd='{cohort}/ibdsegs/ibdends/chr{num}.ibd.gz',
	output:
		cutibd='{cohort}/ibdsegs/ibdends/modified/chr{num}.ibd.gz',
	shell:
		'zcat {input.bigibd} | tail -n +2 | cut -f 1-5,10,11,12 | gzip > {output.cutibd}'

rule filter_ibdends_scan: # applying cutoffs
    input:
        ibd='{cohort}/ibdsegs/ibdends/modified/chr{num}.ibd.gz',
    output:
        fipass='{cohort}/ibdsegs/ibdends/modified/scan/chr{num}.ibd.gz',
    shell:
        'zcat {input.ibd} | java -jar {config[CHANGE][FOLDERS][SOFTWARE]}/{config[CHANGE][PROGRAMS][FILTER]} "I" -8 0.00 {config[FIXED][ISWEEP][SCANCUTOFF]} | gzip > {output.fipass}'

rule count_ibdends_scan: # computing counts over windows
	input:
		filein='{cohort}/ibdsegs/ibdends/modified/scan/chr{num}.ibd.gz',
		mapin='{cohort}/maps/chr{num}.map',
	output:
		fileout='{cohort}/ibdsegs/ibdends/modified/scan/chr{num}.ibd.windowed.tsv.gz',
	script:
		'{config[CHANGE][FOLDERS][SNAKESCRIPTS]}/ibd-windowed.py'

rule filter_ibdends_mom: # applying cutoffs
    input:
        ibd='{cohort}/ibdsegs/ibdends/modified/chr{num}.ibd.gz',
    output:
        fipass='{cohort}/ibdsegs/ibdends/modified/mom/chr{num}.ibd.gz',
    shell:
        'zcat {input.ibd} | java -jar {config[CHANGE][FOLDERS][SOFTWARE]}/{config[CHANGE][PROGRAMS][FILTER]} "I" -8 0.00 {config[FIXED][ISWEEP][MOMCUTOFF]} | gzip > {output.fipass}'

rule count_ibdends_mom: # computing counts over windows
	input:
		filein='{cohort}/ibdsegs/ibdends/modified/mom/chr{num}.ibd.gz',
		mapin='{cohort}/maps/chr{num}.map',
	output:
		fileout='{cohort}/ibdsegs/ibdends/modified/mom/chr{num}.ibd.windowed.tsv.gz',
	script:
		'{config[CHANGE][FOLDERS][SNAKESCRIPTS]}/ibd-windowed.py'

rule scan: # conduct a manhattan scan
	input:
		[macro+'/ibdsegs/ibdends/modified/scan/chr'+str(i)+'.ibd.windowed.tsv.gz' for i in range(low,high+1)],
		[macro+'/ibdsegs/ibdends/modified/mom/chr'+str(i)+'.ibd.windowed.tsv.gz' for i in range(low,high+1)],
	output:
		macro+'/excess.ibd.tsv',
	script:
		'{config[CHANGE][FOLDERS][SNAKESCRIPTS]}/scan-isweep.py'
