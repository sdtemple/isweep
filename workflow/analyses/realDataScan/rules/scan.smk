wildcard_constraints:
	SIMNAME = '\w+',

n=int(float(config['CHANGE']['ISWEEP']['SAMPSIZE']))
ploidy=int(float(config['FIXED']['CANDHAPIBD']['PLOIDY']))
maf1=float(config['FIXED']['CANDHAPIBD']['MINMAF'])
mac1=int(ploidy*n*maf1)
macro=str(config['CHANGE']['FOLDERS']['STUDY'])
low=int(float(str(config['CHANGE']['ISWEEP']['CHRLOW'])))
high=int(float(str(config['CHANGE']['ISWEEP']['CHRHIGH'])))

rule hapibd: # candidate segments from hap-ibd.jar
    input:
        vcf='{config[CHANGE][FOLDERS][STUDY]}/vcfs/chr{num}.vcf.gz',
        map='{config[CHANGE][FOLDERS][STUDY]}/maps/chr{num}.map',
    params:
        minmac=str(mac1),
        out='{config[CHANGE][FOLDERS][STUDY]}/ibdsegs/hapibd/chr${num}',
    output:
        ibd='{config[CHANGE][FOLDERS][STUDY]}/ibdsegs/hapibd/chr${num}.ibd.gz',
        hbd='{config[CHANGE][FOLDERS][STUDY]}/ibdsegs/hapibd/chr${num}.hbd.gz',
        log='{config[CHANGE][FOLDERS][STUDY]}/ibdsegs/hapibd/chr${num}.log',
    shell:
        'java -jar {config[CHANGE][FOLDERS][SOFTWARE]}/{config[CHANGE][PROGRAMS][HAPIBD]} gt={input.vcf} map={input.map} out={params.out} min-seed={config[FIXED][CANDHAPIBD][MINSEED]} min-extend={config[FIXED][CANDHAPIBD][MINEXT]} min-output={config[FIXED][CANDHAPIBD][MINOUT]} min-mac={params.minmac}'

rule ibdends: # ibd-ends.jar
    input:
        vcf='{config[CHANGE][FOLDERS][STUDY]}/vcfs/chr{num}.vcf.gz',
        map='{config[CHANGE][FOLDERS][STUDY]}/maps/chr{num}.map',
        ibd='{config[CHANGE][FOLDERS][STUDY]}/ibdsegs/hapibd/chr${num}.ibd.gz'
    params:
        out='{config[CHANGE][FOLDERS][STUDY]}/ibdsegs/ibdends/chr${num}',
    output:
        ibd='{config[CHANGE][FOLDERS][STUDY]}/ibdsegs/ibdends/chr${num}.ibd.gz',
        hbd='{config[CHANGE][FOLDERS][STUDY]}/ibdsegs/ibdends/chr${num}.hbd.gz',
        log='{config[CHANGE][FOLDERS][STUDY]}/ibdsegs/ibdends/chr${num}.log',
    shell:
        'java -jar {config[CHANGE][FOLDERS][SOFTWARE]}/{config[CHANGE][PROGRAMS][IBDENDS]} gt={input.vcf} ibd={input.ibd} map={input.map} out={params.out} min-maf={config[FIXED][IBDENDS][MINMAF]} quantiles={config[FIXED][IBDENDS][QUANTILES]} nsamples={config[FIXED][IBDENDS][NSAMPLES]} err={config[FIXED][IBDENDS][ERRRATE]}'

rule format_ibdends: # reformatting
	input:
		bigibd='{config[CHANGE][FOLDERS][STUDY]}/ibdsegs/ibdends/chr${num}.ibd.gz',
	output:
		cutibd='{config[CHANGE][FOLDERS][STUDY]}/ibdsegs/ibdends/modified/chr${num}.ibd.gz',
	shell:
		'zcat {input.bigibd} | tail -n +2 | cut -f 1-5,10,11,12 | gzip > {output.cutibd}'

rule filter_ibdends_scan: # applying cutoffs
    input:
        ibd='{config[CHANGE][FOLDERS][STUDY]}/ibdsegs/ibdends/modified/chr${num}.ibd.gz',
    output:
        fipass='{config[CHANGE][FOLDERS][STUDY]}/ibdsegs/ibdends/modified/scan/chr${num}.ibd.gz',
    shell:
        'zcat {input.ibd} | java -jar {config[CHANGE][FOLDERS][SOFTWARE]}/{config[CHANGE][PROGRAMS][FILTER]} "I" -8 0.00 {config[FIXED][ISWEEP][SCANCUTOFF]} | gzip > {output.fipass}'

rule count_ibdends_scan: # computing counts over windows
	input:
		filein='{config[CHANGE][FOLDERS][STUDY]}/ibdsegs/ibdends/modified/scan/chr${num}.ibd.gz',
		mapin='{config[CHANGE][FOLDERS][STUDY]}/maps/chr{num}.map',
	output:
		fileout='{config[CHANGE][FOLDERS][STUDY]}/ibdsegs/ibdends/modified/scan/chr${num}.windowed.tsv.gz',
	script:
		'{config[CHANGE][FOLDERS][SNAKESCRIPTS]}/ibd-windowed.py'

rule filter_ibdends_mom: # applying cutoffs
    input:
        ibd='{config[CHANGE][FOLDERS][STUDY]}/ibdsegs/ibdends/modified/chr${num}.ibd.gz',
    output:
        fipass='{config[CHANGE][FOLDERS][STUDY]}/ibdsegs/ibdends/modified/mom/chr${num}.ibd.gz',
    shell:
        'zcat {input.ibd} | java -jar {config[CHANGE][FOLDERS][SOFTWARE]}/{config[CHANGE][PROGRAMS][FILTER]} "I" -8 0.00 {config[FIXED][ISWEEP][MOMCUTOFF]} | gzip > {output.fipass}'

rule count_ibdends_mom: # computing counts over windows
	input:
		filein='{config[CHANGE][FOLDERS][STUDY]}/ibdsegs/ibdends/modified/mom/chr${num}.ibd.gz',
		mapin='{config[CHANGE][FOLDERS][STUDY]}/maps/chr{num}.map',
	output:
		fileout='{config[CHANGE][FOLDERS][STUDY]}/ibdsegs/ibdends/modified/mom/chr${num}.windowed.tsv.gz',
	script:
		'{config[CHANGE][FOLDERS][SNAKESCRIPTS]}/ibd-windowed.py'

rule scan: # conduct a manhattan scan
	input:
		[macro+'/ibdsegs/ibdends/modified/scan/chr'+i+'.ibd.windowed.tsv.gz' for i in range(low,high+1)],
		[macro+'/ibdsegs/ibdends/modified/mom/chr'+i+'.ibd.windowed.tsv.gz' for i in range(low,high+1)],
	output:
		macro+'/excess.ibd.tsv',
	script:
		'{config[CHANGE][FOLDERS][SNAKESCRIPTS]}/scan-isweep.py'

rule concat_maps: # concatenate maps for ibdne
	input:
		[macro+'/maps/ibdends/modified/scan/chr'+i+'.ibd.gz' for i in range(low,high+1)],
	output:
		allmap='{config[CHANGE][FOLDERS][STUDY]}/maps/chr{config[CHANGE][ISWEEP][CHRLOW]}-{config[CHANGE][ISWEEP][CHRHIGH]}.map',
	shell:
		'for j in $(seq {config[CHANGE][ISWEEP][CHRLOW]} 1 {config[CHANGE][ISWEEP][CHRHIGH]}); do cat {config[CHANGE][FOLDERS][STUDY]}/maps/chr${i}.map >> {output.allmap}; done;'

rule ibdne: # run ibdne
	input:
		allmap='{config[CHANGE][FOLDERS][STUDY]}/maps/chr{config[CHANGE][ISWEEP][CHRLOW]}-{config[CHANGE][ISWEEP][CHRHIGH]}.map',
		ibdfiles=[macro+'/ibdsegs/ibdends/modified/scan/chr'+i+'.ibd.gz' for i in range(low,high+1)],
	output:
		'{config[CHANGE][FOLDERS][STUDY]}/ne/inferredNe.ne',
	params:
		outpre='{config[CHANGE][FOLDERS][STUDY]}/ne/inferredNe',
	shell:
		'zcat {config[CHANGE][FOLDERS][STUDY]}/ibdsegs/ibdends/modified/scan/chr*.ibd.gz | java -jar {config[CHANGE][FOLDERS][SOFTWARE]}/{config[CHANGE][PROGRAMS][IBDNE]} map={input.allmap} mincm={config[FIXED][ISWEEP][SCANCUTOFF]} out={params.outpre}'
