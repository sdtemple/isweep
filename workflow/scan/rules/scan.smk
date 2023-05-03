# conduct ibd calls and scan for isweep

n=int(float(config['CHANGE']['ISWEEP']['SAMPSIZE']))
ploidy=int(float(config['FIXED']['CANDHAPIBD']['PLOIDY']))
maf1=float(config['FIXED']['CANDHAPIBD']['MINMAF'])
mac1=int(ploidy*n*maf1)
macro=str(config['CHANGE']['FOLDERS']['STUDY'])
low=int(float(str(config['CHANGE']['ISWEEP']['CHRLOW'])))
high=int(float(str(config['CHANGE']['ISWEEP']['CHRHIGH'])))
vcffolder=str(config['CHANGE']['EXISTING']['VCFS'])
vcfpre=str(config['CHANGE']['EXISTING']['VCFPRE'])
vcfsuf=str(config['CHANGE']['EXISTING']['VCFSUF'])
concatmap=macro+'/maps/chr'+str(low)+'-'+str(high)+'.map'

rule hapibd: # candidate segments from hap-ibd.jar
    input:
        vcf=vcffolder + '/' + vcfpre + '{num}' + vcfsuf,
        map='{cohort}/maps/chr{num}.map',
        subsample='{cohort}/excludesamples.txt',
    params:
        minmac=str(mac1),
        out='{cohort}/ibdsegs/hapibd/chr{num}',
    output:
        ibd='{cohort}/ibdsegs/hapibd/chr{num}.ibd.gz',
        hbd='{cohort}/ibdsegs/hapibd/chr{num}.hbd.gz',
        log='{cohort}/ibdsegs/hapibd/chr{num}.log',
    resources:
        mem_gb=100
    shell:
        'java -Xmx{config[CHANGE][ISWEEP][XMXMEM]}g -jar {config[CHANGE][FOLDERS][SOFTWARE]}/{config[CHANGE][PROGRAMS][HAPIBD]} gt={input.vcf} map={input.map} out={params.out} min-seed={config[FIXED][CANDHAPIBD][MINSEED]} min-extend={config[FIXED][CANDHAPIBD][MINEXT]} min-output={config[FIXED][CANDHAPIBD][MINOUT]} min-mac={params.minmac} excludesamples={input.subsample}'

rule ibdends: # ibd-ends.jar
    input:
        vcf=vcffolder + '/' + vcfpre + '{num}' + vcfsuf,
        map='{cohort}/maps/chr{num}.map',
        ibd='{cohort}/ibdsegs/hapibd/chr{num}.ibd.gz',
        subsample='{cohort}/excludesamples.txt',
    params:
        out='{cohort}/ibdsegs/ibdends/chr{num}',
    output:
        ibd='{cohort}/ibdsegs/ibdends/chr{num}.ibd.gz',
        log='{cohort}/ibdsegs/ibdends/chr{num}.log',
    resources:
        mem_gb=100
    shell:
        'java -Xmx{config[CHANGE][ISWEEP][XMXMEM]}g -jar {config[CHANGE][FOLDERS][SOFTWARE]}/{config[CHANGE][PROGRAMS][IBDENDS]} gt={input.vcf} ibd={input.ibd} map={input.map} out={params.out} min-maf={config[FIXED][IBDENDS][MINMAF]} quantiles={config[FIXED][IBDENDS][QUANTILES]} nsamples={config[FIXED][IBDENDS][NSAMPLES]} err={config[FIXED][IBDENDS][ERRRATE]} excludesamples={input.subsample}'

rule format_ibdends: # reformatting
    input:
        bigibd='{cohort}/ibdsegs/ibdends/chr{num}.ibd.gz',
    output:
        cutibd='{cohort}/ibdsegs/ibdends/modified/chr{num}.ibd.gz',
    resources:
        mem_gb=10
    shell:
        'zcat {input.bigibd} | tail -n +2 | cut -f 1-5,10,11,12 | gzip > {output.cutibd}'

rule filter_ibdends_scan: # applying cutoffs
    input:
        ibd='{cohort}/ibdsegs/ibdends/modified/chr{num}.ibd.gz',
    output:
        fipass='{cohort}/ibdsegs/ibdends/modified/scan/chr{num}.ibd.gz',
    resources:
        mem_gb=10
    shell:
        'zcat {input.ibd} | java -Xmx{config[CHANGE][ISWEEP][XMXMEM]}g -jar {config[CHANGE][FOLDERS][SOFTWARE]}/{config[CHANGE][PROGRAMS][FILTER]} "I" -8 0.00 {config[FIXED][ISWEEP][SCANCUTOFF]} | gzip > {output.fipass}'

rule count_ibdends_scan: # computing counts over windows
    input:
        filein='{cohort}/ibdsegs/ibdends/modified/scan/chr{num}.ibd.gz',
        mapin='{cohort}/maps/chr{num}.map',
    output:
        fileout='{cohort}/ibdsegs/ibdends/modified/scan/chr{num}.ibd.windowed.tsv.gz',
    resources:
        mem_gb=10
    shell:
        'python {config[CHANGE][FOLDERS][TERMINALSCRIPTS]}/ibd-windowed.py {input.filein} {output.fileout} {input.mapin} {config[FIXED][ISWEEP][BY]} {config[FIXED][ISWEEP][GENOMEEND]} {config[FIXED][ISWEEP][TELOCUT]}'

rule filter_ibdends_mom: # applying cutoffs
    input:
        ibd='{cohort}/ibdsegs/ibdends/modified/chr{num}.ibd.gz',
    output:
        fipass='{cohort}/ibdsegs/ibdends/modified/mom/chr{num}.ibd.gz',
    resources:
        mem_gb=10
    shell:
        'zcat {input.ibd} | java -Xmx{config[CHANGE][ISWEEP][XMXMEM]}g -jar {config[CHANGE][FOLDERS][SOFTWARE]}/{config[CHANGE][PROGRAMS][FILTER]} "I" -8 0.00 {config[FIXED][ISWEEP][MOMCUTOFF]} | gzip > {output.fipass}'

rule count_ibdends_mom: # computing counts over windows
    input:
        filein='{cohort}/ibdsegs/ibdends/modified/mom/chr{num}.ibd.gz',
        mapin='{cohort}/maps/chr{num}.map',
    output:
        fileout='{cohort}/ibdsegs/ibdends/modified/mom/chr{num}.ibd.windowed.tsv.gz',
    resources:
        mem_gb=10
    shell:
        'python {config[CHANGE][FOLDERS][TERMINALSCRIPTS]}/ibd-windowed.py {input.filein} {output.fileout} {input.mapin} {config[FIXED][ISWEEP][BY]} {config[FIXED][ISWEEP][GENOMEEND]} {config[FIXED][ISWEEP][TELOCUT]}'

rule scan: # conduct a manhattan scan
    input:
        [macro+'/ibdsegs/ibdends/modified/scan/chr'+str(i)+'.ibd.windowed.tsv.gz' for i in range(low,high+1)],
        [macro+'/ibdsegs/ibdends/modified/mom/chr'+str(i)+'.ibd.windowed.tsv.gz' for i in range(low,high+1)],
    output:
        macro+'/excess.ibd.tsv',
    resources:
        mem_gb=5
    shell:
        'python {config[CHANGE][FOLDERS][TERMINALSCRIPTS]}/scan-isweep.py {config[CHANGE][FOLDERS][STUDY]} {config[CHANGE][ISWEEP][CHRLOW]} {config[CHANGE][ISWEEP][CHRHIGH]} {config[FIXED][ISWEEP][SCANSIGMA]} {config[FIXED][ISWEEP][TELOSIGMA]}'
