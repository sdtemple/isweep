# preparing regions of interest with browning lab software
n=int(float(config['CHANGE']['ISWEEP']['SAMPSIZE']))
ploidy=int(float(config['FIXED']['HAPIBD']['PLOIDY']))
maf1=float(config['FIXED']['HAPIBD']['MINMAF'])
mac1=int(ploidy*n*maf1)
maf2=float(config['FIXED']['ISWEEP']['MINAF'])
maf=min(maf1,maf2)
macro=str(config['CHANGE']['FOLDERS']['STUDY'])

# make an excludesamples.txt file (just in case)
rule touch_exclude:
    input:
        excludesamples=macro+"/subsample.txt",
    output:
        excludesamples=macro+"/excludesamples.txt",
    shell:
        'touch {output.excludesamples}'

# subset vcf to region of interest
rule subset_vcf: # focus vcf on region of interest
    input:
        folder='{cohort}/roi{roi}/chr{chr}/center{center}/left{left}/right{right}/ready',
        vcf='{config[CHANGE][EXISTING][VCFS]}/{config[CHANGE][EXISTING][VCFPRE]}{chr}{config[CHANGE][EXISTING][VCFSUF]}',
        # vcfin='{cohort}/roi{roi}/chr{chr}/center{center}/left{left}/right{right}/chr{chr}.vcf.gz',
    output:
        subvcf='{cohort}/roi{roi}/chr{chr}/center{center}/left{left}/right{right}/chr{chr}.focused.vcf.gz',
    params:
        qmaf=maf,
        subsample=macro+"/subsample.txt",
    shell: # if chromosome is huge (greater than 10000 Mb), may need to modify the third pipe
        'bash {config[CHANGE][FOLDERS][TERMINALSCRIPTS]}/roi-vcf.sh {wildcards.cohort}/roi{wildcards.roi}/center{wildcards.center}/left{wildcards.left}/right{wildcards.right} {input.vcfin} {output.subvcf} {wildcards.left} {wildcards.right} {params.subsample} {params.qmaf}'

# perform hap-ibd with sequence parameter settings
rule hapibd: # segments from hap-ibd.jar
    input:
        vcf='{cohort}/roi{roi}/chr{chr}/center{center}/left{left}/right{right}/chr{chr}.focused.vcf.gz',
        map='{cohort}/maps/chr{chr}.map',
        excludesamples=macro+"/excludesamples.txt",
    params:
        minmac=str(mac1),
        out='{cohort}/roi{roi}/chr{chr}/center{center}/left{left}/right{right}/',
    output:
        ibd='{cohort}/roi{roi}/chr{chr}/center{center}/left{left}/right{right}/chr{chr}.ibd.gz',
        hbd='{cohort}/roi{roi}/chr{chr}/center{center}/left{left}/right{right}/chr{chr}.hbd.gz',
        log='{cohort}/roi{roi}/chr{chr}/center{center}/left{left}/right{right}/chr{chr}.log',
    shell:
        'java -Xmx{config[CHANGE][ISWEEP][XMXMEM]}g -jar {config[CHANGE][FOLDERS][SOFTWARE]}/{config[CHANGE][PROGRAMS][HAPIBD]} gt={input.vcf} map={input.map} out={params.out} min-seed={config[FIXED][HAPIBD][MINSEED]} min-extend={config[FIXED][HAPIBD][MINEXT]} min-output={config[FIXED][HAPIBD][MINOUT]} min-mac={params.minmac} excludesamples={input.excludesamples}'

# filter outputs of hap-ibd to be about focal location
rule filter_hapibd: # applying focus
    input:
        ibd='{cohort}/roi{roi}/chr{chr}/center{center}/left{left}/right{right}/chr{chr}.ibd.gz',
    output:
        fipass='{cohort}/roi{roi}/chr{chr}/center{center}/left{left}/right{right}/focus.ibd.gz',
    shell: # if chromosome is huge (greater than 10000 Mb), may need to modify the third pipe
        'zcat {input.ibd} | java -Xmx{config[CHANGE][ISWEEP][XMXMEM]}g -jar {config[CHANGE][FOLDERS][SOFTWARE]}/{config[CHANGE][PROGRAMS][FILTER]} "I" 6 0.00 {wildcards.center} | java -Xmx{config[CHANGE][ISWEEP][XMXMEM]}g -jar {config[CHANGE][FOLDERS][SOFTWARE]}/{config[CHANGE][PROGRAMS][FILTER]} "I" 7 {wildcards.center} 10000000000 | gzip > {output.fipass}; rm {input.ibd}'

# filter method of moments to be about focal location
rule filter_mom:
    input:
        ibd='{cohort}/ibdsegs/ibdends/modified/mom/chr{chr}.ibd.gz',
    output:
        fipass='{cohort}/{roi}/chr{chr}/center{center}/left{left}/right{right}/mom.ibd.gz',
    resources:
        mem_gb=50
    shell: # if chromosome is huge (greater than 10000 Mb), may need to modify the third pipe
        'zcat {input.ibd} | java -Xmx{config[CHANGE][ISWEEP][XMXMEM]}g -jar {config[CHANGE][FOLDERS][SOFTWARE]}/{config[CHANGE][PROGRAMS][FILTER]} "I" 6 0.00 {wildcards.center} | java -Xmx{config[CHANGE][ISWEEP][XMXMEM]}g -jar {config[CHANGE][FOLDERS][SOFTWARE]}/{config[CHANGE][PROGRAMS][FILTER]} "I" 7 {wildcards.center} 10000000000 | gzip > {output.fipass}'
