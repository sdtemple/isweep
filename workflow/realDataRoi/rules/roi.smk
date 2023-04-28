# there are no wildcards here

n=int(float(config['CHANGE']['ISWEEP']['SAMPSIZE']))
ploidy=int(float(config['FIXED']['HAPIBD']['PLOIDY']))
maf1=float(config['FIXED']['HAPIBD']['MINMAF'])
mac1=int(ploidy*n*maf1)
largethreads=int(float(str(config['CHANGE']['CLUSTER']['LARGETHREAD'])))

# do this before running snakemake
# mkdir -p {config[CHANGE][FOLDERS][STUDY]}/{config[CHANGE][FOLDERS][ROI]}

rule copy_vcf:
    input:
        scandone='{config[CHANGE][FOLDERS][STUDY]}/excess.regions.tsv',
    output:
        copydone='{config[CHANGE][FOLDERS][STUDY]}/ROI.tsv',
    resources:
        mem_gb='{config[CHANGE][CLUSTER][LARGEMEM]}',
    script:
        '{config[CHANGE][FOLDERS][TERMINALSCRIPTS]}/copy-vcf.py {config[CHANGE][FOLDERS][STUDY]} {config[CHANGE][FOLDERS][ROI]}'

rule subset_vcf: # focus vcf on region of interest
    input:
        vcfin='{cohort}/roi{roi}/chr{chr}/center{center}/left{left}/right{right}/chr{chr}.vcf.gz',
    output:
        subvcf='{cohort}/roi{roi}/chr{chr}/center{center}/left{left}/right{right}/chr{chr}.focused.vcf.gz',
    resources:
        mem_gb='{config[CHANGE][CLUSTER][LARGEMEM]}',
    shell: # if chromosome is huge (greater than 10000 Mb), may need to modify the third pipe
        'bash {config[CHANGE][FOLDERS][TERMINALSCRIPTS]}/roi-vcf.sh {cohort}/{roi}/{center}/{left}/{right} {input.vcfin} {output.vcfout} {left} {right} {chr} {config[FIXED][ISWEEP][MINAF]}; rm {input.vcfin}'

rule hapibd: # segments from hap-ibd.jar
    input:
        vcf='{cohort}/roi{roi}/chr{chr}/center{center}/left{left}/right{right}/chr{chr}.focused.vcf.gz',
        map='{cohort}/maps/chr{chr}.map',
    params:
        minmac=str(mac1),
        out='{cohort}/roi{roi}/chr{chr}/center{center}/left{left}/right{right}/',
    output:
        ibd='{cohort}/roi{roi}/chr{chr}/center{center}/left{left}/right{right}/chr{chr}.ibd.gz',
        hbd='{cohort}/roi{roi}/chr{chr}/center{center}/left{left}/right{right}/chr{chr}.hbd.gz',
        log='{cohort}/roi{roi}/chr{chr}/center{center}/left{left}/right{right}/chr{chr}.log',
    threads:
        largethreads,
    resources:
        mem_gb='{config[CHANGE][CLUSTER][LARGEMEM]}',
    shell:
        'java -Xmx{config[CHANGE][CLUSTER][LARGEMEM]}g -jar {config[CHANGE][FOLDERS][SOFTWARE]}/{config[CHANGE][PROGRAMS][HAPIBD]} gt={input.vcf} map={input.map} out={params.out} min-seed={config[FIXED][HAPIBD][MINSEED]} min-extend={config[FIXED][HAPIBD][MINEXT]} min-output={config[FIXED][HAPIBD][MINOUT]} min-mac={params.minmac}'

rule filter_hapibd: # applying focus
    input:
        ibd='{cohort}/roi{roi}/chr{chr}/center{center}/left{left}/right{right}/chr{chr}.ibd.gz',
    output:
        fipass='{cohort}/roi{roi}/chr{chr}/center{center}/left{left}/right{right}/focus.ibd.gz',
    resources:
        mem_gb='{config[CHANGE][CLUSTER][LARGEMEM]}',
    shell: # if chromosome is huge (greater than 10000 Mb), may need to modify the third pipe
        'zcat {input.ibd} | java -Xmx{config[CHANGE][CLUSTER][LARGEMEM]}g -jar {config[CHANGE][FOLDERS][SOFTWARE]}/{config[CHANGE][PROGRAMS][FILTER]} "I" 6 0.00 {center} | java -Xmx{config[CHANGE][CLUSTER][LARGEMEM]}g -jar {config[CHANGE][FOLDERS][SOFTWARE]}/{config[CHANGE][PROGRAMS][FILTER]} "I" 7 {center} 10000000000 | gzip > {output.fipass}; rm {input.ibd}'

rule filter_mom:
    input:
        '{cohort}/ibdsegs/ibdends/modified/mom/chr{chr}.ibd.gz',
    output:
        fipass='{cohort}/{roi}/chr{chr}/center{center}/left{left}/right{right}/mom.ibd.gz',
    resources:
        mem_gb='{config[CHANGE][CLUSTER][LARGEMEM]}',
    shell: # if chromosome is huge (greater than 10000 Mb), may need to modify the third pipe
        'zcat {input.ibd} | java -Xmx{config[CHANGE][CLUSTER][LARGEMEM]}g -jar {config[CHANGE][FOLDERS][SOFTWARE]}/{config[CHANGE][PROGRAMS][FILTER]} "I" 6 0.00 {center} | java -Xmx{config[CHANGE][CLUSTER][LARGEMEM]}g -jar {config[CHANGE][FOLDERS][SOFTWARE]}/{config[CHANGE][PROGRAMS][FILTER]} "I" 7 {center} 10000000000 | gzip > {output.fipass}'

rule done_filtering:
    input:
        '{cohort}/{roi}/chr{chr}/center{center}/left{left}/right{right}/mom.ibd.gz',
        '{cohort}/roi{roi}/chr{chr}/center{center}/left{left}/right{right}/focus.ibd.gz',
    output:
        '{cohort}/roi{roi}/chr{chr}/center{center}/left{left}/right{right}/touched',
    shell:
        'touch {cohort}/roi{roi}/chr{chr}/center{center}/left{left}/right{right}/touched',
