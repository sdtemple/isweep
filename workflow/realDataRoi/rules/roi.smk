# there are no wildcards here

n=int(float(config['CHANGE']['ISWEEP']['SAMPSIZE']))
ploidy=int(float(config['FIXED']['HAPIBD']['PLOIDY']))
maf1=float(config['FIXED']['HAPIBD']['MINMAF'])
mac1=int(ploidy*n*maf1)
smallthreads=int(float(str(config['CHANGE']['CLUSTER']['SMALLTHREAD'])))
largethreads=int(float(str(config['CHANGE']['CLUSTER']['LARGETHREAD'])))

# do this before running snakemake
# mkdir -p {config[CHANGE][FOLDERS][STUDY]}/{config[CHANGE][FOLDERS][ROI]}

rule subset_vcf: # focus vcf on region of interest
	input:
		vcfin='{cohort}/vcfs/chr{chr}.vcf.gz',
	output:
		subvcf='{cohort}/{config[CHANGE][FOLDERS][ROI]}/chr{chr}.vcf.gz',
	shell:
		'bash {config[CHANGE][FOLDERS][TERMINALSCRIPTS]}/roi-vcf.sh {cohort}/{roi} {input.vcfin} {output.vcfout} {config[CHANGE][ROI][LEFT]} {config[CHANGE][ROI][RIGHT]} {chr} {config[FIXED][ISWEEP][MINAF]}'

rule hapibd: # segments from hap-ibd.jar
    input:
        vcf='{cohort}/{roi}/chr{chr}.vcf.gz',
        map='{cohort}/maps/chr{chr}.map',
    params:
        minmac=str(mac1),
        out='{cohort}/{roi}/chr{chr}',
    output:
        ibd='{cohort}/{roi}/chr{chr}.ibd.gz',
        hbd='{cohort}/{roi}/chr{chr}.hbd.gz',
        log='{cohort}/{roi}/chr{chr}.log',
	threads:
		largethreads,
	resources:
		mem_gb='{config[CHANGE][PROGRAMS][LARGEMEM]}',
    shell:
        'java -Xmx{config[CHANGE][PROGRAMS][LARGEMEM]}g -jar {config[CHANGE][FOLDERS][SOFTWARE]}/{config[CHANGE][PROGRAMS][HAPIBD]} gt={input.vcf} map={input.map} out={params.out} min-seed={config[FIXED][HAPIBD][MINSEED]} min-extend={config[FIXED][HAPIBD][MINEXT]} min-output={config[FIXED][HAPIBD][MINOUT]} min-mac={params.minmac}'

rule filter_hapibd: # applying focus
    input:
        ibd='{cohort}/{roi}/chr{config[CHANGE][ROI][CHRNUM]}.ibd.gz',
    output:
        fipass='{cohort}/{roi}/focus.ibd.gz',
    shell: # if chromosome is huge (greater than 10000 Mb), may need to modify the third pipe
        'zcat {input.ibd} | java -jar {config[CHANGE][FOLDERS][SOFTWARE]}/{config[CHANGE][PROGRAMS][FILTER]} "I" 6 0.00 {config[CHANGE][ROI][CENTER]} | java -jar {config[CHANGE][FOLDERS][SOFTWARE]}/{config[CHANGE][PROGRAMS][FILTER]} "I" 7 {config[CHANGE][ROI][CENTER]} 10000000000 | gzip > {output.fipass}'

rule filter_mom:
	input:
		'{cohort}/ibdsegs/ibdends/modified/mom/chr{config[CHANGE][ROI][CHRNUM]}.ibd.gz',
    output:
        fipass='{cohort}/{config[CHANGE][FOLDERS][ROI]}/mom.ibd.gz',
    shell: # if chromosome is huge (greater than 10000 Mb), may need to modify the third pipe
        'zcat {input.ibd} | java -jar {config[CHANGE][FOLDERS][SOFTWARE]}/{config[CHANGE][PROGRAMS][FILTER]} "I" 6 0.00 {config[CHANGE][ROI][CENTER]} | java -jar {config[CHANGE][FOLDERS][SOFTWARE]}/{config[CHANGE][PROGRAMS][FILTER]} "I" 7 {config[CHANGE][ROI][CENTER]} 10000000000 | gzip > {output.fipass}'
