# there are no wildcards here

n=int(float(config['CHANGE']['ISWEEP']['SAMPSIZE']))
ploidy=int(float(config['FIXED']['HAPIBD']['PLOIDY']))
maf1=float(config['FIXED']['HAPIBD']['MINMAF'])
mac1=int(ploidy*n*maf1)

# do this before running snakemake
# mkdir -p {config[CHANGE][FOLDERS][STUDY]}/{config[CHANGE][FOLDERS][ROI]}

rule subset_vcf: # focus vcf on region of interest
	input:
		vcfin='{config[CHANGE][FOLDERS][STUDY]}/vcfs/chr{config[CHANGE][ROI][CHRNUM]}.vcf.gz',
	output:
		subvcf='{config[CHANGE][FOLDERS][STUDY]}/{config[CHANGE][FOLDERS][ROI]}/chr{config[CHANGE][ROI][CHRNUM]}.vcf.gz',
	shell:
		'bash {config[CHANGE][FOLDERS][TERMINALSCRIPTS]}/roi-vcf.sh {config[CHANGE][FOLDERS][STUDY]}/{config[CHANGE][FOLDERS][ROI]} {input.vcfin} {output.vcfout} {config[CHANGE][ROI][LEFT]} {config[CHANGE][ROI][RIGHT]} {config[CHANGE][ROI][CHRNUM]} {config[FIXED][ISWEEP][MINAF]}'

rule hapibd: # segments from hap-ibd.jar
    input:
        vcf='{config[CHANGE][FOLDERS][STUDY]}/{config[CHANGE][FOLDERS][ROI]}/chr{config[CHANGE][ROI][CHRNUM]}.vcf.gz',
        map='{config[CHANGE][FOLDERS][STUDY]}/maps/chr{config[CHANGE][ROI][CHRNUM]}.map',
    params:
        minmac=str(mac1),
        out='{config[CHANGE][FOLDERS][STUDY]}/{config[CHANGE][FOLDERS][ROI]}/chr{config[CHANGE][ROI][CHRNUM]}',
    output:
        ibd='{config[CHANGE][FOLDERS][STUDY]}/{config[CHANGE][FOLDERS][ROI]}/chr{config[CHANGE][ROI][CHRNUM]}.ibd.gz',
        hbd='{config[CHANGE][FOLDERS][STUDY]}/{config[CHANGE][FOLDERS][ROI]}/chr{config[CHANGE][ROI][CHRNUM]}.hbd.gz',
        log='{config[CHANGE][FOLDERS][STUDY]}/{config[CHANGE][FOLDERS][ROI]}/chr{config[CHANGE][ROI][CHRNUM]}.log',
    shell:
        'java -jar {config[CHANGE][FOLDERS][SOFTWARE]}/{config[CHANGE][PROGRAMS][HAPIBD]} gt={input.vcf} map={input.map} out={params.out} min-seed={config[FIXED][HAPIBD][MINSEED]} min-extend={config[FIXED][HAPIBD][MINEXT]} min-output={config[FIXED][HAPIBD][MINOUT]} min-mac={params.minmac}'

rule filter_hapibd: # applying focus
    input:
        ibd='{config[CHANGE][FOLDERS][STUDY]}/{config[CHANGE][FOLDERS][ROI]}/chr{config[CHANGE][ROI][CHRNUM]}.ibd.gz',
    output:
        fipass='{config[CHANGE][FOLDERS][STUDY]}/{config[CHANGE][FOLDERS][ROI]}/focus.ibd.gz',
    shell: # if chromosome is huge (greater than 10000 Mb), may need to modify the third pipe
        'zcat {input.ibd} | java -jar {config[CHANGE][FOLDERS][SOFTWARE]}/{config[CHANGE][PROGRAMS][FILTER]} "I" 6 0.00 {config[CHANGE][ROI][CENTER]} | java -jar {config[CHANGE][FOLDERS][SOFTWARE]}/{config[CHANGE][PROGRAMS][FILTER]} "I" 7 {config[CHANGE][ROI][CENTER]} 10000000000 | gzip > {output.fipass}'

rule filter_mom:
	input:
		'{config[CHANGE][FOLDERS][STUDY]}/ibdsegs/ibdends/modified/mom/chr${config[CHANGE][ROI][CHRNUM]}.ibd.gz',
    output:
        fipass='{config[CHANGE][FOLDERS][STUDY]}/{config[CHANGE][FOLDERS][ROI]}/mom.ibd.gz',
    shell: # if chromosome is huge (greater than 10000 Mb), may need to modify the third pipe
        'zcat {input.ibd} | java -jar {config[CHANGE][FOLDERS][SOFTWARE]}/{config[CHANGE][PROGRAMS][FILTER]} "I" 6 0.00 {config[CHANGE][ROI][CENTER]} | java -jar {config[CHANGE][FOLDERS][SOFTWARE]}/{config[CHANGE][PROGRAMS][FILTER]} "I" 7 {config[CHANGE][ROI][CENTER]} 10000000000 | gzip > {output.fipass}'
