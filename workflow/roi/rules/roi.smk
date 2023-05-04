# preparing regions of interest
# seth temple, sdtemple@uw.edu
# may 3, 2023

# some inputs, string managements, count sample size
subsamplefile=str(config['CHANGE']['ISWEEP']['SUBSAMPLE'])
macro=str(config['CHANGE']['FOLDERS']['STUDY'])
samplesize=0
with open(macro+'/'+subsamplefile,'r') as f:
    for line in f:
        samplesize+=1
samplesize=str(samplesize)
ploidy=int(float(config['FIXED']['HAPIBD']['PLOIDY']))
maf1=float(config['FIXED']['HAPIBD']['MINMAF'])
mac1=int(ploidy*int(samplesize)*maf1)
maf2=float(config['FIXED']['ISWEEP']['MINAAF'])
maf=min(maf1,maf2)
macro=str(config['CHANGE']['FOLDERS']['STUDY'])
vcffolder=str(config['CHANGE']['EXISTING']['VCFS'])
vcfpre=str(config['CHANGE']['EXISTING']['VCFPRE'])
vcfsuf=str(config['CHANGE']['EXISTING']['VCFSUF'])

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
        folder='{cohort}/{roi}/chr{chr}/center{center}/left{left}/right{right}/empty.txt',
        vcfin=vcffolder + '/' + vcfpre + '{chr}' + vcfsuf,
        subsample=macro+"/subsample.txt",
        excludesamples=macro+"/excludesamples.txt",
    output:
        subvcf='{cohort}/{roi}/chr{chr}/center{center}/left{left}/right{right}/chr.focused.vcf.gz',
    params:
        qmaf=maf,
        chrpre=str(config['CHANGE']['ISWEEP']['CHRPRE']),
    resources:
        mem_gb=10
    shell: # if chromosome is huge (greater than 10000 Mb), may need to modify the third pipe
        """
        tabix -fp vcf {input.vcfin}
        bcftools view {input.vcfin} -r {params.chrpre}{wildcards.chr}:{wildcards.left}-{wildcards.right} -Ob | \
            bcftools view -S {input.subsample} -Ob | \
            bcftools view -q {params.qmaf}:nonmajor -Oz -o {output.subvcf}
        """
        # 'bash {config[CHANGE][FOLDERS][TERMINALSCRIPTS]}/roi-vcf.sh {input.vcfin} {output.subvcf} {wildcards.left} {wildcards.right} {input.subsample} {params.qmaf} {wildcards.chr} {params.chrpre}'

# perform hap-ibd with sequence parameter settings
rule hapibd: # segments from hap-ibd.jar
    input:
        vcf='{cohort}/{roi}/chr{chr}/center{center}/left{left}/right{right}/chr.focused.vcf.gz',
        map='{cohort}/maps/chr{chr}.map',
        excludesamples=macro+"/excludesamples.txt",
    params:
        minmac=str(mac1),
        out='{cohort}/{roi}/chr{chr}/center{center}/left{left}/right{right}/chr',
        soft=str(config['CHANGE']['FOLDERS']['SOFTWARE']),
        prog=str(config['CHANGE']['PROGRAMS']['HAPIBD']),
        xmx=str(config['CHANGE']['ISWEEP']['XMXMEM']),
        minsee=str(config['FIXED']['HAPIBD']['MINSEED']),
        minext=str(config['FIXED']['HAPIBD']['MINEXT']),
        minout=str(config['FIXED']['HAPIBD']['MINOUT']),
    output:
        ibd='{cohort}/{roi}/chr{chr}/center{center}/left{left}/right{right}/chr.ibd.gz',
    resources:
        mem_gb=100
    shell:
        """
        java -Xmx{params.xmx}g -jar {params.soft}/{params.prog} \
            gt={input.vcf} \
            map={input.map} \
            out={params.out} \
            min-seed={params.minsee} \
            min-extend={params.minext} \
            min-output={params.minout} \
            min-mac={params.minmac} \
            excludesamples={input.excludesamples}
        """
        # 'java -Xmx{config[CHANGE][ISWEEP][XMXMEM]}g -jar {config[CHANGE][FOLDERS][SOFTWARE]}/{config[CHANGE][PROGRAMS][HAPIBD]} gt={input.vcf} map={input.map} out={params.out} min-seed={config[FIXED][HAPIBD][MINSEED]} min-extend={config[FIXED][HAPIBD][MINEXT]} min-output={config[FIXED][HAPIBD][MINOUT]} min-mac={params.minmac} excludesamples={input.excludesamples}'

# filter outputs of hap-ibd to be about focal location
rule filter_hapibd: # applying focus
    input:
        ibd='{cohort}/{roi}/chr{chr}/center{center}/left{left}/right{right}/chr.ibd.gz',
    output:
        fipass='{cohort}/{roi}/chr{chr}/center{center}/left{left}/right{right}/focus.ibd.gz',
    params:
        soft=str(config['CHANGE']['FOLDERS']['SOFTWARE']),
        prog=str(config['CHANGE']['PROGRAMS']['FILTER']),
        xmx=str(config['CHANGE']['ISWEEP']['XMXMEM']),
    resources:
        mem_gb=10
    shell: # if chromosome is huge (greater than 10000 Mb), may need to modify the third pipe
        """
        zcat {input.ibd} | \
            java -Xmx{params.xmx}g -jar {params.soft}/{params.prog} "I" 6 0.00 {wildcards.center} | \
            java -Xmx{params.xmx}g -jar {params.soft}/{params.prog} "I" 7 {wildcards.center} 10000000000 | \
            gzip > {output.fipass}; rm {input.ibd}
        """
        # 'zcat {input.ibd} | java -Xmx{config[CHANGE][ISWEEP][XMXMEM]}g -jar {config[CHANGE][FOLDERS][SOFTWARE]}/{config[CHANGE][PROGRAMS][FILTER]} "I" 6 0.00 {wildcards.center} | java -Xmx{config[CHANGE][ISWEEP][XMXMEM]}g -jar {config[CHANGE][FOLDERS][SOFTWARE]}/{config[CHANGE][PROGRAMS][FILTER]} "I" 7 {wildcards.center} 10000000000 | gzip > {output.fipass}; rm {input.ibd}'

# filter method of moments to be about focal location
rule filter_mom:
    input:
        ibd='{cohort}/ibdsegs/ibdends/modified/mom/chr{chr}.ibd.gz',
        ibdpseudo='{cohort}/{roi}/chr{chr}/center{center}/left{left}/right{right}/focus.ibd.gz',
    output:
        fipass='{cohort}/{roi}/chr{chr}/center{center}/left{left}/right{right}/mom.ibd.gz',
        # doneso='{cohort}/{roi}/chr{chr}/center{center}/left{left}/right{right}/done.txt',
    params:
        soft=str(config['CHANGE']['FOLDERS']['SOFTWARE']),
        prog=str(config['CHANGE']['PROGRAMS']['FILTER']),
        xmx=str(config['CHANGE']['ISWEEP']['XMXMEM']),
    resources:
        mem_gb=10
    shell: # if chromosome is huge (greater than 10000 Mb), may need to modify the third pipe
        """
        zcat {wildcards.cohort}/ibdsegs/ibdends/modified/mom/chr{wildcards.chr}.ibd.gz | \
            java -Xmx{params.xmx}g -jar {params.soft}/{params.prog} "I" 6 0.00 {wildcards.center} | \
            java -Xmx{params.xmx}g -jar {params.soft}/{params.prog} "I" 7 {wildcards.center} 10000000000 | \
            gzip > {output.fipass}
        rm {input.ibd}
        """
        # 'zcat {wildcards.cohort}/ibdsegs/ibdends/modified/mom/chr{wildcards.chr}.ibd.gz | java -Xmx{config[CHANGE][ISWEEP][XMXMEM]}g -jar {config[CHANGE][FOLDERS][SOFTWARE]}/{config[CHANGE][PROGRAMS][FILTER]} "I" 6 0.00 {wildcards.center} | java -Xmx{config[CHANGE][ISWEEP][XMXMEM]}g -jar {config[CHANGE][FOLDERS][SOFTWARE]}/{config[CHANGE][PROGRAMS][FILTER]} "I" 7 {wildcards.center} 10000000000 | gzip > {output.fipass}'
