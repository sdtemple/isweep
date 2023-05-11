wildcard_constraints:
	SIMNAME = '\w+',

### focus VCF file for short IBD ###

rule focus_region:
    input:
        vcfin='{macro}/{micro}/{seed}/large.chr1.vcf.gz',
    output:
        vcfout='{macro}/{micro}/{seed}/small.chr1.vcf.gz',
    params:
        left=str(config['FIXED']['SIMULATE']['LEFT']),
        right=str(config['FIXED']['SIMULATE']['RIGHT']),
        maf=str(config['CHANGE']['SIMULATE']['MSPMAF']),
        folder='{macro}/{micro}/{seed}',
    shell: # to bgz and back is being consertative
        """
        gunzip -c {input.vcfin} | bgzip  > {params.folder}/chrtemp.vcf.bgz
        tabix -fp vcf {params.folder}/chrtemp.vcf.bgz
        bcftools view {params.folder}/chrtemp.vcf.bgz \
            -r 1:{params.left}-{params.right} \
            -Oz -o {output.vcfout}
		rm {params.folder}/chrtemp.vcf.bgz
        """

### short IBD calls for focus region ###

n=int(float(config['CHANGE']['SIMULATE']['SAMPSIZE']))
ploidy=int(float(config['FIXED']['HAPIBD']['PLOIDY']))
maf3=float(config['FIXED']['HAPIBD']['MINMAF'])
mac3=int(ploidy*n*maf3)

rule small_hapibd:
    input:
        vcf='{macro}/{micro}/{seed}/small.chr1.vcf.gz',
        map=str(config["CHANGE"]["FOLDERS"]["MACRO"]) + '/uniform.map',
    params:
        minmac=str(mac3),
        out='{macro}/{micro}/{seed}/small.chr1',
        soft=str(config['CHANGE']['FOLDERS']['SOFTWARE']),
        prog=str(config['CHANGE']['PROGRAMS']['HAPIBD']),
        minsee=str(config['FIXED']['HAPIBD']['MINSEED']),
        minext=str(config['FIXED']['HAPIBD']['MINEXT']),
        minout=str(config['FIXED']['HAPIBD']['MINOUT']),
    output:
        ibd='{macro}/{micro}/{seed}/small.chr1.ibd.gz',
        hbd='{macro}/{micro}/{seed}/small.chr1.hbd.gz',
        log='{macro}/{micro}/{seed}/small.chr1.log',
    resources:
        mem_gb='{config[CHANGE][CLUSTER][LARGEMEM]}',
    shell:
        """
        java -Xmx{config[CHANGE][CLUSTER][LARGEMEM]}g -jar {params.soft}/{params.prog} \
            gt={input.vcf} \
            map={input.map} \
            out={params.out} \
            min-seed={params.minsee} \
            min-extend={params.minext} \
            min-output={params.minout} \
            min-mac={params.minmac}
        """
        # 'java -Xmx{config[CHANGE][CLUSTER][LARGEMEM]}g -jar {config[CHANGE][FOLDERS][SOFTWARE]}/{config[CHANGE][PROGRAMS][HAPIBD]} gt={input.vcf} map={input.map} out={params.out} min-seed={config[FIXED][HAPIBD][MINSEED]} min-extend={config[FIXED][HAPIBD][MINEXT]} min-output={config[FIXED][HAPIBD][MINOUT]} min-markers={config[FIXED][HAPIBD][MINMARK]} min-mac={params.minmac}'

### long IBD ###

rule long_ibd:
    input:
        ibd='{macro}/{micro}/{seed}/scan.chr1.ibd.gz',
		ibdwin='{macro}/{micro}/{seed}/scan.chr1.windowed.tsv.gz',
    output:
        ibd='{macro}/{micro}/{seed}/long.chr1.ibd.gz',
    params:
        soft=str(config['CHANGE']['FOLDERS']['SOFTWARE']),
        prog=str(config['CHANGE']['PROGRAMS']['FILTER']),
        # loc=str(config['FIXED']['SIMULATE']['LOC']),
		script=str(config['CHANGE']['FOLDERS']['TERMINALSCRIPTS'])+'/where-mode-ibd.py',
    resources:
        mem_gb='{config[CHANGE][CLUSTER][LARGEMEM]}'
    shell:
        """
		themode=$(python ${{params.script}} ${{input.ibdwin}})
        zcat {input.ibd} | \
            java -Xmx{config[CHANGE][CLUSTER][LARGEMEM]}g -jar {params.soft}/{params.prog} \
            "I" 6 0.00 $themode | \
            java -Xmx{config[CHANGE][CLUSTER][LARGEMEM]}g -jar {params.soft}/{params.prog} \
            "I" 7 $themode 10000000000 | \
            gzip > {output.ibd}
        """
        # """
        # zcat {input.ibd} | \
        #     java -Xmx{config[CHANGE][CLUSTER][LARGEMEM]}g -jar {params.soft}/{params.prog} \
        #     "I" 6 0.00 {params.loc} | \
        #     java -Xmx{config[CHANGE][CLUSTER][LARGEMEM]}g -jar {params.soft}/{params.prog} \
        #     "I" 7 {params.loc} 10000000000 | \
        #     gzip > {output.ibd}
        # """

### small IBD ###

rule small_ibd:
    input:
        ibd='{macro}/{micro}/{seed}/small.chr1.ibd.gz',
		ibdwin='{macro}/{micro}/{seed}/scan.chr1.windowed.tsv.gz',
    output:
        ibd='{macro}/{micro}/{seed}/short.chr1.ibd.gz',
    params:
        soft=str(config['CHANGE']['FOLDERS']['SOFTWARE']),
        prog=str(config['CHANGE']['PROGRAMS']['FILTER']),
        # loc=str(config['FIXED']['SIMULATE']['LOC']),
		script=str(config['CHANGE']['FOLDERS']['TERMINALSCRIPTS'])+'/where-mode-ibd.py',
    resources:
        mem_gb='{config[CHANGE][CLUSTER][LARGEMEM]}'
    shell:
        """
		themode=$(python ${{params.script}} ${{input.ibdwin}})
        zcat {input.ibd} | \
            java -Xmx{config[CHANGE][CLUSTER][LARGEMEM]}g -jar {params.soft}/{params.prog} \
            "I" 6 0.00 $themode | \
            java -Xmx{config[CHANGE][CLUSTER][LARGEMEM]}g -jar {params.soft}/{params.prog} \
            "I" 7 $themode 10000000000 | \
            gzip > {output.ibd}
        """
        # """
        # zcat {input.ibd} | \
        #     java -Xmx{config[CHANGE][CLUSTER][LARGEMEM]}g -jar {params.soft}/{params.prog} \
        #     "I" 6 0.00 {params.loc} | \
        #     java -Xmx{config[CHANGE][CLUSTER][LARGEMEM]}g -jar {params.soft}/{params.prog} \
        #     "I" 7 {params.loc} 10000000000 | \
        #     gzip > {output.ibd}
        # """
