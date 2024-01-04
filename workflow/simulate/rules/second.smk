### refining location of analysis

# inputs
n=int(float(config['CHANGE']['SIMULATE']['SAMPSIZE']))
ploidy=int(float(config['FIXED']['SIMULATE']['PLOIDY']))
maf3=float(config['FIXED']['HAPIBD']['MINMAF'])
mac3=int(ploidy*n*maf3)

### filter vcf ###

rule second_region:
    input:
        vcfin='{macro}/{micro}/{seed}/large.chr1.vcf.gz',
        locus='{macro}/{micro}/{seed}/first.pos.txt',
    output:
        vcfout='{macro}/{micro}/{seed}/second.chr1.vcf.gz',
    params:
        folder='{macro}/{micro}/{seed}',
        pm=str(config['FIXED']['SIMULATE']['BUFFER']),
        script=str(config['CHANGE']['FOLDERS']['TERMINALSCRIPTS'])+'/lines.py',
    shell: # to bgz and back is being consertative
        """
        thecenter=$(python {params.script} {input.locus} 1 2)
        gunzip -c {input.vcfin} | bgzip  > {params.folder}/chrtemp.vcf.bgz
        tabix -fp vcf {params.folder}/chrtemp.vcf.bgz
        left=$(python -c "out = $thecenter - {params.pm} ; print(out)")
        right=$(python -c "out = $thecenter + {params.pm} ; print(out)")
        bcftools view {params.folder}/chrtemp.vcf.bgz \
            -r 1:${{left}}-${{right}} \
            -Oz -o {output.vcfout}
        rm {params.folder}/chrtemp.vcf.bgz
        """

### call hap-ibd ###

rule second_hapibd:
    input:
        vcf='{macro}/{micro}/{seed}/second.chr1.vcf.gz',
        map=str(config["CHANGE"]["FOLDERS"]["MACRO"]) + '/uniform.map',
    params:
        minmac=str(mac3),
        out='{macro}/{micro}/{seed}/second.chr1',
        soft=str(config['CHANGE']['FOLDERS']['SOFTWARE']),
        prog=str(config['CHANGE']['PROGRAMS']['HAPIBD']),
        minsee=str(config['FIXED']['HAPIBD']['MINSEED']),
        minext=str(config['FIXED']['HAPIBD']['MINEXT']),
        minout=str(config['FIXED']['HAPIBD']['MINOUT']),
    output:
        ibd='{macro}/{micro}/{seed}/second.chr1.ibd.gz',
        hbd='{macro}/{micro}/{seed}/second.chr1.hbd.gz',
        log='{macro}/{micro}/{seed}/second.chr1.log',
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

### filter ibd file ###

rule second_filt:
    input:
        ibd='{macro}/{micro}/{seed}/second.chr1.ibd.gz',
        locus='{macro}/{micro}/{seed}/first.pos.txt',
    output:
        ibd='{macro}/{micro}/{seed}/second.filt.chr1.ibd.gz',
    params:
        soft=str(config['CHANGE']['FOLDERS']['SOFTWARE']),
        prog=str(config['CHANGE']['PROGRAMS']['FILTER']),
        script=str(config['CHANGE']['FOLDERS']['TERMINALSCRIPTS'])+'/lines.py',
    resources:
        mem_gb='{config[CHANGE][CLUSTER][LARGEMEM]}'
    shell:
        """
        thecenter=$(python {params.script} {input.locus} 1 2)
        zcat {input.ibd} | \
            java -Xmx{config[CHANGE][CLUSTER][LARGEMEM]}g -jar {params.soft}/{params.prog} \
            "I" 6 0.00 $thecenter | \
            java -Xmx{config[CHANGE][CLUSTER][LARGEMEM]}g -jar {params.soft}/{params.prog} \
            "I" 7 $thecenter 10000000000 | \
            gzip > {output.ibd}
        """

### rank snps ###

rule second_rank:
    input:
        short='{macro}/{micro}/{seed}/second.filt.chr1.ibd.gz',
        vcf='{macro}/{micro}/{seed}/second.chr1.vcf.gz',
    output:
        fileout='{macro}/{micro}/{seed}/second.ranks.tsv.gz',
    params:
        scripts=str(config['CHANGE']['FOLDERS']['TERMINALSCRIPTS']),
        diameter=str(config['FIXED']['ISWEEP']['DIAMETER']),
        q1=str(config['FIXED']['ISWEEP']['MINAAF']),
        rulesigma=str(config['FIXED']['ISWEEP']['GROUPCUTOFF']),
    shell:
        """
        python {params.scripts}/rank.py \
            {input.short} \
            {input.vcf} \
            {output.fileout} \
            {params.diameter} \
            {params.q1} \
            {params.rulesigma}
        """
        
### write outliers ###

rule second_outlier:
    input:
        short='{macro}/{micro}/{seed}/second.filt.chr1.ibd.gz',
    output:
        fileout='{macro}/{micro}/{seed}/second.outliers.txt',
    params:
        scripts=str(config['CHANGE']['FOLDERS']['TERMINALSCRIPTS']),
        diameter=str(config['FIXED']['ISWEEP']['DIAMETER']),
        rulesigma=str(config['FIXED']['ISWEEP']['GROUPCUTOFF']),
    shell:
        """
        python {params.scripts}/outliers.py \
            {input.short} \
            {wildcards.macro}/{wildcards.micro}/{wildcards.seed} \
            {params.diameter} \
            {params.rulesigma}
        touch {output.fileout}
        """

rule ibd_entropy:
	input:
		filein='{macro}/{micro}/second.outliers.txt',
	output:
		fileout='{macro}/{micro}/ibd.entropy.tsv',
	params:
		scripts=str(config['CHANGE']['FOLDERS']['TERMINALSCRIPTS']),
		samplesize=str(n*ploidy),
	shell:
		"""
		python {params.scripts}/ibd-entropy.py \
			{wildcards.macro}/{wildcards.micro} \
			{output.fileout} \
			{params.samplesize}
		"""