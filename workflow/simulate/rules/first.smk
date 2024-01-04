### initializing regions of interest
### some localizing

# inputs
n=int(float(config['CHANGE']['SIMULATE']['SAMPSIZE']))
ploidy=int(float(config['FIXED']['SIMULATE']['PLOIDY']))
maf3=float(config['FIXED']['HAPIBD']['MINMAF'])
mac3=int(ploidy*n*maf3)

### modal location

rule first_mode:
    input:
        ibd='{macro}/{micro}/{seed}/scan.chr1.windowed.tsv.gz',
    output:
        ibd='{macro}/{micro}/{seed}/first.mode.txt',
    params:
        soft=str(config['CHANGE']['FOLDERS']['SOFTWARE']),
        prog=str(config['CHANGE']['PROGRAMS']['FILTER']),
        script=str(config['CHANGE']['FOLDERS']['TERMINALSCRIPTS'])+'/where-mode-ibd.py',
    resources:
        mem_gb='{config[CHANGE][CLUSTER][LARGEMEM]}'
    shell:
        """
        python {params.script} {input.ibd} {output.ibd}
        """

### filter vcf ###

rule first_region:
    input:
        vcfin='{macro}/{micro}/{seed}/large.chr1.vcf.gz',
        ibdwin='{macro}/{micro}/{seed}/first.mode.txt',
    output:
        vcfout='{macro}/{micro}/{seed}/first.chr1.vcf.gz',
    params:
        folder='{macro}/{micro}/{seed}',
        pm=str(config['FIXED']['SIMULATE']['BUFFER']),
        script=str(config['CHANGE']['FOLDERS']['TERMINALSCRIPTS'])+'/lines.py',
    shell: # to bgz and back is being consertative
        """
        thecenter=$(python {params.script} {input.ibdwin} 1 2)
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

rule first_hapibd:
    input:
        vcf='{macro}/{micro}/{seed}/first.chr1.vcf.gz',
        map=str(config["CHANGE"]["FOLDERS"]["MACRO"]) + '/uniform.map',
    params:
        minmac=str(mac3),
        out='{macro}/{micro}/{seed}/first.chr1',
        soft=str(config['CHANGE']['FOLDERS']['SOFTWARE']),
        prog=str(config['CHANGE']['PROGRAMS']['HAPIBD']),
        minsee=str(config['FIXED']['HAPIBD']['MINSEED']),
        minext=str(config['FIXED']['HAPIBD']['MINEXT']),
        minout=str(config['FIXED']['HAPIBD']['MINOUT']),
    output:
        ibd='{macro}/{micro}/{seed}/first.chr1.ibd.gz',
        hbd='{macro}/{micro}/{seed}/first.chr1.hbd.gz',
        log='{macro}/{micro}/{seed}/first.chr1.log',
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

rule first_filt:
    input:
        ibd='{macro}/{micro}/{seed}/first.chr1.ibd.gz',
        ibdwin='{macro}/{micro}/{seed}/first.mode.txt',
    output:
        ibd='{macro}/{micro}/{seed}/first.filt.chr1.ibd.gz',
    params:
        soft=str(config['CHANGE']['FOLDERS']['SOFTWARE']),
        prog=str(config['CHANGE']['PROGRAMS']['FILTER']),
        script=str(config['CHANGE']['FOLDERS']['TERMINALSCRIPTS'])+'/lines.py',
    resources:
        mem_gb='{config[CHANGE][CLUSTER][LARGEMEM]}'
    shell:
        """
        thecenter=$(python {params.script} {input.ibdwin} 1 2)
        zcat {input.ibd} | \
            java -Xmx{config[CHANGE][CLUSTER][LARGEMEM]}g -jar {params.soft}/{params.prog} \
            "I" 6 0.00 $thecenter | \
            java -Xmx{config[CHANGE][CLUSTER][LARGEMEM]}g -jar {params.soft}/{params.prog} \
            "I" 7 $thecenter 10000000000 | \
            gzip > {output.ibd}
        """

### rank snps ###

rule first_rank:
    input:
        short='{macro}/{micro}/{seed}/first.filt.chr1.ibd.gz',
        vcf='{macro}/{micro}/{seed}/first.chr1.vcf.gz',
    output:
        fileout='{macro}/{micro}/{seed}/first.ranks.tsv.gz',
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

rule first_score:
    input:
        snps='{macro}/{micro}/{seed}/first.ranks.tsv.gz',
    output:
        loci='{macro}/{micro}/{seed}/first.pos.txt',
    params:
        windowsize=str(config['FIXED']['ISWEEP']['WINSIZE']),
        windowstep=str(config['FIXED']['ISWEEP']['WINSTEP']),
        qrng=str(config['FIXED']['ISWEEP']['QRANGE']),
        maxspace=str(config['FIXED']['ISWEEP']['MAXSPACING']),
        scripts=str(config['CHANGE']['FOLDERS']['TERMINALSCRIPTS']),
        folderout='{macro}/{micro}/{seed}',
    shell:
        """
        python {params.scripts}/site.py \
            {input.snps} \
            {params.folderout} \
            0 \
            1 \
            {params.windowsize} \
            {params.windowstep} \
            {params.qrng} \
            {params.maxspace}
        """
