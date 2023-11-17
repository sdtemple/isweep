wildcard_constraints:
	SIMNAME = '\w+',

### focus VCF file for short IBD ###

rule focus_region:
    input:
        vcfin='{macro}/{micro}/{seed}/large.chr1.vcf.gz',
    output:
        vcfout='{macro}/{micro}/{seed}/small.chr1.vcf.gz',
    params:
        loc=str(config['FIXED']['SIMULATE']['LOC']),
        pm=str(config['FIXED']['SIMULATE']['BUFFER']),
        folder='{macro}/{micro}/{seed}',
    shell: # to bgz and back is being consertative
        """
        gunzip -c {input.vcfin} | bgzip > {params.folder}/chrtemp.vcf.bgz
        tabix -fp vcf {params.folder}/chrtemp.vcf.bgz
        left=$(python -c "out = {params.loc} - {params.pm} ; print(out)")
        right=$(python -c "out = {params.loc} + {params.pm} ; print(out)")
        bcftools view {params.folder}/chrtemp.vcf.bgz \
            -r 1:${{left}}-${{right}} \
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
        loc=str(config['FIXED']['SIMULATE']['LOC']),
    resources:
        mem_gb='{config[CHANGE][CLUSTER][LARGEMEM]}'
    shell:
        """
        zcat {input.ibd} | \
            java -Xmx{config[CHANGE][CLUSTER][LARGEMEM]}g -jar {params.soft}/{params.prog} \
            "I" 6 0.00 {params.loc} | \
            java -Xmx{config[CHANGE][CLUSTER][LARGEMEM]}g -jar {params.soft}/{params.prog} \
            "I" 7 {params.loc} 10000000000 | \
            gzip > {output.ibd}
        """

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
        loc=str(config['FIXED']['SIMULATE']['LOC']),
    resources:
        mem_gb='{config[CHANGE][CLUSTER][LARGEMEM]}'
    shell:
        """
        zcat {input.ibd} | \
            java -Xmx{config[CHANGE][CLUSTER][LARGEMEM]}g -jar {params.soft}/{params.prog} \
            "I" 6 0.00 {params.loc} | \
            java -Xmx{config[CHANGE][CLUSTER][LARGEMEM]}g -jar {params.soft}/{params.prog} \
            "I" 7 {params.loc} 10000000000 | \
            gzip > {output.ibd}
        """

# mean, median, mode against true loc

rule focus_region_mean:
    input:
        vcfin='{macro}/{micro}/{seed}/large.chr1.vcf.gz',
        ibdwin='{macro}/{micro}/{seed}/scan.chr1.windowed.tsv.gz',
    output:
        vcfout='{macro}/{micro}/{seed}/small.chr1.mean.vcf.gz',
    params:
        folder='{macro}/{micro}/{seed}',
        pm=str(config['FIXED']['SIMULATE']['BUFFER']),
        script=str(config['CHANGE']['FOLDERS']['TERMINALSCRIPTS'])+'/test/where-mean-ibd.py',
    shell: # to bgz and back is being consertative
        """
        themean=$(python {params.script} {input.ibdwin})
        gunzip -c {input.vcfin} | bgzip  > {params.folder}/chrtemp.mean.v2.vcf.bgz
        tabix -fp vcf {params.folder}/chrtemp.mean.v2.vcf.bgz
        left=$(python -c "out = $themean - {params.pm} ; print(out)")
        right=$(python -c "out = $themean + {params.pm} ; print(out)")
        bcftools view {params.folder}/chrtemp.mean.v2.vcf.bgz \
            -r 1:${{left}}-${{right}} \
            -Oz -o {output.vcfout}
        rm {params.folder}/chrtemp.mean.v2.vcf.bgz
        """

rule focus_region_median:
    input:
        vcfin='{macro}/{micro}/{seed}/large.chr1.vcf.gz',
        ibdwin='{macro}/{micro}/{seed}/scan.chr1.windowed.tsv.gz',
    output:
        vcfout='{macro}/{micro}/{seed}/small.chr1.median.vcf.gz',
    params:
        folder='{macro}/{micro}/{seed}',
        pm=str(config['FIXED']['SIMULATE']['BUFFER']),
        script=str(config['CHANGE']['FOLDERS']['TERMINALSCRIPTS'])+'/test/where-median-ibd.py',
    shell: # to bgz and back is being consertative
        """
        themedian=$(python {params.script} {input.ibdwin})
        gunzip -c {input.vcfin} | bgzip  > {params.folder}/chrtemp.median.v2.vcf.bgz
        tabix -fp vcf {params.folder}/chrtemp.median.v2.vcf.bgz
        left=$(python -c "out = $themedian - {params.pm} ; print(out)")
        right=$(python -c "out = $themedian + {params.pm} ; print(out)")
        bcftools view {params.folder}/chrtemp.median.v2.vcf.bgz \
            -r 1:${{left}}-${{right}} \
            -Oz -o {output.vcfout}
        rm {params.folder}/chrtemp.median.v2.vcf.bgz
        """

rule focus_region_mode:
    input:
        vcfin='{macro}/{micro}/{seed}/large.chr1.vcf.gz',
        ibdwin='{macro}/{micro}/{seed}/scan.chr1.windowed.tsv.gz',
    output:
        vcfout='{macro}/{micro}/{seed}/small.chr1.mode.vcf.gz',
    params:
        folder='{macro}/{micro}/{seed}',
        pm=str(config['FIXED']['SIMULATE']['BUFFER']),
        script=str(config['CHANGE']['FOLDERS']['TERMINALSCRIPTS'])+'/test/where-mode-ibd.py',
    shell: # to bgz and back is being consertative
        """
        themode=$(python {params.script} {input.ibdwin})
        gunzip -c {input.vcfin} | bgzip  > {params.folder}/chrtemp.mode.v2.vcf.bgz
        tabix -fp vcf {params.folder}/chrtemp.mode.v2.vcf.bgz
        left=$(python -c "out = $themode - {params.pm} ; print(out)")
        right=$(python -c "out = $themode + {params.pm} ; print(out)")
        bcftools view {params.folder}/chrtemp.mode.v2.vcf.bgz \
            -r 1:${{left}}-${{right}} \
            -Oz -o {output.vcfout}
        rm {params.folder}/chrtemp.mode.v2.vcf.bgz
        """

rule small_hapibd_mode:
    input:
        vcf='{macro}/{micro}/{seed}/small.chr1.mode.vcf.gz',
        map=str(config["CHANGE"]["FOLDERS"]["MACRO"]) + '/uniform.map',
    params:
        minmac=str(mac3),
        out='{macro}/{micro}/{seed}/small.chr1.mode',
        soft=str(config['CHANGE']['FOLDERS']['SOFTWARE']),
        prog=str(config['CHANGE']['PROGRAMS']['HAPIBD']),
        minsee=str(config['FIXED']['HAPIBD']['MINSEED']),
        minext=str(config['FIXED']['HAPIBD']['MINEXT']),
        minout=str(config['FIXED']['HAPIBD']['MINOUT']),
    output:
        ibd='{macro}/{micro}/{seed}/small.chr1.mode.ibd.gz',
        hbd='{macro}/{micro}/{seed}/small.chr1.mode.hbd.gz',
        log='{macro}/{micro}/{seed}/small.chr1.mode.log',
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

rule small_hapibd_mean:
    input:
        vcf='{macro}/{micro}/{seed}/small.chr1.mean.vcf.gz',
        map=str(config["CHANGE"]["FOLDERS"]["MACRO"]) + '/uniform.map',
    params:
        minmac=str(mac3),
        out='{macro}/{micro}/{seed}/small.chr1.mean',
        soft=str(config['CHANGE']['FOLDERS']['SOFTWARE']),
        prog=str(config['CHANGE']['PROGRAMS']['HAPIBD']),
        minsee=str(config['FIXED']['HAPIBD']['MINSEED']),
        minext=str(config['FIXED']['HAPIBD']['MINEXT']),
        minout=str(config['FIXED']['HAPIBD']['MINOUT']),
    output:
        ibd='{macro}/{micro}/{seed}/small.chr1.mean.ibd.gz',
        hbd='{macro}/{micro}/{seed}/small.chr1.mean.hbd.gz',
        log='{macro}/{micro}/{seed}/small.chr1.mean.log',
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

rule small_hapibd_median:
    input:
        vcf='{macro}/{micro}/{seed}/small.chr1.median.vcf.gz',
        map=str(config["CHANGE"]["FOLDERS"]["MACRO"]) + '/uniform.map',
    params:
        minmac=str(mac3),
        out='{macro}/{micro}/{seed}/small.chr1.median',
        soft=str(config['CHANGE']['FOLDERS']['SOFTWARE']),
        prog=str(config['CHANGE']['PROGRAMS']['HAPIBD']),
        minsee=str(config['FIXED']['HAPIBD']['MINSEED']),
        minext=str(config['FIXED']['HAPIBD']['MINEXT']),
        minout=str(config['FIXED']['HAPIBD']['MINOUT']),
    output:
        ibd='{macro}/{micro}/{seed}/small.chr1.median.ibd.gz',
        hbd='{macro}/{micro}/{seed}/small.chr1.median.hbd.gz',
        log='{macro}/{micro}/{seed}/small.chr1.median.log',
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

rule long_ibd_mean:
    input:
        ibd='{macro}/{micro}/{seed}/scan.chr1.ibd.gz',
        ibdwin='{macro}/{micro}/{seed}/scan.chr1.windowed.tsv.gz',
    output:
        ibd='{macro}/{micro}/{seed}/long.chr1.mean.ibd.gz',
    params:
        soft=str(config['CHANGE']['FOLDERS']['SOFTWARE']),
        prog=str(config['CHANGE']['PROGRAMS']['FILTER']),
        script=str(config['CHANGE']['FOLDERS']['TERMINALSCRIPTS'])+'/test/where-mean-ibd.py',
    resources:
        mem_gb='{config[CHANGE][CLUSTER][LARGEMEM]}'
    shell:
        """
        themean=$(python {params.script} {input.ibdwin})
        zcat {input.ibd} | \
            java -Xmx{config[CHANGE][CLUSTER][LARGEMEM]}g -jar {params.soft}/{params.prog} \
            "I" 6 0.00 ${{themean}} | \
            java -Xmx{config[CHANGE][CLUSTER][LARGEMEM]}g -jar {params.soft}/{params.prog} \
            "I" 7 ${{themean}} 10000000000 | \
            gzip > {output.ibd}
        """

rule long_ibd_median:
    input:
        ibd='{macro}/{micro}/{seed}/scan.chr1.ibd.gz',
        ibdwin='{macro}/{micro}/{seed}/scan.chr1.windowed.tsv.gz',
    output:
        ibd='{macro}/{micro}/{seed}/long.chr1.median.ibd.gz',
    params:
        soft=str(config['CHANGE']['FOLDERS']['SOFTWARE']),
        prog=str(config['CHANGE']['PROGRAMS']['FILTER']),
        script=str(config['CHANGE']['FOLDERS']['TERMINALSCRIPTS'])+'/test/where-median-ibd.py',
    resources:
        mem_gb='{config[CHANGE][CLUSTER][LARGEMEM]}'
    shell:
        """
        themedian=$(python {params.script} {input.ibdwin})
        zcat {input.ibd} | \
            java -Xmx{config[CHANGE][CLUSTER][LARGEMEM]}g -jar {params.soft}/{params.prog} \
            "I" 6 0.00 ${{themedian}} | \
            java -Xmx{config[CHANGE][CLUSTER][LARGEMEM]}g -jar {params.soft}/{params.prog} \
            "I" 7 ${{themedian}} 10000000000 | \
            gzip > {output.ibd}
        """

rule long_ibd_mode:
    input:
        ibd='{macro}/{micro}/{seed}/scan.chr1.ibd.gz',
        ibdwin='{macro}/{micro}/{seed}/scan.chr1.windowed.tsv.gz',
    output:
        ibd='{macro}/{micro}/{seed}/long.chr1.mode.ibd.gz',
    params:
        soft=str(config['CHANGE']['FOLDERS']['SOFTWARE']),
        prog=str(config['CHANGE']['PROGRAMS']['FILTER']),
        script=str(config['CHANGE']['FOLDERS']['TERMINALSCRIPTS'])+'/test/where-mode-ibd.py',
    resources:
        mem_gb='{config[CHANGE][CLUSTER][LARGEMEM]}'
    shell:
        """
        themode=$(python {params.script} {input.ibdwin})
        zcat {input.ibd} | \
            java -Xmx{config[CHANGE][CLUSTER][LARGEMEM]}g -jar {params.soft}/{params.prog} \
            "I" 6 0.00 ${{themode}} | \
            java -Xmx{config[CHANGE][CLUSTER][LARGEMEM]}g -jar {params.soft}/{params.prog} \
            "I" 7 ${{themode}} 10000000000 | \
            gzip > {output.ibd}
        """

rule small_ibd_mode:
    input:
        ibd='{macro}/{micro}/{seed}/small.chr1.mode.ibd.gz',
        ibdwin='{macro}/{micro}/{seed}/scan.chr1.windowed.tsv.gz',
    output:
        ibd='{macro}/{micro}/{seed}/short.chr1.mode.ibd.gz',
    params:
        soft=str(config['CHANGE']['FOLDERS']['SOFTWARE']),
        prog=str(config['CHANGE']['PROGRAMS']['FILTER']),
        script=str(config['CHANGE']['FOLDERS']['TERMINALSCRIPTS'])+'/test/where-mode-ibd.py',
    resources:
        mem_gb='{config[CHANGE][CLUSTER][LARGEMEM]}'
    shell:
        """
        themode=$(python {params.script} {input.ibdwin})
        zcat {input.ibd} | \
            java -Xmx{config[CHANGE][CLUSTER][LARGEMEM]}g -jar {params.soft}/{params.prog} \
            "I" 6 0.00 $themode | \
            java -Xmx{config[CHANGE][CLUSTER][LARGEMEM]}g -jar {params.soft}/{params.prog} \
            "I" 7 $themode 10000000000 | \
            gzip > {output.ibd}
        """

rule small_ibd_mean:
    input:
        ibd='{macro}/{micro}/{seed}/small.chr1.mean.ibd.gz',
        ibdwin='{macro}/{micro}/{seed}/scan.chr1.windowed.tsv.gz',
    output:
        ibd='{macro}/{micro}/{seed}/short.chr1.mean.ibd.gz',
    params:
        soft=str(config['CHANGE']['FOLDERS']['SOFTWARE']),
        prog=str(config['CHANGE']['PROGRAMS']['FILTER']),
        script=str(config['CHANGE']['FOLDERS']['TERMINALSCRIPTS'])+'/test/where-mean-ibd.py',
    resources:
        mem_gb='{config[CHANGE][CLUSTER][LARGEMEM]}'
    shell:
        """
        themean=$(python {params.script} {input.ibdwin})
        zcat {input.ibd} | \
            java -Xmx{config[CHANGE][CLUSTER][LARGEMEM]}g -jar {params.soft}/{params.prog} \
            "I" 6 0.00 $themean | \
            java -Xmx{config[CHANGE][CLUSTER][LARGEMEM]}g -jar {params.soft}/{params.prog} \
            "I" 7 $themean 10000000000 | \
            gzip > {output.ibd}
        """

rule small_ibd_median:
    input:
        ibd='{macro}/{micro}/{seed}/small.chr1.median.ibd.gz',
        ibdwin='{macro}/{micro}/{seed}/scan.chr1.windowed.tsv.gz',
    output:
        ibd='{macro}/{micro}/{seed}/short.chr1.median.ibd.gz',
    params:
        soft=str(config['CHANGE']['FOLDERS']['SOFTWARE']),
        prog=str(config['CHANGE']['PROGRAMS']['FILTER']),
        script=str(config['CHANGE']['FOLDERS']['TERMINALSCRIPTS'])+'/test/where-median-ibd.py',
    resources:
        mem_gb='{config[CHANGE][CLUSTER][LARGEMEM]}'
    shell:
        """
        themedian=$(python {params.script} {input.ibdwin})
        zcat {input.ibd} | \
            java -Xmx{config[CHANGE][CLUSTER][LARGEMEM]}g -jar {params.soft}/{params.prog} \
            "I" 6 0.00 $themedian | \
            java -Xmx{config[CHANGE][CLUSTER][LARGEMEM]}g -jar {params.soft}/{params.prog} \
            "I" 7 $themedian 10000000000 | \
            gzip > {output.ibd}
        """

### narrowing in ###

# rank polymorphisms in focus region
rule init_rank:
    input:
        short='{macro}/{micro}/{seed}/short.chr1.ibd.gz',
        vcf='{macro}/{micro}/{seed}/small.chr1.vcf.gz',
    output:
        fileout='{macro}/{micro}/{seed}/init.ranks.tsv.gz',
    params:
        scripts=str(config['CHANGE']['FOLDERS']['TERMINALSCRIPTS']),
        diameter=str(config['FIXED']['ISWEEP']['DIAMETER']),
        q1=str(config['FIXED']['ISWEEP']['MINMAXAAF']),
        rulesigma=str(config['FIXED']['ISWEEP']['RULESIGMA']),
    shell:
        """
        python {params.scripts}/rank-isweep.py \
            {input.short} \
            {input.vcf} \
            {output.fileout} \
            {params.diameter} \
            {params.q1} \
            {params.rulesigma}
        """

rule qscore:
    input:
        snps='{macro}/{micro}/{seed}/init.ranks.tsv.gz',
        conv='{macro}/{micro}/{seed}/scan.chr1.windowed.tsv.gz',
    output:
        wins='{macro}/{micro}/{seed}/qtable.tsv.gz',
        loci='{macro}/{micro}/{seed}/qscore.pos.txt',
    params:
        windowsize=str(config['FIXED']['ISWEEP']['WINSIZE']),
        windowstep=str(config['FIXED']['ISWEEP']['WINSTEP']),
        windowcol=str(config['FIXED']['ISWEEP']['WINCOL']),
        freqcol=str(config['FIXED']['ISWEEP']['FREQCOL']),
        qrng=str(config['FIXED']['ISWEEP']['QRANGE']),
        maxspace=str(config['FIXED']['ISWEEP']['MAXSPACING'])
        scripts=str(config['CHANGE']['FOLDERS']['TERMINALSCRIPTS']),
        folderout='{macro}/{micro}/{seed}',
    shell:
        """
        python {params.scripts}/refine-site.py \
            {input.snps} \
            {input.conv} \
            {params.folderout} \
            0 \
            BPWINDOW \
            CMWINDOW \
            {params.windowsize} \
            {params.windowstep} \
            {params.windowcol} \
            {params.freqcol} \
            {params.qrng} \
            {params.maxspace}
        """

rule qscore_ibd:
    input:
        ibd='{macro}/{micro}/{seed}/small.chr1.ibd.gz',
        locus='{macro}/{micro}/{seed}/qscore.pos.txt',
    output:
        ibd='{macro}/{micro}/{seed}/short.chr1.qscore.ibd.gz',
    params:
        soft=str(config['CHANGE']['FOLDERS']['SOFTWARE']),
        prog=str(config['CHANGE']['PROGRAMS']['FILTER']),
        script=str(config['CHANGE']['FOLDERS']['TERMINALSCRIPTS'])+'/first-line.py',
    resources:
        mem_gb='{config[CHANGE][CLUSTER][LARGEMEM]}'
    shell:
        """
        themean=$(python {params.script} {input.locus})
        zcat {input.ibd} | \
            java -Xmx{config[CHANGE][CLUSTER][LARGEMEM]}g -jar {params.soft}/{params.prog} \
            "I" 6 0.00 ${{themean}} | \
            java -Xmx{config[CHANGE][CLUSTER][LARGEMEM]}g -jar {params.soft}/{params.prog} \
            "I" 7 ${{themean}} 10000000000 | \
            gzip > {output.ibd}
        """

# rank polymorphisms in focus region
rule qscore_rank:
    input:
        short='{macro}/{micro}/{seed}/short.chr1.qscore.ibd.gz',
        vcf='{macro}/{micro}/{seed}/small.chr1.vcf.gz',
    output:
        fileout='{macro}/{micro}/{seed}/qscore.ranks.tsv.gz',
    params:
        scripts=str(config['CHANGE']['FOLDERS']['TERMINALSCRIPTS']),
        diameter=str(config['FIXED']['ISWEEP']['DIAMETER']),
        q1=str(config['FIXED']['ISWEEP']['MINMAXAAF']),
        rulesigma=str(config['FIXED']['ISWEEP']['RULESIGMA']),
    shell:
        """
        python {params.scripts}/rank-isweep.py \
            {input.short} \
            {input.vcf} \
            {output.fileout} \
            {params.diameter} \
            {params.q1} \
            {params.rulesigma}
        """
