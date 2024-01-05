# zooming in on a region of interest

# some inputs, string managements, count sample size
subsamplefile=str(config['CHANGE']['ISWEEP']['SUBSAMPLE'])
cohort=str(config['CHANGE']['FOLDERS']['STUDY'])
samplesize=0
with open(cohort+'/'+subsamplefile,'r') as f:
    for line in f:
        samplesize+=1
ploidy=2
# ploidy=int(float(config['FIXED']['HAPIBD']['PLOIDY']))
maf3=float(config['FIXED']['HAPIBD']['MINMAF'])
mac3=int(ploidy*samplesize*maf3)
maf=float(config['FIXED']['ISWEEP']['MINAAF'])

# subset vcf to region of interest
rule first_region: # focus vcf on region of interest
    input:
        locus='{cohort}/{hit}/locus.txt',
        subsample="{cohort}/subsample.txt",
        excludesamples="{cohort}/excludesamples.txt",
    output:
        subvcf='{cohort}/{hit}/first.focused.vcf.gz',
    params:
        qmaf=maf,
        chrpre=str(config['CHANGE']['ISWEEP']['CHRPRE']),
        scripts=str(config['CHANGE']['FOLDERS']['TERMINALSCRIPTS']),
        vcfs=str(config['CHANGE']['EXISTING']['VCFS']),
        vcfprefix=str(config['CHANGE']['EXISTING']['VCFPRE']),
        vcfsuffix=str(config['CHANGE']['EXISTING']['VCFSUF']),
    resources:
        mem_gb='{config[CHANGE][ISWEEP][XMXMEM]}',
    shell: # if chromosome is huge (greater than 10000 Mb), may need to modify the third pipe
        """
        chr=$(python {params.scripts}/lines.py {input.locus} 2 2)
        left=$(python {params.scripts}/lines.py {input.locus} 4 2)
        right=$(python {params.scripts}/lines.py {input.locus} 5 2)
        vcf={params.vcfs}/{params.vcfprefix}${{chr}}{params.vcfsuffix}
        bcftools view ${{vcf}} -r {params.chrpre}${{chr}}:${{left}}-${{right}} -Ob | \
            bcftools view -S {input.subsample} -Ob | \
            bcftools view -q {params.qmaf}:nonmajor -Oz -o {output.subvcf}
        """

### call hap-ibd ###

rule first_hapibd:
    input:
        vcf='{cohort}/{hit}/first.focused.vcf.gz',
        locus='{cohort}/{hit}/locus.txt',
    params:
        minmac=str(mac3),
        out='{cohort}/{hit}/first',
        soft=str(config['CHANGE']['FOLDERS']['SOFTWARE']),
        prog=str(config['CHANGE']['PROGRAMS']['HAPIBD']),
        minsee=str(config['FIXED']['HAPIBD']['MINSEED']),
        minext=str(config['FIXED']['HAPIBD']['MINEXT']),
        minout=str(config['FIXED']['HAPIBD']['MINOUT']),
        scripts=str(config['CHANGE']['FOLDERS']['TERMINALSCRIPTS']),
    output:
        ibd='{cohort}/{hit}/first.ibd.gz',
        hbd='{cohort}/{hit}/first.hbd.gz',
        log='{cohort}/{hit}/first.log',
    resources:
        mem_gb='{config[CHANGE][ISWEEP][XMXMEM]}',
    shell:
        """
        chr=$(python {params.scripts}/lines.py {input.locus} 2 2)
        java -Xmx{config[CHANGE][ISWEEP][XMXMEM]}g -jar {params.soft}/{params.prog} \
            gt={input.vcf} \
            map={wildcards.cohort}/maps/chr${{chr}}.map \
            out={params.out} \
            min-seed={params.minsee} \
            min-extend={params.minext} \
            min-output={params.minout} \
            min-mac={params.minmac}
        """

### filter ibd file ###

rule first_filt:
    input:
        ibd='{cohort}/{hit}/first.ibd.gz',
        locus='{cohort}/{hit}/locus.txt',
    output:
        ibd='{cohort}/{hit}/first.filt.ibd.gz',
    params:
        soft=str(config['CHANGE']['FOLDERS']['SOFTWARE']),
        prog=str(config['CHANGE']['PROGRAMS']['FILTER']),
        scripts=str(config['CHANGE']['FOLDERS']['TERMINALSCRIPTS']),
    resources:
        mem_gb='{config[CHANGE][ISWEEP][XMXMEM]}'
    shell:
        """
        thecenter=$(python {params.scripts}/lines.py {input.locus} 3 2)
        zcat {input.ibd} | \
            java -Xmx{config[CHANGE][ISWEEP][XMXMEM]}g -jar {params.soft}/{params.prog} \
            "I" 6 0.00 $thecenter | \
            java -Xmx{config[CHANGE][ISWEEP][XMXMEM]}g -jar {params.soft}/{params.prog} \
            "I" 7 $thecenter 10000000000 | \
            gzip > {output.ibd}
        """

### rank snps ###

rule first_rank:
    input:
        short='{cohort}/{hit}/first.filt.ibd.gz',
        vcf='{cohort}/{hit}/first.focused.vcf.gz',
    output:
        fileout='{cohort}/{hit}/first.ranks.tsv.gz',
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
        snps='{cohort}/{hit}/first.ranks.tsv.gz',
    output:
        '{cohort}/{hit}/first.pos.txt',
        '{cohort}/{hit}/first.qs.tsv.gz',
        '{cohort}/{hit}/first.snp.png',
        '{cohort}/{hit}/first.score.png',
    params:
        windowsize=str(config['FIXED']['ISWEEP']['WINSIZE']),
        windowstep=str(config['FIXED']['ISWEEP']['WINSTEP']),
        qrng=str(config['FIXED']['ISWEEP']['QRANGE']),
        maxspace=str(config['FIXED']['ISWEEP']['MAXSPACING']),
        scripts=str(config['CHANGE']['FOLDERS']['TERMINALSCRIPTS']),
        folderout='{cohort}/{hit}',
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
