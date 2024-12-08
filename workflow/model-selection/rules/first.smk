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
        vcfs=str(config['CHANGE']['EXISTING']['VCFS']),
        vcfprefix=str(config['CHANGE']['EXISTING']['VCFPRE']),
        vcfsuffix=str(config['CHANGE']['EXISTING']['VCFSUF']),
    resources:
        mem_gb='{config[CHANGE][ISWEEP][XMXMEM]}',
    shell: # if chromosome is huge (greater than 10000 Mb), may need to modify the third pipe
        """
        chr=$(python ../../scripts/lines.py {input.locus} 2 2)
        left=$(python ../../scripts/lines.py {input.locus} 4 2)
        right=$(python ../../scripts/lines.py {input.locus} 5 2)
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
        minsee=str(config['FIXED']['HAPIBD']['MINSEED']),
        minext=str(config['FIXED']['HAPIBD']['MINEXT']),
        minout=str(config['FIXED']['HAPIBD']['MINOUT']),
    output:
        ibd='{cohort}/{hit}/first.ibd.gz',
        hbd='{cohort}/{hit}/first.hbd.gz',
        log='{cohort}/{hit}/first.log',
    resources:
        mem_gb='{config[CHANGE][ISWEEP][XMXMEM]}',
    shell:
        """
        chr=$(python ../../scripts/lines.py {input.locus} 2 2)
        java -Xmx{config[CHANGE][ISWEEP][XMXMEM]}g -jar ../../software/hap-ibd.jar \
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
    resources:
        mem_gb='{config[CHANGE][ISWEEP][XMXMEM]}'
    shell:
        """
        thecenter=../../scripts/lines.py {input.locus} 3 2)
        python ../../scripts/filter-lines.py {input.ibd} \
            {wildcards.cohort}/{wildcards.hit}/intermediate.ibd.gz \
            --column_index 6 \
            --upper_bound $thecenter \
            --complement 0
        python ../../scripts/filter-lines.py \
            {wildcards.cohort}/{wildcards.hit}/intermediate.ibd.gz \
            {output.ibd} \
            --column_index 7 \
            --lower_bound $thecenter \
            --upper_bound 10000000000 \
            --complement 0
        rm {wildcards.cohort}/{wildcards.hit}/intermediate.ibd.gz
        """

### rank snps ###

rule first_rank:
    input:
        short='{cohort}/{hit}/first.filt.ibd.gz',
        vcf='{cohort}/{hit}/first.focused.vcf.gz',
    output:
        fileout='{cohort}/{hit}/first.ranks.tsv.gz',
    params:
        diameter=str(config['FIXED']['ISWEEP']['DIAMETER']),
        q1=str(config['FIXED']['ISWEEP']['MINAAF']),
        rulesigma=str(config['FIXED']['ISWEEP']['GROUPCUTOFF']),
    shell:
        """
        python ../../scripts/rank.py \
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
        folderout='{cohort}/{hit}',
    shell:
        """
        python ../../scripts/site.py \
            {input.snps} \
            {params.folderout} \
            0 \
            1 \
            {params.windowsize} \
            {params.windowstep} \
            {params.qrng} \
            {params.maxspace}
        """
