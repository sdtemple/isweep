# zooming in on a region of interest

# some inputs, string managements, count sample size
cohort=str(config['CHANGE']['FOLDERS']['STUDY'])
samplesize=0
with open(cohort+'/subsample.txt','r') as f:
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
        chr=$(python ../../scripts/utilities/lines.py {input.locus} 2 2)
        left=$(python ../../scripts/utilities/lines.py {input.locus} 4 2)
        right=$(python ../../scripts/utilities/lines.py {input.locus} 5 2)
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
        chr=$(python ../../scripts/utilities/lines.py {input.locus} 2 2)
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
    shell:
        """
        thecenter=$(python ../../scripts/utilities/lines.py {input.locus} 3 2)
        python ../../scripts/utilities/filter-lines.py \
            --input_file {input.ibd} \
            --output_file {wildcards.cohort}/{wildcards.hit}/intermediate.ibd.gz \
            --column_index 6 \
            --upper_bound $thecenter \
            --complement 0
        python ../../scripts/utilities/filter-lines.py \
            --input_file {wildcards.cohort}/{wildcards.hit}/intermediate.ibd.gz \
            --output_file {output.ibd} \
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
        python ../../scripts/model/rank.py \
            --input_ibd_file {input.short} \
            --input_vcf_file {input.vcf} \
            --output_file {output.fileout} \
            --graph_diameter {params.diameter} \
            --group_cutoff {params.rulesigma} \
            --lowest_freq {params.q1} \
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
        python ../../scripts/model/site.py \
            --input_snp_file {input.snps} \
            --output_folder {params.folderout} \
            --window_index 0 \
            --freq_index 1 \
            --window_size {params.windowsize} \
            --window_step {params.windowstep} \
            --quantile_sum {params.qrng} \
            --max_spacing {params.maxspace}
        """
