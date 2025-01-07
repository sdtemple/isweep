# Zooming in on a region of interest

maf3=float(config['fixed']['hap_ibd']['min_minor_allele_frequency'])
mac3=int(ploidy*samplesize*maf3)

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
        chrpre=str(config['change']['files']['chromosome_prefix']),
        vcfs=str(config['change']['files']['vcfs']),
        vcfprefix=str(config['change']['files']['vcf_prefix']),
        vcfsuffix=str(config['change']['files']['vcf_suffix']),
    resources:
        mem_gb='{config[change][isweep][xmx_mem]}',
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

rule first_hap_ibd:
    input:
        vcf='{cohort}/{hit}/first.focused.vcf.gz',
        locus='{cohort}/{hit}/locus.txt',
    params:
        minmac=str(mac3),
        out='{cohort}/{hit}/first',
        minsee=str(config['fixed']['hap_ibd']['min_seed']),
        minext=str(config['fixed']['hap_ibd']['min_extend']),
        minout=str(config['fixed']['hap_ibd']['min_output']),
    output:
        ibd='{cohort}/{hit}/first.ibd.gz',
        hbd='{cohort}/{hit}/first.hbd.gz',
        log='{cohort}/{hit}/first.log',
    resources:
        mem_gb='{config[change][isweep][xmx_mem]}',
    shell:
        """
        chr=$(python ../../scripts/utilities/lines.py {input.locus} 2 2)
        java -Xmx{config[change][isweep][xmx_mem]}g -jar ../../software/hap-ibd.jar \
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
        diameter=diameter,
        q1=maf,
        rulesigma=group_cutoff,
        ploidy=str(ploidy),
    shell:
        """
        python ../../scripts/model/rank.py \
            --input_ibd_file {input.short} \
            --input_vcf_file {input.vcf} \
            --output_file {output.fileout} \
            --graph_diameter {params.diameter} \
            --group_cutoff {params.rulesigma} \
            --lowest_freq {params.q1} \
            --ploidy {params.ploidy}
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
        windowsize=str(config['fixed']['isweep']['window_size']),
        windowstep=str(config['fixed']['isweep']['window_step']),
        qrng=str(config['fixed']['isweep']['num_in_quantile_sum']),
        maxspace=str(config['fixed']['isweep']['max_spacing']),
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
