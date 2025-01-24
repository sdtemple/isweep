#s Finding excess IBD groups

localrules: gini_impurity

# subset vcf to region of interest
rule second_region: # focus vcf on region of interest
    input:
        locus='{cohort}/{hit}/locus.txt',
        focus='{cohort}/{hit}/first.pos.txt',
        vcf='{cohort}/{hit}/first.focused.vcf.gz',
        subsample=cohort+"/subsample.txt",
    output:
        subvcf='{cohort}/{hit}/second.focused.vcf.gz',
    params:
        qmaf=maf,
        chrpre=str(config['change']['files']['chromosome_prefix']),
        pm=str(config['fixed']['isweep']['base_pair_buffer']),
    shell:
        """
        chr=$(python ../../scripts/utilities/lines.py {input.locus} 2 2)
        center=$(python ../../scripts/utilities/lines.py {input.focus} 1 2)
        left=$(python -c "out = $center - {params.pm} ; print(max(out,2))")
        right=$(python -c "out = $center + {params.pm} ; print(out)")
        tabix -fp vcf {input.vcf}
        bcftools view {input.vcf} -r {params.chrpre}${{chr}}:${{left}}-${{right}} -Ob | \
            bcftools view -S {input.subsample} -Ob | \
            bcftools view -q {params.qmaf}:nonmajor -Oz -o {output.subvcf}
        """

### filter ibd file ###

rule second_filt:
    input:
        ibd='{cohort}/{hit}/first.ibd.gz',
        locus='{cohort}/{hit}/first.pos.txt',
    output:
        ibd='{cohort}/{hit}/second.filt.ibd.gz',
    shell:
        """
        thecenter=$(python ../../scripts/utilities/lines.py {input.locus} 1 2)
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

rule second_rank:
    input:
        short='{cohort}/{hit}/second.filt.ibd.gz',
        vcf='{cohort}/{hit}/second.focused.vcf.gz',
    output:
        fileout='{cohort}/{hit}/second.ranks.tsv.gz',
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

### write outliers ###

rule second_outlier:
    input:
        short='{cohort}/{hit}/second.filt.ibd.gz',
    output:
        fileout='{cohort}/{hit}/second.outliers.txt',
        out1='{cohort}/{hit}/outlier1.txt',
    params:
        diameter=diameter,
        rulesigma=group_cutoff,
    shell:
        """
        python ../../scripts/model/outliers.py \
            --input_ibd_file {input.short} \
            --output_folder {wildcards.cohort}/{wildcards.hit} \
            --graph_diameter {params.diameter} \
            --group_cutoff {params.rulesigma}
        touch {output.fileout}
        touch {output.out1}
        """

### write entropy ###

rule gini_impurity:
	input:
		filein='{cohort}/{hit}/outlier1.txt',
	output:
		fileout='{cohort}/{hit}/ibd.gini.tsv',
	params:
		samplesizep=str(samplesize_ploidy),
	shell:
		"""
		python ../../scripts/model/ibd-gini-entropy.py \
			--input_folder {wildcards.cohort}/{wildcards.hit} \
			--output_file {output.fileout} \
			--sample_size {params.samplesizep}
		"""
