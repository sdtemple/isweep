### finding ibd groups

# some inputs, string managements, count sample size
cohort=str(config['CHANGE']['FOLDERS']['STUDY'])
ploidy=2
# ploidy=int(float(config['FIXED']['HAPIBD']['PLOIDY']))
maf=float(config['FIXED']['ISWEEP']['MINAAF'])

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
        chrpre=str(config['CHANGE']['ISWEEP']['CHRPRE']),
        pm=str(config['FIXED']['ISWEEP']['PM']),
    shell: # if chromosome is huge (greater than 10000 Mb), may need to modify the third pipe
        """
        chr=$(python ../../scripts/lines.py {input.locus} 2 2)
        center=$(python ../../scripts/lines.py {input.focus} 1 2)
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
        thecenter=$(python ../../scripts/lines.py {input.locus} 1 2)
        python ../../scripts/filter-lines.py \
            --input_file {input.ibd} \
            --output_file {wildcards.cohort}/{wildcards.hit}/intermediate.ibd.gz \
            --column_index 6 \
            --upper_bound $thecenter \
            --complement 0
        python ../../scripts/filter-lines.py \
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
        diameter=str(config['FIXED']['ISWEEP']['DIAMETER']),
        q1=str(config['FIXED']['ISWEEP']['MINAAF']),
        rulesigma=str(config['FIXED']['ISWEEP']['GROUPCUTOFF']),
    shell:
        """
        python ../../scripts/rank.py \
            --ibd_file {input.short} \
            --vcf {input.vcf} \
            --file_out {output.fileout} \
            --graph_diameter {params.diameter} \
            --group_cutoff {params.rulesigma} \
            --lowest_freq {params.q1} \
        """

### write outliers ###

rule second_outlier:
    input:
        short='{cohort}/{hit}/second.filt.ibd.gz',
    output:
        fileout='{cohort}/{hit}/second.outliers.txt',
        out1='{cohort}/{hit}/outlier1.txt',
    params:
        diameter=str(config['FIXED']['ISWEEP']['DIAMETER']),
        rulesigma=str(config['FIXED']['ISWEEP']['GROUPCUTOFF']),
    shell:
        """
        python ../../scripts/outliers.py \
            --ibd_file {input.short} \
            --folder_out {wildcards.cohort}/{wildcards.hit} \
            --graph_diameter {params.diameter} \
            --group_cutoff {params.rulesigma}
        touch {output.fileout}
        touch {output.out1}
        """