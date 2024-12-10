##### haplotype analysis #####

### estimating selection coefficients

rule third_hap:
    input:
        rankin='{macro}/{micro}/{seed}/second.ranks.tsv.gz',
        outlied='{macro}/{micro}/{seed}/second.outliers.txt',
    output:
        lociout='{macro}/{micro}/{seed}/third.best.hap.txt',
        happng='{macro}/{micro}/{seed}/third.hap.png',
        snppng='{macro}/{micro}/{seed}/third.snp.png',
    params:
        windowsize=str(config['FIXED']['ISWEEP']['WINSIZE']),
        windowstep=str(config['FIXED']['ISWEEP']['WINSTEP']),
        freqsize=str(config['FIXED']['ISWEEP']['FREQSIZE']),
        freqstep=str(config['FIXED']['ISWEEP']['FREQSTEP']),
        numsnp=str(config['FIXED']['ISWEEP']['NUMSNP']),
        lowbnd=str(config['FIXED']['ISWEEP']['MINAAF']),
    shell:
        """
        python ../../scripts/model/haplotypes.py \
            --input_snp_file {input.rankin} \
            --output_folder {wildcards.macro}/{wildcards.micro}/{wildcards.seed} \
            --window_index 0 \
            --freq_index 1 \
            --score_index -1 \
            --freq_size {params.freqsize} \
            --freq_step {params.freqstep} \
            --window_size {params.windowsize} \
            --window_step {params.windowstep} \
            --num_snp {params.numsnp} \
            --lowest_freq {params.lowbnd}
        """

rule third_hap_ibd:
    input:
        ibd='{macro}/{micro}/{seed}/scan.chr1.ibd.gz',
        locus='{macro}/{micro}/{seed}/third.best.hap.txt',
    output:
        ibd='{macro}/{micro}/{seed}/third.chr1.hap.ibd.gz',
    params:
        mlecutoff=str(config['FIXED']['ISWEEP']['MLECUTOFF']),
    shell:
        """
        thecenter=$(python ../../scripts/lines.py {input.locus} 1 2)
        python ../../scripts/utilities/filter-lines.py \
            --input_file {input.ibd} \
            --output_file {wildcards.macro}/{wildcards.micro}/{wildcards.seed}/intermediate.ibd.gz \
            --column_index 6 \
            --upper_bound $thecenter \
            --complement 0
        python ../../scripts/utilities/filter-lines.py \
            --input_file {wildcards.macro}/{wildcards.micro}/{wildcards.seed}/intermediate.ibd.gz \
            --output_file {wildcards.macro}/{wildcards.micro}/{wildcards.seed}/intermediate2.ibd.gz \
            --column_index 7 \
            --lower_bound $thecenter \
            --upper_bound 10000000000 \
            --complement 0
        python ../../scripts/utilities/filter-lines.py \
            --input_file {wildcards.macro}/{wildcards.micro}/{wildcards.seed}/intermediate2.ibd.gz \
            --output_file {output.ibd} \
            --upper_bound {params.mlecutoff}
        rm {wildcards.macro}/{wildcards.micro}/{wildcards.seed}/intermediate.ibd.gz
        rm {wildcards.macro}/{wildcards.micro}/{wildcards.seed}/intermediate2.ibd.gz
        """

rule third_hap_infer:
    input:
        long='{macro}/{micro}/{seed}/third.chr1.hap.ibd.gz',
        freq='{macro}/{micro}/{seed}/third.best.hap.txt',
    output:
        fileout='{macro}/{micro}/{seed}/results.hap.tsv',
    params:
        nboot=str(config['FIXED']['ISWEEP']['NBOOT']),
        mlecutoff=str(config['FIXED']['ISWEEP']['MLECUTOFF']),
        n=str(config['CHANGE']['SIMULATE']['SAMPSIZE']),
        ploidy=str(config['FIXED']['SIMULATE']['PLOIDY']),
        longNe=str(config['CHANGE']['SIMULATE']['iNe']),
        model='a',
        alpha=0.05,
    shell:
        """
        ibdest=$(zcat {input.long} | wc -l)
        freqest=$(python ../../scripts/utilities/lines.py {input.freq} 2 2)
        python ../../scripts/model/estimate-norm.py \
            --output_file {output.fileout} \
            --ibd_count ${{ibdest}} \
            --p_est ${{freqest}} \
            --num_bootstraps {params.nboot} \
            --ibd_cutoff {params.mlecutoff} \
            --sample_size {params.n} \
            --Ne_est {params.longNe} \
            --model {params.model} \
            --alpha {params.alpha} \
            --ploidy {params.ploidy} \
        """

##### snp analysis #####

rule third_snp:
    input:
        rankin='{macro}/{micro}/{seed}/second.ranks.tsv.gz',
        outlied='{macro}/{micro}/{seed}/second.outliers.txt',
    output:
        lociout='{macro}/{micro}/{seed}/third.best.snp.txt',
    params:
        lowbnd=str(config['FIXED']['ISWEEP']['MINAAF']),
    shell:
        """
        python ../../scripts/model/snp.py \
            --input_snp_file {input.rankin} \
            --output_file {output.lociout} \
            --lowest_freq {params.lowbnd}
        """

rule third_snp_ibd:
    input:
        ibd='{macro}/{micro}/{seed}/scan.chr1.ibd.gz',
        locus='{macro}/{micro}/{seed}/third.best.snp.txt',
    output:
        ibd='{macro}/{micro}/{seed}/third.chr1.snp.ibd.gz',
    params:
        mlecutoff=str(config['FIXED']['ISWEEP']['MLECUTOFF']),
    resources:
        mem_gb='{config[CHANGE][CLUSTER][LARGEMEM]}'
    shell:
        """
        thecenter=$(python ../../scripts/utilities/lines.py {input.locus} 1 2)
        python ../../scripts/utilities/filter-lines.py \
            --input_file {input.ibd} \
            --output_file {wildcards.macro}/{wildcards.micro}/{wildcards.seed}/intermediate3.ibd.gz \
            --column_index 6 \
            --upper_bound $thecenter \
            --complement 0
        python ../../scripts/utilities/filter-lines.py \
            --input_file {wildcards.macro}/{wildcards.micro}/{wildcards.seed}/intermediate3.ibd.gz \
            --output_file {wildcards.macro}/{wildcards.micro}/{wildcards.seed}/intermediate4.ibd.gz \
            --column_index 7 \
            --lower_bound $thecenter \
            --upper_bound 10000000000 \
            --complement 0
        python ../../scripts/utilities/filter-lines.py \
            --input_file {wildcards.macro}/{wildcards.micro}/{wildcards.seed}/intermediate4.ibd.gz \
            --output_file {output.ibd} \
            --upper_bound {params.mlecutoff}
        rm {wildcards.macro}/{wildcards.micro}/{wildcards.seed}/intermediate3.ibd.gz
        rm {wildcards.macro}/{wildcards.micro}/{wildcards.seed}/intermediate4.ibd.gz
        """

rule third_snp_infer:
    input:
        long='{macro}/{micro}/{seed}/third.chr1.snp.ibd.gz',
        freq='{macro}/{micro}/{seed}/third.best.snp.txt',
    output:
        fileout='{macro}/{micro}/{seed}/results.snp.tsv',
    params:
        nboot=str(config['FIXED']['ISWEEP']['NBOOT']),
        mlecutoff=str(config['FIXED']['ISWEEP']['MLECUTOFF']),
        n=str(config['CHANGE']['SIMULATE']['SAMPSIZE']),
        ploidy=str(config['FIXED']['SIMULATE']['PLOIDY']),
        longNe=str(config['CHANGE']['SIMULATE']['iNe']),
        model='a',
        alpha=0.05,
    shell:
        """
        ibdest=$(zcat {input.long} | wc -l)
        freqest=$(python ../../scripts/utilities/lines.py {input.freq} 2 2)
        python ../../scripts/model/estimate-norm.py \
            --output_file {output.fileout} \
            --ibd_count ${{ibdest}} \
            --p_est ${{freqest}} \
            --num_bootstraps {params.nboot} \
            --ibd_cutoff {params.mlecutoff} \
            --sample_size {params.n} \
            --Ne_est {params.longNe} \
            --model {params.model} \
            --alpha {params.alpha} \
            --ploidy {params.ploidy} \
        """

