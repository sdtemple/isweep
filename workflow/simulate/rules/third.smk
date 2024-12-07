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
        scripts=str(config['CHANGE']['FOLDERS']['TERMINALSCRIPTS']),
    shell:
        """
        python {params.scripts}/haplotypes.py \
            {input.rankin} \
            {wildcards.macro}/{wildcards.micro}/{wildcards.seed} \
            0 \
            1 \
            -1 \
            {params.freqsize} \
            {params.freqstep} \
            {params.windowsize} \
            {params.windowstep} \
            {params.numsnp} \
            {params.lowbnd}
        """

rule third_hap_ibd:
    input:
        ibd='{macro}/{micro}/{seed}/scan.chr1.ibd.gz',
        locus='{macro}/{micro}/{seed}/third.best.hap.txt',
    output:
        ibd='{macro}/{micro}/{seed}/third.chr1.hap.ibd.gz',
    params:
        scripts=str(config['CHANGE']['FOLDERS']['TERMINALSCRIPTS']),
        mlecutoff=str(config['FIXED']['ISWEEP']['MLECUTOFF']),
    resources:
        mem_gb='{config[CHANGE][CLUSTER][LARGEMEM]}'
    shell:
        """
        thecenter=$(python {params.scripts}/lines.py {input.best} 1 2)
        python {params.scripts}/filter-lines.py \
            {input.ibd} \
            {wildcards.cohort}/{wildcards.hit}/intermediate.ibd.gz \
            --column_index 6 \
            --upper_bound $thecenter \
            --complement 0
        python {params.scripts}/filter-lines.py \
            {wildcards.cohort}/{wildcards.hit}/intermediate.ibd.gz \
            {wildcards.cohort}/{wildcards.hit}/intermediate2.ibd.gz \
            --column_index 7 \
            --lower_bound $thecenter \
            --upper_bound 10000000000 \
            --complement 0
        python {params.scripts}/filter-lines.py \
            {wildcards.cohort}/{wildcards.hit}/intermediate2.ibd.gz \
            {output.ibd} \
            --upper_bound {params.mlecutoff}
        rm {wildcards.cohort}/{wildcards.hit}/intermediate.ibd.gz
        rm {wildcards.cohort}/{wildcards.hit}/intermediate2.ibd.gz
        """

rule third_hap_infer:
    input:
        long='{macro}/{micro}/{seed}/third.chr1.hap.ibd.gz',
        freq='{macro}/{micro}/{seed}/third.best.hap.txt',
    output:
        fileout='{macro}/{micro}/{seed}/results.hap.tsv',
    params:
        scripts=str(config['CHANGE']['FOLDERS']['TERMINALSCRIPTS']),
        nboot=str(config['FIXED']['ISWEEP']['NBOOT']),
        mlecutoff=str(config['FIXED']['ISWEEP']['MLECUTOFF']),
        n=str(config['CHANGE']['SIMULATE']['SAMPSIZE']),
        ploidy=str(config['FIXED']['SIMULATE']['PLOIDY']),
        effdemo=str(config['CHANGE']['SIMULATE']['iNe']),
    shell:
        """
        ibdest=$(zcat {input.long} | wc -l)
        freqest=$(python {params.scripts}/lines.py {input.freq} 2 2)
        python {params.scripts}/estimate.py \
            {output.fileout} \
            ${{ibdest}} \
            ${{freqest}} \
            {params.nboot} \
            {params.mlecutoff} \
            {params.n} \
            {params.effdemo} \
            {params.ploidy}
        """

##### snp analysis #####

rule third_snp:
    input:
        rankin='{macro}/{micro}/{seed}/second.ranks.tsv.gz',
        outlied='{macro}/{micro}/{seed}/second.outliers.txt',
    output:
        lociout='{macro}/{micro}/{seed}/third.best.snp.txt',
    params:
        scripts=str(config['CHANGE']['FOLDERS']['TERMINALSCRIPTS']),
        lowbnd=str(config['FIXED']['ISWEEP']['MINAAF']),
    shell:
        """
        python {params.scripts}/snp.py \
            {input.rankin} \
            {output.lociout} \
            {params.lowbnd}
        """

rule third_snp_ibd:
    input:
        ibd='{macro}/{micro}/{seed}/scan.chr1.ibd.gz',
        locus='{macro}/{micro}/{seed}/third.best.snp.txt',
    output:
        ibd='{macro}/{micro}/{seed}/third.chr1.snp.ibd.gz',
    params:
        scripts=str(config['CHANGE']['FOLDERS']['TERMINALSCRIPTS']),
        mlecutoff=str(config['FIXED']['ISWEEP']['MLECUTOFF']),
    resources:
        mem_gb='{config[CHANGE][CLUSTER][LARGEMEM]}'
    shell:
        """
        thecenter=$(python {params.scripts}/lines.py {input.best} 1 2)
        python {params.scripts}/filter-lines.py \
            {input.ibd} \
            {wildcards.cohort}/{wildcards.hit}/intermediate.ibd.gz \
            --column_index 6 \
            --upper_bound $thecenter \
            --complement 0
        python {params.scripts}/filter-lines.py \
            {wildcards.cohort}/{wildcards.hit}/intermediate.ibd.gz \
            {wildcards.cohort}/{wildcards.hit}/intermediate2.ibd.gz \
            --column_index 7 \
            --lower_bound $thecenter \
            --upper_bound 10000000000 \
            --complement 0
        python {params.scripts}/filter-lines.py \
            {wildcards.cohort}/{wildcards.hit}/intermediate2.ibd.gz \
            {output.ibd} \
            --upper_bound {params.mlecutoff}
        rm {wildcards.cohort}/{wildcards.hit}/intermediate.ibd.gz
        rm {wildcards.cohort}/{wildcards.hit}/intermediate2.ibd.gz
        """

rule third_snp_infer:
    input:
        long='{macro}/{micro}/{seed}/third.chr1.snp.ibd.gz',
        freq='{macro}/{micro}/{seed}/third.best.snp.txt',
    output:
        fileout='{macro}/{micro}/{seed}/results.snp.tsv',
    params:
        scripts=str(config['CHANGE']['FOLDERS']['TERMINALSCRIPTS']),
        nboot=str(config['FIXED']['ISWEEP']['NBOOT']),
        mlecutoff=str(config['FIXED']['ISWEEP']['MLECUTOFF']),
        n=str(config['CHANGE']['SIMULATE']['SAMPSIZE']),
        ploidy=str(config['FIXED']['SIMULATE']['PLOIDY']),
        effdemo=str(config['CHANGE']['SIMULATE']['iNe']),
    shell:
        """
        ibdest=$(zcat {input.long} | wc -l)
        freqest=$(python {params.scripts}/lines.py {input.freq} 2 2)
        python {params.scripts}/estimate.py \
            {output.fileout} \
            ${{ibdest}} \
            ${{freqest}} \
            {params.nboot} \
            {params.mlecutoff} \
            {params.n} \
            {params.effdemo} \
            {params.ploidy}
        """

