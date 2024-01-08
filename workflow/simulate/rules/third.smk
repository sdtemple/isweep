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
        soft=str(config['CHANGE']['FOLDERS']['SOFTWARE']),
        prog=str(config['CHANGE']['PROGRAMS']['FILTER']),
        script=str(config['CHANGE']['FOLDERS']['TERMINALSCRIPTS'])+'/lines.py',
        mlecutoff=str(config['FIXED']['ISWEEP']['MLECUTOFF']),
    resources:
        mem_gb='{config[CHANGE][CLUSTER][LARGEMEM]}'
    shell:
        """
        thecenter=$(python {params.script} {input.locus} 1 2)
        zcat {input.ibd} | \
            java -Xmx{config[CHANGE][CLUSTER][LARGEMEM]}g -jar {params.soft}/{params.prog} \
            "I" 6 0.00 ${{thecenter}} | \
            java -Xmx{config[CHANGE][CLUSTER][LARGEMEM]}g -jar {params.soft}/{params.prog} \
            "I" 7 ${{thecenter}} 10000000000 | \
            java -Xmx{config[CHANGE][CLUSTER][LARGEMEM]}g -jar {params.soft}/{params.prog} \
            "I" -8 0.00 {params.mlecutoff} | \
            gzip > {output.ibd}
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
        soft=str(config['CHANGE']['FOLDERS']['SOFTWARE']),
        prog=str(config['CHANGE']['PROGRAMS']['FILTER']),
        script=str(config['CHANGE']['FOLDERS']['TERMINALSCRIPTS'])+'/lines.py',
        mlecutoff=str(config['FIXED']['ISWEEP']['MLECUTOFF']),
    resources:
        mem_gb='{config[CHANGE][CLUSTER][LARGEMEM]}'
    shell:
        """
        thecenter=$(python {params.script} {input.locus} 1 2)
        zcat {input.ibd} | \
            java -Xmx{config[CHANGE][CLUSTER][LARGEMEM]}g -jar {params.soft}/{params.prog} \
            "I" 6 0.00 ${{thecenter}} | \
            java -Xmx{config[CHANGE][CLUSTER][LARGEMEM]}g -jar {params.soft}/{params.prog} \
            "I" 7 ${{thecenter}} 10000000000 | \
            java -Xmx{config[CHANGE][CLUSTER][LARGEMEM]}g -jar {params.soft}/{params.prog} \
            "I" -8 0.00 {params.mlecutoff} | \
            gzip > {output.ibd}
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

