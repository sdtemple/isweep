wildcard_constraints:
	SIMNAME = '\w+',

##### haplotype analysis #####

rule third_haplotypes:
    input:
        rankin='{macro}/{micro}/{seed}/second.ranks.tsv.gz',
    output:
        lociout='{macro}/{micro}/{seed}/third.best.txt',
        happng='{macro}/{micro}/{seed}/third.hap.png',
        snppng='{macro}/{micro}/{seed}/third.snp.png',
    params:
        windowsize=str(config['FIXED']['ISWEEP']['WINSIZE']),
        windowstep=str(config['FIXED']['ISWEEP']['WINSTEP']),
        freqsize=str(config['FIXED']['ISWEEP']['FREQSIZE']),
        freqstep=str(config['FIXED']['ISWEEP']['FREQSTEP']),
        numsnp=str(config['FIXED']['ISWEEP']['NUMSNP']),
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
            {params.numsnp}
        """

rule third_ibd:
    input:
        ibd='{macro}/{micro}/{seed}/scan.chr1.ibd.gz',
        locus='{macro}/{micro}/{seed}/third.best.txt',
    output:
        ibd='{macro}/{micro}/{seed}/third.chr1.ibd.gz',
    params:
        soft=str(config['CHANGE']['FOLDERS']['SOFTWARE']),
        prog=str(config['CHANGE']['PROGRAMS']['FILTER']),
        script=str(config['CHANGE']['FOLDERS']['TERMINALSCRIPTS'])+'/lines.py',
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
            gzip > {output.ibd}
        """

rule third_infer:
    input:
        long='{macro}/{micro}/{seed}/third.chr1.ibd.gz',
        freq='{macro}/{micro}/{seed}/third.best.txt',
    output:
        fileout='{macro}/{micro}/{seed}/results.tsv',
    params:
        scripts=str(config['CHANGE']['FOLDERS']['TERMINALSCRIPTS']),
        nboot=str(config['FIXED']['ISWEEP']['NBOOT']),
        cm=str(config['FIXED']['ISWEEP']['MOMCUTOFF']),
        n=str(config['CHANGE']['SIMULATE']['SAMPSIZE']),
        ploidy=str(config['FIXED']['HAPIBD']['PLOIDY']),
        effdemo=str(config['CHANGE']['SIMULATE']['iNe']),
    shell:
        """
        freqest=$(python {params.scripts}/lines.py {input.freq} 2 2)
        python {params.scripts}/estimate.py \
            {input.long} \
            {output.fileout} \
            ${{freqest}} \
            {params.nboot} \
            {params.cm} \
            {params.n} \
            {params.effdemo} \
            {params.ploidy}
        """
