wildcard_constraints:
	SIMNAME = '\w+',

##### haplotype analysis #####

rule third_haplotypes:
    input:
        rankin='{macro}/{micro}/{seed}/second.ranks.tsv.gz',
        ibdwin='{macro}/{micro}/{seed}/scan.chr1.windowed.tsv.gz',
    output:
        lociout='{macro}/{micro}/{seed}/third.pos.txt',
        freqout='{macro}/{micro}/{seed}/third.freq.txt',
        happng='{macro}/{micro}/{seed}/third.hap.png',
        snppng='{macro}/{micro}/{seed}/third.snp.png',
    params:
        windowsize=str(config['FIXED']['ISWEEP']['WINSIZE']),
        windowstep=str(config['FIXED']['ISWEEP']['WINSTEP']),
        windowcol=str(config['FIXED']['ISWEEP']['WINCOL']),
        freqsize=str(config['FIXED']['ISWEEP']['FREQSIZE']),
        freqstep=str(config['FIXED']['ISWEEP']['FREQSTEP']),
        freqcol=str(config['FIXED']['ISWEEP']['FREQCOL']),
        scorecol=str(config['FIXED']['ISWEEP']['SCORECOL']),
        numsnp=str(config['FIXED']['ISWEEP']['NUMSNP']),
        scripts=str(config['CHANGE']['FOLDERS']['TERMINALSCRIPTS']),
    shell:
        """
        python {params.scripts}/haplotypes.py \
            {input.rankin} \
            {input.ibdwin} \
            {wildcards.macro}/{wildcards.micro}/{wildcards.seed} \
            0 \
            BPWINDOW \
            CMWINDOW \
            {params.freqsize} \
            {params.freqstep} \
            {params.freqcol} \
            {params.windowsize} \
            {params.windowstep} \
            {params.windowcol} \
            {params.scorecol} \
            {params.numsnp}
        """

rule third_ibd:
    input:
        ibd='{macro}/{micro}/{seed}/scan.chr1.ibd.gz',
        locus='{macro}/{micro}/{seed}/third.pos.txt',
    output:
        ibd='{macro}/{micro}/{seed}/third.chr1.ibd.gz',
    params:
        soft=str(config['CHANGE']['FOLDERS']['SOFTWARE']),
        prog=str(config['CHANGE']['PROGRAMS']['FILTER']),
        script=str(config['CHANGE']['FOLDERS']['TERMINALSCRIPTS'])+'/first-line.py',
    resources:
        mem_gb='{config[CHANGE][CLUSTER][LARGEMEM]}'
    shell:
        """
        thecenter=$(python {params.script} {input.locus})
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
        freq='{macro}/{micro}/{seed}/third.freq.txt',
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
        freqest=$(python {params.scripts}/first-line.py {input.freq})
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
