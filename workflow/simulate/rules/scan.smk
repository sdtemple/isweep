# looking for regions of excess ibd

# inputs
n=int(float(config['CHANGE']['SIMULATE']['SAMPSIZE']))
ploidy=int(float(config['FIXED']['SIMULATE']['PLOIDY']))
maf1=float(config['FIXED']['CANDHAPIBD']['MINMAF'])
mac1=int(ploidy*n*maf1)

### ibd-ends from candidate segments ###

rule hapibd:
    input:
        vcf='{macro}/{micro}/{seed}/large.chr1.vcf.gz',
        map=str(config['CHANGE']["FOLDERS"]["MACRO"]) + '/uniform.map',
    output:
        ibd='{macro}/{micro}/{seed}/large.chr1.hapibd.candidate.ibd.gz',
    params:
        minmac=str(mac1),
        out='{macro}/{micro}/{seed}/large.chr1.hapibd.candidate',
        soft=str(config['CHANGE']['FOLDERS']['SOFTWARE']),
        prog=str(config['CHANGE']['PROGRAMS']['HAPIBD']),
        minsee=str(config['FIXED']['CANDHAPIBD']['MINSEED']),
        minext=str(config['FIXED']['CANDHAPIBD']['MINEXT']),
        minout=str(config['FIXED']['CANDHAPIBD']['MINOUT']),
    resources:
        mem_gb='{config[CHANGE][CLUSTER][LARGEMEM]}'
    shell:
        """
        java -Xmx{config[CHANGE][CLUSTER][LARGEMEM]}g -jar {params.soft}/{params.prog} \
            gt={input.vcf} \
            map={input.map} \
            out={params.out} \
            min-seed={params.minsee} \
            min-extend={params.minext} \
            min-output={params.minout} \
            min-mac={params.minmac}
        """

rule ibdends:
    input:
        vcf='{macro}/{micro}/{seed}/large.chr1.vcf.gz',
        map=str(config['CHANGE']["FOLDERS"]["MACRO"]) + '/uniform.map',
        ibd='{macro}/{micro}/{seed}/large.chr1.hapibd.candidate.ibd.gz'
    params:
        out='{macro}/{micro}/{seed}/large.chr1.ibdends',
        soft=str(config['CHANGE']['FOLDERS']['SOFTWARE']),
        prog=str(config['CHANGE']['PROGRAMS']['IBDENDS']),
        maf=str(config['FIXED']['IBDENDS']['MINMAF']),
        qua=str(config['FIXED']['IBDENDS']['QUANTILES']),
        sam=str(config['FIXED']['IBDENDS']['NSAMPLES']),
        err=str(config['FIXED']['IBDENDS']['ERRRATE']),
    output:
        ibd='{macro}/{micro}/{seed}/large.chr1.ibdends.ibd.gz',
    resources:
        mem_gb='{config[CHANGE][CLUSTER][LARGEMEM]}'
    shell:
        """
        java -Xmx{config[CHANGE][CLUSTER][LARGEMEM]}g -jar {params.soft}/{params.prog} \
            gt={input.vcf} \
            ibd={input.ibd} \
            map={input.map} \
            out={params.out} \
            min-maf={params.maf} \
            quantiles={params.qua} \
            nsamples={params.sam} \
            err={params.err}
        """

rule format_ibdends:
    input:
        bigibd='{macro}/{micro}/{seed}/large.chr1.ibdends.ibd.gz',
    output:
        cutibd='{macro}/{micro}/{seed}/large.chr1.ibdends.modified.ibd.gz',
    shell:
        'zcat {input.bigibd} | tail -n +2 | cut -f 1-5,10,11,12 | gzip > {output.cutibd}'

rule filter_ibdends:
    input:
        ibd='{macro}/{micro}/{seed}/large.chr1.ibdends.modified.ibd.gz',
    output:
        fipass='{macro}/{micro}/{seed}/scan.chr1.ibd.gz',
    params:
        soft=str(config['CHANGE']['FOLDERS']['SOFTWARE']),
        prog=str(config['CHANGE']['PROGRAMS']['FILTER']),
        momcut=str(config['FIXED']['ISWEEP']['SCANCUTOFF']),
    resources:
        mem_gb='{config[CHANGE][CLUSTER][LARGEMEM]}'
    shell:
        """
        zcat {input.ibd} | \
            java -Xmx{config[CHANGE][CLUSTER][LARGEMEM]}g -jar {params.soft}/{params.prog} "I" -8 0.00 {params.momcut} | \
            gzip > {output.fipass}
        """

rule count_ibdends:
    input:
        filein='{macro}/{micro}/{seed}/scan.chr1.ibd.gz',
        mapin=str(config['CHANGE']["FOLDERS"]["MACRO"]) + '/uniform.map',
    output:
        fileout='{macro}/{micro}/{seed}/scan.chr1.windowed.tsv.gz',
    params:
        scripts=str(config['CHANGE']['FOLDERS']['TERMINALSCRIPTS']),
        by=str(config['FIXED']['ISWEEP']['BY']),
        ge=str(config['FIXED']['ISWEEP']['GENOMEEND']),
        tc=str(config['FIXED']['ISWEEP']['TELOCUT']),
    shell:
        """
        python {params.scripts}/ibd-windowed.py \
            {input.filein} \
            {output.fileout} \
            {input.mapin} \
            {params.by} \
            {params.ge} \
            {params.tc}
        """
