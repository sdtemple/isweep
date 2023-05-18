# conduct ibd calls and scan for isweep
# seth temple, sdtemple@uw.edu
# may 3, 2023

# inputs, string management, count sample size, make mac
subsamplefile=str(config['CHANGE']['ISWEEP']['SUBSAMPLE'])
macro=str(config['CHANGE']['FOLDERS']['STUDY'])
samplesize=0
with open(subsamplefile,'r') as f:
# with open(macro+'/'+subsamplefile,'r') as f:
    for line in f:
        samplesize+=1
# samplesize=str(samplesize)
ploidy=int(float(config['FIXED']['CANDHAPIBD']['PLOIDY']))
maf1=float(config['FIXED']['CANDHAPIBD']['MINMAF'])
mac1=int(ploidy*samplesize*maf1)
macro=str(config['CHANGE']['FOLDERS']['STUDY'])
low=int(float(str(config['CHANGE']['ISWEEP']['CHRLOW'])))
high=int(float(str(config['CHANGE']['ISWEEP']['CHRHIGH'])))
vcffolder=str(config['CHANGE']['EXISTING']['VCFS'])
vcfpre=str(config['CHANGE']['EXISTING']['VCFPRE'])
vcfsuf=str(config['CHANGE']['EXISTING']['VCFSUF'])
concatmap=macro+'/maps/chr'+str(low)+'-'+str(high)+'.map'

rule hapibd: # candidate segments from hap-ibd.jar
    input:
        vcf=vcffolder + '/' + vcfpre + '{num}' + vcfsuf,
        map='{cohort}/maps/chr{num}.map',
        excludesamples='{cohort}/excludesamples.txt',
    output:
        ibd='{cohort}/ibdsegs/hapibd/chr{num}.ibd.gz',
    params:
        minmac=str(mac1),
        out='{cohort}/ibdsegs/hapibd/chr{num}',
        soft=str(config['CHANGE']['FOLDERS']['SOFTWARE']),
        prog=str(config['CHANGE']['PROGRAMS']['HAPIBD']),
        xmx=str(config['CHANGE']['ISWEEP']['XMXMEM']),
        minsee=str(config['FIXED']['CANDHAPIBD']['MINSEED']),
        minext=str(config['FIXED']['CANDHAPIBD']['MINEXT']),
        minout=str(config['FIXED']['CANDHAPIBD']['MINOUT']),
    resources:
        mem_gb=100
    shell:
        """
        java -Xmx{params.xmx}g -jar {params.soft}/{params.prog} \
            gt={input.vcf} \
            map={input.map} \
            out={params.out} \
            min-seed={params.minsee} \
            min-extend={params.minext} \
            min-output={params.minout} \
            min-mac={params.minmac} \
            excludesamples={input.excludesamples}
        """

rule ibdends: # ibd-ends.jar
    input:
        vcf=vcffolder + '/' + vcfpre + '{num}' + vcfsuf,
        map='{cohort}/maps/chr{num}.map',
        ibd='{cohort}/ibdsegs/hapibd/chr{num}.ibd.gz',
        subsample='{cohort}/excludesamples.txt',
    output:
        ibd='{cohort}/ibdsegs/ibdends/chr{num}.ibd.gz',
    params:
        out='{cohort}/ibdsegs/ibdends/chr{num}',
        soft=str(config['CHANGE']['FOLDERS']['SOFTWARE']),
        prog=str(config['CHANGE']['PROGRAMS']['IBDENDS']),
        xmx=str(config['CHANGE']['ISWEEP']['XMXMEM']),
        maf=str(config['FIXED']['IBDENDS']['MINMAF']),
        qua=str(config['FIXED']['IBDENDS']['QUANTILES']),
        sam=str(config['FIXED']['IBDENDS']['NSAMPLES']),
        err=str(config['FIXED']['IBDENDS']['ERRRATE']),
    resources:
        mem_gb=100
    shell:
        """
        java -Xmx{params.xmx}g -jar {params.soft}/{params.prog} \
            gt={input.vcf} \
            ibd={input.ibd} \
            map={input.map} \
            out={params.out} \
            min-maf={params.maf} \
            quantiles={params.qua} \
            nsamples={params.sam} \
            err={params.err} \
            excludesamples={input.subsample}
        """

rule format_ibdends: # reformatting
    input:
        bigibd='{cohort}/ibdsegs/ibdends/chr{num}.ibd.gz',
    output:
        cutibd='{cohort}/ibdsegs/ibdends/modified/chr{num}.ibd.gz',
    resources:
        mem_gb=10
    shell:
        'zcat {input.bigibd} | tail -n +2 | cut -f 1-5,10,11,12 | gzip > {output.cutibd}'

rule filter_ibdends_scan: # applying cutoffs
    input:
        ibd='{cohort}/ibdsegs/ibdends/modified/chr{num}.ibd.gz',
    output:
        fipass='{cohort}/ibdsegs/ibdends/modified/scan/chr{num}.ibd.gz',
    params:
        soft=str(config['CHANGE']['FOLDERS']['SOFTWARE']),
        prog=str(config['CHANGE']['PROGRAMS']['FILTER']),
        xmx=str(config['CHANGE']['ISWEEP']['XMXMEM']),
        scancut=str(config['FIXED']['ISWEEP']['SCANCUTOFF']),
    resources:
        mem_gb=10
    shell:
        """
        zcat {input.ibd} | \
            java -Xmx{params.xmx}g -jar {params.soft}/{params.prog} "I" -8 0.00 {params.scancut} | \
            gzip > {output.fipass}
        """

rule count_ibdends_scan: # computing counts over windows
    input:
        filein='{cohort}/ibdsegs/ibdends/modified/scan/chr{num}.ibd.gz',
        mapin='{cohort}/maps/chr{num}.map',
    output:
        fileout='{cohort}/ibdsegs/ibdends/modified/scan/chr{num}.ibd.windowed.tsv.gz',
    params:
        scripts=str(config['CHANGE']['FOLDERS']['TERMINALSCRIPTS']),
        by=str(config['FIXED']['ISWEEP']['BY']),
        ge=str(config['FIXED']['ISWEEP']['GENOMEEND']),
        tc=str(config['FIXED']['ISWEEP']['TELOCUT']),
    resources:
        mem_gb=10
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

rule filter_ibdends_mom: # applying cutoffs
    input:
        ibd='{cohort}/ibdsegs/ibdends/modified/chr{num}.ibd.gz',
    output:
        fipass='{cohort}/ibdsegs/ibdends/modified/mom/chr{num}.ibd.gz',
    params:
        soft=str(config['CHANGE']['FOLDERS']['SOFTWARE']),
        prog=str(config['CHANGE']['PROGRAMS']['FILTER']),
        xmx=str(config['CHANGE']['ISWEEP']['XMXMEM']),
        momcut=str(config['FIXED']['ISWEEP']['MOMCUTOFF']),
    resources:
        mem_gb=10
    shell:
        """
        zcat {input.ibd} | \
            java -Xmx{params.xmx}g -jar {params.soft}/{params.prog} "I" -8 0.00 {params.momcut} | \
            gzip > {output.fipass}
        """

rule count_ibdends_mom: # computing counts over windows
    input:
        filein='{cohort}/ibdsegs/ibdends/modified/mom/chr{num}.ibd.gz',
        mapin='{cohort}/maps/chr{num}.map',
    output:
        fileout='{cohort}/ibdsegs/ibdends/modified/mom/chr{num}.ibd.windowed.tsv.gz',
    params:
        scripts=str(config['CHANGE']['FOLDERS']['TERMINALSCRIPTS']),
        by=str(config['FIXED']['ISWEEP']['BY']),
        ge=str(config['FIXED']['ISWEEP']['GENOMEEND']),
        tc=str(config['FIXED']['ISWEEP']['TELOCUT']),
    resources:
        mem_gb=10
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

rule scan: # conduct a manhattan scan
    input:
        [macro+'/ibdsegs/ibdends/modified/scan/chr'+str(i)+'.ibd.windowed.tsv.gz' for i in range(low,high+1)],
        [macro+'/ibdsegs/ibdends/modified/mom/chr'+str(i)+'.ibd.windowed.tsv.gz' for i in range(low,high+1)],
    output:
        macro+'/excess.ibd.tsv',
    params:
        scripts=str(config['CHANGE']['FOLDERS']['TERMINALSCRIPTS']),
        folder=macro,
        chrlow=str(low),
        chrhigh=str(high),
        scansigma=str(config['FIXED']['ISWEEP']['SCANSIGMA']),
        telosigma=str(config['FIXED']['ISWEEP']['TELOSIGMA']),
    resources:
        mem_gb=10
    shell:
        """
        python {params.scripts}/scan-isweep.py \
            {params.folder} \
            {params.chrlow} \
            {params.chrhigh} \
            {params.scansigma} \
            {params.telosigma}
        """

rule ibdne:
    input:
        ibds=[macro+'/ibdsegs/ibdends/modified/scan/chr'+str(i)+'.ibd.gz' for i in range(low,high+1)],
        maps=[macro+'/maps/chr'+str(i)+'.map' for i in range(low,high+1)],
    output:
        ne=macro+'/ibdne.ne',
    params:
        chrlow=str(low),
        chrhigh=str(high),
        prefix=macro+'/ibdsegs/ibdends/modified/scan/',
        folder=macro,
        soft=str(config['CHANGE']['FOLDERS']['SOFTWARE']),
        prog=str(config['CHANGE']['PROGRAMS']['IBDNE']),
        xmx=str(config['CHANGE']['ISWEEP']['XMXMEM']),
        j='j',
    resources:
        mem_gb=100
    shell:
        """
        for {params.j} in $(seq {params.chrlow} 1 {params.chrhigh}); \
        do cat {params.folder}/maps/chr${j}.map \
        >> {params.folder}/maps/chr{params.chrlow}-{params.chrhigh}.map ; \
        done;
        for {params.j} in $(seq {params.chrlow} 1 {params.chrhigh}); \
        do zcat {params.prefix}chr${j}.ibd.gz \
        >> {params.prefix}chrall.ibd ; \
        done;
        cat {params.prefix}chrall.ibd | \
        java -Xmx{params.xmx}g -jar {params.soft}/{params.prog} \
        map={params.folder}/maps/chr{params.chrlow}-{params.chrhigh}.map \
        out={params.folder}/ibdne \
        rm {params.prefix}chrall.ibd
        """
