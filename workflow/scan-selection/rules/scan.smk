# conduct ibd calls and scan for isweep

localrules: plot_histogram

import os

# inputs, string management, count sample size, make mac
subsamplefile=str(config['CHANGE']['ISWEEP']['SUBSAMPLE'])
macro=str(config['CHANGE']['FOLDERS']['STUDY'])
samplesize=0
with open(subsamplefile,'r') as f:
    for line in f:
        samplesize+=1
# samplesize=str(samplesize)
ploidy=2
# ploidy=int(float(config['FIXED']['CANDHAPIBD']['PLOIDY']))
maf1=float(config['FIXED']['CANDHAPIBD']['MINMAF'])
mac1=int(ploidy*samplesize*maf1)
macro=str(config['CHANGE']['FOLDERS']['STUDY'])
low=int(float(str(config['CHANGE']['ISWEEP']['CHRLOW'])))
high=int(float(str(config['CHANGE']['ISWEEP']['CHRHIGH'])))
vcffolder=str(config['CHANGE']['EXISTING']['VCFS'])
vcfpre=str(config['CHANGE']['EXISTING']['VCFPRE'])
vcfsuf=str(config['CHANGE']['EXISTING']['VCFSUF'])

chromsstr = []
with open(macro+'/chromosome-sizes-kept.tsv','r') as f:
    f.readline()
    for line in f:
        chromsstr.append(line.strip().split('\t')[0])
chroms = [int(i) for i in chromsstr]

mgb=int(float(str(config['CHANGE']['ISWEEP']['XMXMEM'])))

rule hapibd: # candidate segments from hap-ibd.jar
    input:
        vcf=vcffolder + '/' + vcfpre + '{num}' + vcfsuf,
        mapfile='{cohort}/maps/chr{num}.map',
        excludesamples='{cohort}/excludesamples.txt',
    output:
        ibd='{cohort}/ibdsegs/hapibd/chr{num}.ibd.gz',
    params:
        minmac=str(mac1),
        out='{cohort}/ibdsegs/hapibd/chr{num}',
        xmx=str(config['CHANGE']['ISWEEP']['XMXMEM']),
        minsee=str(config['FIXED']['CANDHAPIBD']['MINSEED']),
        minext=str(config['FIXED']['CANDHAPIBD']['MINEXT']),
        minout=str(config['FIXED']['CANDHAPIBD']['MINOUT']),
    resources:
        mem_gb=mgb+1,
    shell:
        """
        java -Xmx{params.xmx}g -jar ../../software/hap-ibd.jar \
            gt={input.vcf} \
            map={input.mapfile} \
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
        mapfile='{cohort}/maps/chr{num}.map',
        ibd='{cohort}/ibdsegs/hapibd/chr{num}.ibd.gz',
        subsample='{cohort}/excludesamples.txt',
    output:
        ibd='{cohort}/ibdsegs/ibdends/chr{num}.ibd.gz',
    params:
        out='{cohort}/ibdsegs/ibdends/chr{num}',
        xmx=str(config['CHANGE']['ISWEEP']['XMXMEM']),
        maf=str(config['FIXED']['IBDENDS']['MINMAF']),
        qua=str(config['FIXED']['IBDENDS']['QUANTILES']),
        sam=str(config['FIXED']['IBDENDS']['NSAMPLES']),
        err=str(config['CHANGE']['IBDENDS']['ERRRATE']),
    resources:
        mem_gb=mgb+1,
    shell:
        """
        java -Xmx{params.xmx}g -jar ../../software/ibd-ends.jar \
            gt={input.vcf} \
            ibd={input.ibd} \
            map={input.mapfile} \
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
        cutibd='{cohort}/ibdsegs/ibdends/chr{num}.formatted.ibd.gz',
    shell:
        '''
        zcat {input.bigibd} | \
            tail -n +2 | \
            cut -f 1-5,10,11,12 | \
            gzip > {output.cutibd}
        '''

rule filter_ibdends_scan: # applying cutoffs
    input:
        ibd='{cohort}/ibdsegs/ibdends/chr{num}.formatted.ibd.gz',
    output:
        fipass='{cohort}/ibdsegs/ibdends/scan/chr{num}.ibd.gz',
    params:
        scancut=str(config['FIXED']['ISWEEP']['SCANCUTOFF']),
    shell:
        """
        python ../../scripts/utilities/filter-lines.py \
            --input_file {input.ibd} \
            --output_file {output.fipass} \
            --upper_bound {params.scancut}
        """

rule count_ibdends_scan: # computing counts over windows
    input:
        filein='{cohort}/ibdsegs/ibdends/scan/chr{num}.ibd.gz',
        mapin='{cohort}/maps/chr{num}.trimmed.map',
    output:
        fileout='{cohort}/ibdsegs/ibdends/scan/chr{num}.ibd.windowed.tsv.gz',
    shell:
        """
        python ../../scripts/utilities/count-ibd.py \
            --input_ibd_file {input.filein} \
            --input_map_file {input.mapin} \
            --output_file {output.fileout} \
            --start 5 \
            --end 6 \
        """

rule filter_ibdends_mle: # applying cutoffs
    input:
        ibd='{cohort}/ibdsegs/ibdends/chr{num}.formatted.ibd.gz',
    output:
        fipass='{cohort}/ibdsegs/ibdends/mle/chr{num}.ibd.gz',
    params:
        mlecut=str(config['FIXED']['ISWEEP']['MLECUTOFF']),
    shell:
        """
        python ../../scripts/utilities/filter-lines.py \
            --input_file {input.ibd} \
            --output_file {output.fipass} \
            --upper_bound {params.mlecut}
        """

rule count_ibdends_mle: # computing counts over windows
    input:
        filein='{cohort}/ibdsegs/ibdends/mle/chr{num}.ibd.gz',
        mapin='{cohort}/maps/chr{num}.trimmed.map',
    output:
        fileout='{cohort}/ibdsegs/ibdends/mle/chr{num}.ibd.windowed.tsv.gz',
    shell:
        """
        python ../../scripts/utilities/count-ibd.py \
            --input_ibd_file {input.filein} \
            --input_map_file {input.mapin} \
            --output_file {output.fileout} \
            --start 5 \
            --end 6 \
        """

rule scan: # conduct a manhattan scan
    input:
        [macro+'/ibdsegs/ibdends/scan/chr'+str(i)+'.ibd.windowed.tsv.gz' for i in chroms],
        [macro+'/ibdsegs/ibdends/mle/chr'+str(i)+'.ibd.windowed.tsv.gz' for i in chroms],
    output:
        scandata=macro+'/scan.ibd.tsv',
    params:
        study=macro,
        folder='/ibdsegs/ibdends/scan',
        chrlow=str(low),
        chrhigh=str(high),
        telosigma=str(config['FIXED']['ISWEEP']['TELOSIGMA']),
    shell:
        """
        python ../../scripts/scan/scan.py \
            --input_study {params.study} \
            --input_folder {params.folder} \
            --output_all_genome_file scan.ibd.tsv \
            --chr_low {params.chrlow} \
            --chr_high {params.chrhigh} \
            --input_prefix chr \
            --input_suffix .ibd.windowed.tsv.gz \
            --outlier_cutoff {params.telosigma} \
        """

rule plot_histogram:
    input:
        scandata=macro+'/scan.ibd.tsv'
    output:
        histogram=macro+'/zhistogram.png'
    shell:
        """
        python ../../scripts/plotting/plot-histogram.py \
            --input_file {input.scandata} \
            --output_file {output.histogram} \
            --chr_high 100 \
            --statistic Z \
            --xlabel z-score \
            --xupp 6. \
        """

