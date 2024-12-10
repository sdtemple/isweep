# process the isweep scan
# form table of roi

macro=str(config['CHANGE']['FOLDERS']['STUDY'])
low=int(float(str(config['CHANGE']['ISWEEP']['CHRLOW'])))
high=int(float(str(config['CHANGE']['ISWEEP']['CHRHIGH'])))

subsamplefile=str(config['CHANGE']['ISWEEP']['SUBSAMPLE'])
macro=str(config['CHANGE']['FOLDERS']['STUDY'])
samplesize=0
with open(subsamplefile,'r') as f:
    for line in f:
        samplesize+=1

numsims = int(str(config['CHANGE']['ISWEEP']['SIMS']))

rule plot_scan:
    input:
        ibd=macro+'/scan.ibd.tsv',
    output:
        png=macro+'/scan.png',
    params:
        prefix=macro+'/scan',
        samplesize=str(samplesize),
        numsims=str(numsims),
        chrlow=str(low),
        chrhigh=str(high),
    shell:
        """
        python ../../scripts/plotting/plot-scan-pipeline.py \
            --input_file {input.ibd} \
            --output_prefix {params.prefix} \
            --sample_size {params.samplesize} \
            --ploidy 2 \
            --heuristic_cutoff 4 \
            --num_sims {params.numsims} \
            --chr_low {params.chrlow} \
            --chr_high {params.chrhigh} \
            --statistic COUNT \
            --ylabel 'IBD rate' \
            --xlabel 'Position along genome' \
        """

rule excess_region: # concatenate regions of excess IBD
    input:
        filein=macro+'/excess.ibd.tsv',
    output:
        fileout=macro+'/excess.region.ibd.tsv',
    params:
        cMgap=str(config['FIXED']['ISWEEP']['CMGAP']),
    shell:
        '''
        python ../../scripts/scan/excess-region.py \
            --input_file {input.filein} \
            --output_file {output.fileout} \
            --max_cM_gap {params.cMgap}
        '''

rule make_roi_table:
    input:
        filein=macro+'/excess.region.ibd.tsv',
    output:
        fileout=macro+'/roi.tsv',
    params:
        folder=macro+'/ibdsegs/ibdends/scan',
        mbbuf=str(config['FIXED']['ISWEEP']['MBBUF']),
        cmcover=str(config['FIXED']['ISWEEP']['CMCOVER']),
        cmsmall=str(config['FIXED']['ISWEEP']['CMSMALL']),
    shell:
        """
        python ../../scripts/scan/make-roi-table.py \
            --input_file {input.filein} \
            --output_file {output.fileout} \
            --input_prefix {params.folder}/chr \
            --input_suffix .ibd.windowed.tsv.gz \
            --cM_cover {params.cmcover} \
            --cM_small {params.cmsmall} \
            --Mb_buffer {params.mbbuf}
        """
