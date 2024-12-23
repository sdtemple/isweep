# Summarize the core findings of the selection scan.
# The output is a regions of interest tabular file.

localrules: plot_scan, make_roi_table, excess_region

subsamplefile=str(config['change']['files']['cases'])
samplesize=0
with open(subsamplefile,'r') as f:
    for line in f:
        samplesize+=1

numsims = int(str(config['change']['isweep']['num_sims']))

rule plot_scan:
    input:
        ibd=macro+'/scan.case.ibd.tsv',
    output:
        png=macro+'/scan.case.control.png',
    params:
        prefix=macro+'/scan.case.control',
        samplesize=str(samplesize),
        numsims=str(numsims),
        chrlow=str(low),
        chrhigh=str(high),
    shell:
        """
        python ../../scripts/plotting/plot-scan-case-pipeline.py \
            --input_file {input.ibd} \
            --output_prefix {params.prefix} \
            --chr_low {params.chrlow} \
            --chr_high {params.chrhigh} \
            --statistic ZDIFFZ \
            --ylabel 'Standardized rate difference' \
            --xlabel 'Position along genome' \
            --yupp 8 \
        """

rule excess_region: # concatenate regions of excess IBD
    input:
        filein=macro+'/excess.case.ibd.tsv',
    output:
        fileout=macro+'/excess.case.region.ibd.tsv',
    params:
        cMgap=str(config['fixed']['isweep']['cm_gap']),
    shell:
        '''
        python ../../scripts/scan/excess-region.py \
            --input_file {input.filein} \
            --output_file {output.fileout} \
            --max_cM_gap {params.cMgap} \
            --statistic ZDIFFZ \
        '''

# work on this
rule make_roi_table:
    input:
        filein=macro+'/excess.case.region.ibd.tsv',
        filein2=macro+'/scan.case.ibd.tsv',
    output:
        fileout=macro+'/roi.case.tsv',
    params:
        mbbuf=str(config['fixed']['isweep']['mega_base_buffer']),
        cmcover=str(config['fixed']['isweep']['covered_cm_region']),
        cmsmall=str(config['fixed']['isweep']['small_cm_region']),
    shell:
        """
        python ../../scripts/scan/make-roi-table.py \
            --input_excess_file {input.filein} \
            --input_scan_file {input.filein2} \
            --output_file {output.fileout} \
            --cM_cover {params.cmcover} \
            --cM_small {params.cmsmall} \
            --Mb_buffer {params.mbbuf} \
            --statistic ZDIFFZ \
            --sweep 0 \
        """
