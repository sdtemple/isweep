# process the isweep scan
# form table of roi

macro=str(config['CHANGE']['FOLDERS']['STUDY'])
low=int(float(str(config['CHANGE']['ISWEEP']['CHRLOW'])))
high=int(float(str(config['CHANGE']['ISWEEP']['CHRHIGH'])))

# rule plot_manhattan:
#     input:
#     output:
#     params:
#     shell:
#         """
#         python ../../scripts/plot-manhattan-pipeline.py \
#         """

rule excess_region: # concatenate regions of excess IBD
    input:
        filein=macro+'/excess.ibd.tsv',
    output:
        fileout=macro+'/excess.region.ibd.tsv',
    params:
        cMgap=str(config['FIXED']['ISWEEP']['CMGAP']),
    shell:
        '''
        python ../../scripts/excess-region.py \
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
        python ../../scripts/make-roi-table.py \
            --input_file {input.filein} \
            --output_file {output.fileout} \
            --input_prefix {params.folder}/chr \
            --input_suffix .ibd.windowed.tsv.gz \
            --cM_cover {params.cmcover} \
            --cM_small {params.cmsmall} \
            --Mb_buffer {params.mbbuf}
        """
