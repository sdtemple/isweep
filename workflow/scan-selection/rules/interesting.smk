# process the isweep scan
# form table of roi

macro=str(config['CHANGE']['FOLDERS']['STUDY'])
low=int(float(str(config['CHANGE']['ISWEEP']['CHRLOW'])))
high=int(float(str(config['CHANGE']['ISWEEP']['CHRHIGH'])))

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
            --file_input {input.filein} \
            --file_output {output.fileout} \
            --min_cm_region {params.cMgap}
        '''

rule make_roi_table:
    input:
        filein=macro+'/excess.region.ibd.tsv',
    output:
        fileout=macro+'/roi.tsv',
    params:
        folder=macro+'/ibdsegs/ibdends/scan'
        mbbuf=str(config['FIXED']['ISWEEP']['MBBUF']),
        cmcover=str(config['FIXED']['ISWEEP']['CMCOVER']),
        cmsmall=str(config['FIXED']['ISWEEP']['CMSMALL']),
    shell:
        """
        python ../../scripts/make-roi-table.py \
            --file_input {input.filein} \
            --file_output {output.fileout} \
            --file_prefix {params.folder}/chr \
            --file_suffix .ibd.windowed.tsv.gz \
            --cmcover {params.cmcover} \
            --cmsmall {params.cmsmall} \
            --mbbuffer {params.mbbuf}
        """