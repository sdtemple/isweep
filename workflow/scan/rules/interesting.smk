# process the isweep scan
# form table of roi
# seth temple, sdtemple@uw.edu
# may 3, 2023

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
        scripts=str(config['CHANGE']['FOLDERS']['TERMINALSCRIPTS']),
    shell:
        'python {params.scripts}/excess-region.py {input.filein} {output.fileout} {params.cMgap}'

rule make_roi_table:
    input:
        maps=[macro+'/maps/chr'+str(i)+'.map' for i in range(low,high+1)],
        filein=macro+'/excess.region.ibd.tsv',
    output:
        fileout=macro+'/roi.tsv',
    params:
        study=macro,
        mapfolder=macro+'/maps',
        scripts=str(config['CHANGE']['FOLDERS']['TERMINALSCRIPTS']),
        mbbuf=str(config['FIXED']['ISWEEP']['MBBUF']),
    shell:
        """
        mkdir -p {params.study}/plots
        mkdir -p {params.study}/stats
        python {params.scripts}/make-roi-table.py {params.study} {params.mbbuf} {params.mapfolder}
        """
