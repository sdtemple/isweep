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
        cMgap='{config[FIXED][ISWEEP][CMGAP]}',
    shell:
        'python {config[CHANGE][FOLDERS][TERMINALSCRIPTS]}/excess-region.py {input.filein} {output.fileout} {config[FIXED][ISWEEP][CMGAP]}'

rule make_roi_table:
    input:
        filein=macro+'/excess.region.ibd.tsv',
    output:
        fileout=macro+'/roi.tsv',
    shell:
        'mkdir -p {config[CHANGE][FOLDERS][STUDY]}/plots; mkdir -p {config[CHANGE][FOLDERS][STUDY]}/stats; python {config[CHANGE][FOLDERS][TERMINALSCRIPTS]}/make-roi-table.py {config[CHANGE][FOLDERS][STUDY]}'
