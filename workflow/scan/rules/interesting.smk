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
    resources:
        mem_gb=5
    script:
        '{config[CHANGE][FOLDERS][SNAKESCRIPTS]}/excess-region.py'

# check this
rule concat_windowed:
    input:
        [macro+'/ibdsegs/ibdends/modified/scan/chr'+str(i)+'.ibd.windowed.tsv.gz' for i in range(low,high+1)],
    output:
        allwindows=macro+'/ibd.windowed.tsv',
    resources:
        mem_gb=5
    shell:
        'for j in $(seq {config[CHANGE][ISWEEP][CHRLOW]} 1 {config[CHANGE][ISWEEP][CHRLOW]}); do zcat {config[CHANGE][FOLDERS][STUDY]}/ibdsegs/ibdends/modified/scan/chr${j}.ibd.windowed.tsz.gz | tail -n +2 >> {output.allwindows} ; done;'

# work this
# check this
rule make_roi_table:
    input:
        windowed=macro+'/ibd.windowed.tsv',
        excessregion=macro+'/excess.region.ibd.tsv',
    output:
        roitable=macro+'/roi.table.tsv',
    resources:
        mem_gb=5
    script:
        '{config[CHANGE][FOLDERS][SNAKESCRIPTS]}/make-roi-table.py'

# check this
rule make_roi_folders:
    input:
        roitable=macro+'/roi.table.tsv',
    output:
        macro+'/roifolders.done',
    script:
        '{config[CHANGE][FOLDERS][SNAKESCRIPTS]}/make-roi-folders.py'
