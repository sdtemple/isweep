# setup map files

macro=str(config['CHANGE']['FOLDERS']['STUDY'])
low=int(float(str(config['CHANGE']['ISWEEP']['CHRLOW'])))
high=int(float(str(config['CHANGE']['ISWEEP']['CHRHIGH'])))

# check this
rule write_map:
    input:
        oldmap='{config[CHANGE][EXISTING][MAPS]}/{config[CHANGE][EXISTING][MAPPRE]}{num}{config[CHANGE][EXISTING][MAPSUF]},
    output:
        newmap='{config[CHANGE][EXISTING][MAPS]}/maps/chr{num}.map',
    resources:
        mem_gb=5
    shell:
        'cp {input.oldmap} {output.newmap}'
