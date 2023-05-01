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

# check this
rule concat_map:
    input:
        [macro+'/maps/chr'+str(i)+'.map' for i in range(low,high+1)],
    output:
        allmap='{config[CHANGE][FOLDERS][STUDY]}/maps/chr{config[CHANGE][ISWEEP][CHRLOW]}-{config[CHANGE][ISWEEP][CHRHIGH]}.map',
    resources:
        mem_gb=5
    shell:
        'for j in $(seq {config[CHANGE][ISWEEP][CHRLOW]} 1 {config[CHANGE][ISWEEP][CHRHIGH]}); do cat {config[CHANGE][FOLDERS][STUDY]}/maps/chr${j}.map >> {output.allmap} ; done;'
