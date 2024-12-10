# make genetic map with specified cM step size
rule interpolate_map:
    input:
        mapfile='{cohort}/maps/chr{num}.map',
    output:
        mapfile='{cohort}/maps/chr{num}.interpolated.map',
    params:
        stepsize=str(config['CHANGE']['ISWEEP']['CMSTEPSIZE']),
    shell:
        '''
        python ../../scripts/interpolate-map.py \
            --input_map_file {input.mapfile} \
            --output_map_file {output.mapfile} \
            --cM_step_size {params.stepsize} \
        '''

# cut the telomeres by length equal to scan detection threshold
rule trim_telomeres:
    input:
        mapfile='{cohort}/maps/chr{num}.interpolated.map',
    output:
        mapfile='{cohort}/maps/chr{num}.trimmed.map',
    params:
        trim=str(config['FIXED']['ISWEEP']['SCANCUTOFF']),
    shell:
        '''
        python ../../scripts/trim-telomere-map.py \
            --input_map_file {input.mapfile} \
            --output_map_file {output.mapfile} \
            --cM_telo_trim {params.trim} \
        '''