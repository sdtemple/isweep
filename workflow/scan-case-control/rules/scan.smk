

rule count_ibdends_case: # computing counts over windows
    input:
        filein='{cohort}/ibdsegs/ibdends/scan/chr{num}.ibd.gz',
        mapin='{cohort}/maps/chr{num}.trimmed.map',
    output:
        fileout='{cohort}/ibdsegs/ibdends/scan/chr{num}.case.ibd.windowed.tsv.gz',
    params:
        cases=str(config['CHANGE']['ISWEEP']['CASES']),
    shell:
        """
        python ../../scripts/utilities/count-ibd-case.py \
            --input_ibd_file {input.filein} \
            --input_map_file {input.mapin} \
            --input_case_file {params.cases} \
            --output_file {output.fileout} \
            --start 5 \
            --end 6 \
            --ind1 0 \
            --ind2 2 \
        """

rule scan:
    input:
        pass
    output:
        pass
    params:
        pass
    shell:
        """
        echo 'Hello world'
        """

rule plot_histogram:
    input:
        pass
    output:
        pass
    params:
        pass
    shell:
        """
        echo 'Hello world'
        """

rule plot_scan:
    input:
        pass
    output:
        pass
    params:
        pass
    shell:
        """
        echo 'Hello world'
        """