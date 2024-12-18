# implement family-wise error rate adjustments
# to determine genome-wide significance levels
# for the case-control mapping version

macro=str(config['CHANGE']['FOLDERS']['STUDY'])

pval = float(str(config['CHANGE']['ISWEEP']['CONFLEVEL']))
# you should choose one and stick with it. no p hacking.
stepsize = float(str(config['CHANGE']['ISWEEP']['CMSTEPSIZE']))
stepsize /= 100 # in morgans
genomesize = float(str(config['CHANGE']['ISWEEP']['GENOMESIZE']))
genomesize /= 100 # in morgans
chrlow = int(str(config['CHANGE']['ISWEEP']['CHRLOW']))
chrhigh = int(str(config['CHANGE']['ISWEEP']['CHRHIGH']))
numchr = chrhigh - chrlow + 1
chrsize = genomesize / numchr
covlen = float(str(config['FIXED']['ISWEEP']['AUTOCOVLEN']))
covlen= int(covlen / 100 / stepsize)
numsims = int(str(config['CHANGE']['ISWEEP']['SIMS']))

rule count_ibdends_case: # computing counts over windows for case and controls
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

# implement the siegmund and yakir (2007) method
# with log linear modeled autocovariance
# merged analytical_method, scan, and significance rules together
rule analytical_method:
    input:
        [macro+'/ibdsegs/ibdends/scan/chr'+str(i)+'.case.ibd.windowed.tsv.gz' for i in range(chrlow,chrhigh+1)],
    output:
        testing=macro+'/fwer.analytical.case.tsv',
        autocov0=macro+'/fwer.autocovariance.control.tsv',
        autocov1=macro+'/fwer.autocovariance.case.tsv',
        autocov=macro+'/fwer.autocovariance.diff.tsv',
        crosscov=macro+'/fwer.crosscovariance.diff.tsv',
        scan=macro+'/scan.case.ibd.tsv',
        excess=macro+'/excess.case.ibd.tsv',
    params:
       chrlow=str(chrlow),
       chrhigh=str(chrhigh),
       chrsize=str(chrsize),
       covlen=str(covlen),
       pval=str(pval),
       stepsize=str(stepsize),
       initcut=str(config['FIXED']['ISWEEP']['TELOSIGMA']),
       pre=macro+'/ibdsegs/ibdends/scan/chr',
       autocov=macro+'/fwer.autocovariance'
    shell:
        """
        python ../../scripts/scan/multiple-testing-analytical-case.py \
            --output_testing_file {output.testing} \
            --output_autocov_prefix {params.autocov} \
            --output_crosscov_file {output.crosscov} \
            --output_scan_file {output.scan} \
            --output_excess_file {output.excess} \
            --input_prefix {params.pre}\
            --input_suffix .case.ibd.windowed.tsv.gz \
            --chr_low {params.chrlow} \
            --chr_high {params.chrhigh} \
            --chr_average_size {params.chrsize} \
            --cM_step_size {params.stepsize} \
            --autocovariance_steps {params.covlen} \
            --confidence_level {params.pval} \
            --outlier_cutoff {params.initcut} \
            --counts_column0 COUNT0 \
            --counts_column1 COUNT1 \
        """

# # maybe implement later
# rule simulation_method:
#     input:
#         pass
#     output:
#         pass
#     params:
#         pass
#     shell:
#         """
#         echo 'Hello world'
#         """

rule plot_autocovariance:
    input:
        testing=macro+'/fwer.analytical.case.tsv',
        autocov=macro+'/fwer.autocovariance.diff.tsv',
    output:
        figure=macro+'/autocovariance.diff.png',
    params:
        prefix=macro+'/autocovariance.diff'
    shell:
        """
        python ../../scripts/plotting/plot-autocovariance.py \
            --input_autocov_file {input.autocov} \
            --input_analytical_file {input.testing} \
            --output_prefix {params.prefix} \
            --theta_type 'estimated-theta:' \
            --title 'Difference in IBD rates' \
        """

rule plot_autocovariance_control:
    input:
        testing=macro+'/fwer.analytical.case.tsv',
        autocov=macro+'/fwer.autocovariance.control.tsv',
    output:
        figure=macro+'/autocovariance.control.png',
    params:
        prefix=macro+'/autocovariance.control'
    shell:
        """
        python ../../scripts/plotting/plot-autocovariance.py \
            --input_autocov_file {input.autocov} \
            --input_analytical_file {input.testing} \
            --output_prefix {params.prefix} \
            --theta_type 'estimated-theta0:' \
            --title 'Control IBD rates' \
        """

rule plot_autocovariance_case:
    input:
        testing=macro+'/fwer.analytical.case.tsv',
        autocov=macro+'/fwer.autocovariance.case.tsv',
    output:
        figure=macro+'/autocovariance.case.png',
    params:
        prefix=macro+'/autocovariance.case'
    shell:
        """
        python ../../scripts/plotting/plot-autocovariance.py \
            --input_autocov_file {input.autocov} \
            --input_analytical_file {input.testing} \
            --output_prefix {params.prefix} \
            --theta_type 'estimated-theta1:' \
            --title 'Case IBD rates' \
        """

rule plot_histogram_diff:
    input:
        scandata=macro+'/scan.case.ibd.tsv'
    output:
        histogram=macro+'/zhistogram.diff.png'
    shell:
        """
        python ../../scripts/plotting/plot-histogram.py \
            --input_file {input.scandata} \
            --output_file {output.histogram} \
            --chr_high 100 \
            --statistic ZDIFFZ \
            --xlabel z-score \
            --xupp 6. \
            --title 'Difference of IBD rates' \
        """

rule plot_histogram_case:
    input:
        scandata=macro+'/scan.case.ibd.tsv'
    output:
        histogram=macro+'/zhistogram.case.png'
    shell:
        """
        python ../../scripts/plotting/plot-histogram.py \
            --input_file {input.scandata} \
            --output_file {output.histogram} \
            --chr_high 100 \
            --statistic Z1 \
            --xlabel z-score \
            --xupp 6. \
            --title 'Case IBD rates' \
        """

rule plot_histogram_control:
    input:
        scandata=macro+'/scan.case.ibd.tsv'
    output:
        histogram=macro+'/zhistogram.control.png'
    shell:
        """
        python ../../scripts/plotting/plot-histogram.py \
            --input_file {input.scandata} \
            --output_file {output.histogram} \
            --chr_high 100 \
            --statistic Z0 \
            --xlabel z-score \
            --xupp 6. \
            --title 'Control IBD rates' \
        """