# Implement the multiple-testing corrections
# to determine genome-wide significance levels
# that control the family-wise error rate.
# This is the case-control scan version.

import pandas as pd

localrules: \
    plot_autocovariance, \
    plot_autocovariance_case, \
    plot_histogram_control, \
    plot_histogram_case, \
    plot_histogram_diff, \
    plot_autocovariance_control

macro=str(config['change']['files']['study'])

pval = float(str(config['change']['isweep']['confidence_level']))
# you should choose one and stick with it. no p hacking.
stepsize = float(str(config['change']['isweep']['step_size_cm']))
stepsize /= 100 # in morgans
chromosome_sizes = pd.read_csv(macro+'/chromosome-sizes-kept.tsv',sep='\t')
genomesize = chromosome_sizes['CMSIZE'].sum()
chroms2 = chromosome_sizes.CHROM.astype(int).tolist()
telocutting = float(str(config['change']['isweep']['scan_cutoff']))
numchr = chromosome_sizes.shape[0]
genomesize -= numchr * telocutting * 2
genomesize /= 100 # in morgans
chrsize = genomesize / numchr
covlen = float(str(config['fixed']['isweep']['auto_covariance_length']))
covlen= int(covlen / 100 / stepsize)
numsims = int(str(config['change']['isweep']['num_sims']))

rule count_ibdends_case: # computing counts over windows for case and controls
    input:
        filein='{cohort}/ibdsegs/ibdends/scan/chr{num}.ibd.gz',
        mapin='{cohort}/maps/chr{num}.trimmed.map',
    output:
        fileout='{cohort}/ibdsegs/ibdends/scan/chr{num}.case.ibd.windowed.tsv.gz',
    params:
        cases=str(config['change']['files']['cases']),
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
        [macro+'/ibdsegs/ibdends/scan/chr'+str(i)+'.case.ibd.windowed.tsv.gz' for i in chroms2],
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
       initcut=str(config['fixed']['isweep']['outlier_cutoff']),
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

rule simulation_method:
    input:
        testing=macro+'/fwer.analytical.case.tsv',
    output:
        testing=macro+'/fwer.simulation.case.txt',
    params:
        numsims=str(config['change']['isweep']['num_sims']),
    shell:
        """
        python ../../scripts/scan/multiple-testing-simulation-case-pipeline.py \
            --input_file {input.testing} \
            --output_file {output.testing} \
            --num_sims {params.numsims} \
        """

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