# implement family-wise error rate adjustments
# to determine genome-wide significance levels

import pandas as pd

macro=str(config['CHANGE']['FOLDERS']['STUDY'])

pval = float(str(config['CHANGE']['ISWEEP']['CONFLEVEL']))
# you should choose one and stick with it. no p hacking.
stepsize = float(str(config['CHANGE']['ISWEEP']['CMSTEPSIZE']))
stepsize /= 100 # in morgans
chromosome_sizes = pd.read_csv(macro+'/chromosome-sizes-kept.tsv',sep='\t')
genomesize = chromosome_sizes['CMSIZE'].sum()
numchr = chromosome_sizes.shape[0]
telocutting = float(str(config['FIXED']['ISWEEP']['SCANCUTOFF']))
genomesize -= numchr * telocutting * 2
genomesize /= 100 # in morgans
chrlow = int(str(config['CHANGE']['ISWEEP']['CHRLOW']))
chrhigh = int(str(config['CHANGE']['ISWEEP']['CHRHIGH']))
chrsize = genomesize / numchr
covlen = float(str(config['FIXED']['ISWEEP']['AUTOCOVLEN']))
covlen= int(covlen / 100 / stepsize)
numsims = int(str(config['CHANGE']['ISWEEP']['SIMS'])) 

# implement the siegmund and yakir (2007) method
# with log linear modeled autocovariance
rule analytical_method:
    input:
        scan=macro+'/scan.ibd.tsv',
    output:
        testing=macro+'/fwer.analytical.tsv',
        autocov=macro+'/fwer.autocovariance.tsv',
    params:
       chrlow=str(chrlow),
       chrhigh=str(chrhigh),
       chrsize=str(chrsize),
       covlen=str(covlen),
       pval=str(pval),
       stepsize=str(stepsize),
       initcut=str(config['FIXED']['ISWEEP']['TELOSIGMA']),
       pre=macro+'/ibdsegs/ibdends/scan/chr',
    shell:
        """
        python ../../scripts/scan/multiple-testing-analytical.py \
            --output_testing_file {output.testing} \
            --output_autocov_file {output.autocov} \
            --input_prefix {params.pre}\
            --input_suffix .ibd.windowed.tsv.gz \
            --chr_low {params.chrlow} \
            --chr_high {params.chrhigh} \
            --chr_average_size {params.chrsize} \
            --cM_step_size {params.stepsize} \
            --autocovariance_steps {params.covlen} \
            --confidence_level {params.pval} \
            --outlier_cutoff {params.initcut} \
            --counts_column COUNT \
        """

# simulate ornstein-uhlenbeck processes to get significance level
rule simulation_method:
    input:
        analytical=macro+'/fwer.analytical.tsv',
    output:
        simulation=macro+'/fwer.simulation.txt',
    params:
       numsims=str(config['CHANGE']['ISWEEP']['SIMS']),
    shell:
        """
        python ../../scripts/scan/multiple-testing-simulation-pipeline.py \
            --input_file {input.analytical} \
            --output_file {output.simulation} \
            --num_sims {params.numsims} \
        """

rule significance:
    input:
        analytical=macro+'/fwer.analytical.tsv',
        simulation=macro+'/fwer.simulation.txt',
        scan=macro+'/scan.ibd.tsv',
    output:
        scan=macro+'/scan.modified.ibd.tsv',
        excess=macro+'/excess.ibd.tsv',
    shell:
        """
        python ../../scripts/scan/significance.py \
            --input_ibd_file {input.scan} \
            --input_analytical_file {input.analytical} \
            --input_simulation_file {input.simulation} \
            --output_modified_file {output.scan} \
            --output_excess_file {output.excess} \
        """

rule plot_autocovariance:
    input:
        testing=macro+'/fwer.analytical.tsv',
        autocov=macro+'/fwer.autocovariance.tsv',
    output:
        figure=macro+'/autocovariance.png',
    params:
        prefix=macro+'/autocovariance'
    shell:
        """
        python ../../scripts/plotting/plot-autocovariance.py \
            --input_autocov_file {input.autocov} \
            --input_analytical_file {input.testing} \
            --output_prefix {params.prefix} \
        """