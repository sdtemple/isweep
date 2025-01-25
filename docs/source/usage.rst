Usage
=====

You should always enter the isweep environment before running workflows: ``mamba activate isweep``.

You should always run the workflows on cluster nodes with ``nohup snakemake [...] --cluster "sbatch [...]" &`` to manage submissions in the background.

The penultimate rule of all workflows is to copy the YAML parameters file to your analysis folder, supporting reproducibility and logging.

I regularly use these ``snakemake`` options:

* ``-s``: point to the right Snakefile
* ``--configfile``: point to your parameters file
* ``--jobs``: how many jobs can queue at once
* ``--cluster "[sbatch ...]"``
* ``-n``: dry run to see what the workflow will run
* ``--latency-wait 200``
* ``--keep-going``
* ``-c1``

I regularly use these SLURM options:

* ``--cpus-per-task``: you should max out the CPUs on a node
* ``--mem``: you should almost max out the memory on a node
* ``-e {rule}`` and ``-o {rule}``: many jobs will be run, so you should define file locations for stdout and stderr
* ``--job-name={rule}``
* ``--mail-type=END`` and ``mail-user``: careful to not send yourself too many emails
* ``--partition``

.. _selection-scan:

Selection scan
------------

The ``worfklow/scan-selection`` implements the IBD rate selection scan with two multiple-testing corrections. You should use the `Snakefile-scan.smk` file as input to the ``-s`` option.

Recipe YAML files to modify are ``sequence.yaml`` and ``array.yaml`` for WGS and SNP array data, respectively. There is a hierarchy of ``change`` versus ``fixed`` parameters, where ``change`` you should modify for your dataset and ``fixed`` you should reach out for advice.

The parameters are:

* Many parameters under ``files`` determine where your data is and where you want outputs to be.
* You can use ``chromosome_low`` and ``chromosome_high`` to determine a range of such to study. All chromosome ``.vcf.gz`` and ``.map`` must be numbered.
* ``subsample``: text file with sample IDs in VCF files, which can/likely is a subset of larger consortium dataset
* ``ibd_ends:error_rate``: set this from estimated error in pilot study of your smallest chromosomes (log files from ibd-ends software)
* ``ploidy``: if your ploidy is not 1 or 2, see :ref:`ploidy`
* ``step_size_cm``: you perform a hypothesis every X.XX centiMorgans
* ``scan_cutoff``: minimum length of detected IBD segments (recommended >= 2.0 or >= 3.0)
* ``confidence_level``: the family-wise error rate you want to control (e.g., 0.05)
* ``num_sims``: simulations to derive one of the multiple-testing corrections
* ``chromosome_exclude``: empty, or a text file with chromosome numbers to not analyze

Chromosomes smaller that the centiMorgan parameter ``fixed:isweep:chromosome_size_cutoff`` will not be analyzed. You cannot reliably estimate the autocovariance decay of IBD rates in mini chromosomes.

The outputs are:

* ``scan.png``: the IBD rates along autosomes and the estimated significance thresholds
* ``scan.modified.ibd.tsv``: the difference in IBD rates along the autosomes in tabular format
* ``autocovariance.png```: model fit for autocovariances of IBD rates
* ``zhistogram.png``: an empirical distribution of IBD rates (should look Gaussian)
* ``roi.tsv``: summary table of the genome-wide significant loci
* ``fwer.analytical.txt``: details about parameters and estimates for multiple testing
* ``chromosome-sizes-kept.tsv``: numbers and cM sizes of chromosomes analyzed
* ``ibdsegs/``: folder and subfolders with detected IBD segments

The multiple-testing corrections are valid asymptotically (Temple and Thompson, 2024+). You can look at the IBD rate histogram to visually assess such. Be wary of IBD rates being zero truncated in small samples.

There is a multiprocessing version using ``Snakefile-scan-mp.smk``, which may only be useful in enormous human biobanks.

.. _hard-sweeps:

Modeling hard sweeps
------------

The ``worfklow/model-selection`` estimates frequencies, locations, and selection coefficients of loci detected in the :ref:selection-scan. This workflow must be run after the selection scan. You should use the `Snakefile-roi.smk` file as input to the ``-s`` option.

The recipe YAML file to modify is ``sweep.yaml``. There is a hierarchy of ``change`` versus ``fixed`` parameters, where ``change`` you should modify for your dataset and ``fixed`` you should reach out for advice.

The parameters are:

* Many parameters under ``files`` determine where your data is and where you want outputs to be.
* ``regions_of_interest``: these are the loci to analyse. The default are those GW significant in the scan. You can delete some, or rename the GW significant "hits".
* ``chromosome_prefix``: this is the name ``chr`` or blank that you see when you run ``bcftools query -f "%CHROM\n" chr.vcf.gz | head``.
* ``ploidy``: if your ploidy is not 1 or 2, see :ref:`ploidy`
* ``Ne``: an estimate of recent effective population sizes (IBDNe text file format)

You can change the genic selection model in ``roi.tsv`` to "a" for additive, "m" for multiplicative, "d" for dominance, and "r" for recessive. You can also change alpha, which determines the (1-alpha) percent confidence intervals.

The script ``scripts/run-ibdne.sh`` runs IBDNe, which is good for populations with exponential growth. You may want to consider another Ne estimator as well.

.. code-block:: shell

   sbatch [...] run-ibdne.sh [ibdne-jar] [memory-in-Gb] [main-folder-of-study] [path-to-subfolder-with-ibd-data] [chromosome_low] [chromosome_high] [output_file] [random_seed]

The outputs are:

* ``summary.hap.norm.tsv``: sweep model estimates for best haplotype-based analysis and Gaussian confidence intervals for selection coefficient
* ``summary.snp.norm.tsv``: sweep model estimates for best SNP-based analysis and Gaussian confidence intervals for selection coefficient
* ``hit*/second.ranks.tsv.gz``: alleles with putative evidence for selection (or strong correlation with a selected allele)
* ``hit*/outlier*.txt``: files with sample haplotype IDs in clusters on excess IBD sharing

The Gaussian bootstrap intervals are valid asymptotically (Temple and Thompson, 2024+). You can uncomment lines in `rule all` of the `Snakefile-roi.smk` to get percentile-based bootstrap intervals.

There is a multiprocessing version using ``Snakefile-scan-mp.smk``, which may only be useful in enormous human biobanks.

.. _case-control-scan:

Case-control scan
------------

The ``worfklow/scan-case-control`` implements the difference in IBD rates scan with two multiple-testing corrections. You should use the `Snakefile-case.smk` file as input to the ``-s`` option.

You must run this workflow after the selection scan workflow (where the IBD segments are detected). You should scrutinize the results to see if strong selection confounds your case-control study.

The recipe YAML file to modify is ``case.yaml``. The parameters are nearly all the same as in :ref:_selection-scan. The ``case`` parameter is a two-column text file with sample IDs and binary phenotypes.

The outputs have the same nomenclature as in the selection scan workflow, but ``.case.`` and ``.control.`` is inserted in file names:

* ``scan.case.control.png``: the standardized difference in IBD rates along autosomes and the estimated significance thresholds
* ``scan.case.ibd.tsv``: the difference in IBD rates along the autosomes in tabular format 
* ``roi.case.tsv``: summary table of the genome-wide significant loci
* ``fwer.analytical.case.txt``: details about parameters and estimates for multiple testing

The multiple-testing corrections are valid asymptotically (Temple and Thompson, 2024+). You can look at the IBD rate histogram to visually assess such. Be wary of IBD rates being zero truncated in small samples. 

There is a multiprocessing version using ``Snakefile-case-mp.smk``, which may only be useful in enormous human biobanks.

.. _prepare:

Preliminary material
------------

Words

Ploidy
------------

Words

