Usage
=====

You should always enter the isweep environment before running workflows: ``mamba activate isweep``.

You should always run the workflows on cluster nodes with ``nohup snakemake [...] --cluster "sbatch [...]" &`` to manage submissions in the background.

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

The ``worfklow/scan-selection`` implements the IBD rate selection scan with two multiple-testing corrections. You should use the `Snakefile-scan.smk` as input to the ``-s`` option.

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
* ``autocovariance.png```: model fit for autocovariances of IBD rates
* ```zhistogram.png``: an empirical distribution of IBD rates (should look Gaussian)
* ``roi.tsv``: summary table of the genome-wide significant loci
* ``fwer.analytical.txt``: details about parameters and estimates for multiple testing
* ``chromosome-sizes-kept.tsv``: numbers and cM sizes of chromosomes analyzed
* ``ibdsegs/``: folder and subfolders with detected IBD segments

The multiple-testing corrections are valid asymptotically (Temple and Thompson, 2024+).

There is a multiprocessing version using ``Snakefile-scan-mp.smk``, which may only be useful in enormous human biobanks.

.. _hard-sweeps:

Modeling hard sweeps
------------

Words

.. _case-control-scan:

Case-control scan
------------

Words

.. _prepare:

Preliminary material
------------

Words

Ploidy
------------

Words

