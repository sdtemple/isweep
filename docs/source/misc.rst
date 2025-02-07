Miscellaneous
=====

Installing a fast package manager
##############

I like to use mamba from miniforge as my package manager.

.. code-block:: shell

    wget https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-x86_64.sh
    bash Miniforge3-Linux-x86_64.sh
    mamba

If the mamba command does not work,

1. ``vim .bashrc`` 
2. Put in a line and save ``alias mamba="/path/to/miniforge3/bin/mamba``
3. ``source .bashrc``
4. Sign out and sign back in of terminal


.. _snakemake-options:

Snakemake options
##############

I regularly use these options.

* ``-s``: point to the right Snakefile
* ``--configfile``: point to your parameters file
* ``--jobs``: how many jobs can queue at once
* ``--cluster "[sbatch ...]"``
* ``-n``: dry run to see what the workflow will run
* ``--latency-wait 200``
* ``--keep-going``
* ``-c1``

.. _slurm-options:

Slurm options
##############

I regularly use these options.

* ``--cpus-per-task``: you should max out the CPUs on a node
* ``--mem``: you should almost max out the memory on a node
* ``-e {rule}`` and ``-o {rule}``: many jobs will be run, so you should define file locations for stdout and stderr
* ``--job-name={rule}``
* ``--mail-type=END`` and ``mail-user``: careful to not send yourself too many emails
* ``--partition``

.. _testing-workflows:

Testing workflows
##############

There is a `Zenodo repository <https://zenodo.org/records/14744019>`_ with some small simulated data examples. The ``README.md`` file gives commands for how I made the test data, and how to run the workflows. Note that the YAML and commands are relative to my local paths, so you will need to modify the paths to your data.

You can first locally run the hard sweep estimation and case-control scan with the data I left in the Zenodo repository.

You can run the phasing and LAI workflow locally, and it will finish in less than 10 minutes, but I recommend using a cluster.

You should run the selection scan with cluster resources.

.. note::

   You are encouraged to test out the workflows if you have >= 500 samples.

.. _ploidy-extension:

Ploidy extension
##############

VCF files with more than 1 or 2 ploidy are minimally supported. The cheat code is to treat them like haploid VCFs for the software using ``scripts/utilities/ploidy-conversion.py``. Let sample 1 have the genotype 0|0|0|1. The script will convert this into 4 samples with a haplotype index appended and the genotypes 0, 0, 0, 1.

You may also need to ``bgzip`` the output files again, and ``tabix -fp bcf``. The output works immediately in the selection scan, but ``tabix``-ing can fail for the hard sweep estimation.

For nondiploidy, you should set ploidy to be 1 in all configuration files. For modeling hard sweeps, you should make sure that your Ne file is scaled by the ploidy. For example, if your Ne file is w.r.t. the number of tetraploids, you should multiply the discrete Ne's by 4. Moreover, the sweep model will assume the formulas for haploid genic selection.

I am not an expert in nondiploidy. This cheat code may not be reasonable for your data.

Other considerations
##############

* There is limited statistical power in the selection scan with high cM length thresholds (>= 4.0).
* For humans, using pedigree-based maps like the deCODE map are crucial for accurate IBD segment detection. Non-pedigree based maps may be suitable in non-humans, as long as the estimated recombination rates are accurate enough for IBD segment detection.
* The p values assume the null model in the scans. If the histograms are far from Gaussian, you should not trust the p values.
* The null model is that there is a genome-wide mean IBD rate. If there are apparently two or more subsets of chromosomes with a different mean IBD rate, you should run such subsets separately using ``chromosome_exclude`` in the YAML file.
* Be cautious about interpretation of results near centromeres, where IBD segment detection is difficult.
* You could analyze recombining sex chromosomes solo, but estimates of the genome-wide significance level will be noisy. You should give the chromosome a pseudo number, e.g., human chromosome X as chr23.
* You can use ``scripts/plotting/plot-sweep.py`` to make figures like those in Temple, Waples, and Browning (2024). The file assumes you use Gaussian-based intervals (``scripts/model/estimate-norm.py``).
* Parameters for Mb buffer, window and haplotype sizes and steps in ``sequence.yaml`` and ``sweep.yaml`` are based on 1 cM ~~ 1 Mb. You may want to scale these accordingly if your species has a very different recombination rate.

.. note::

   The Temple and Thompson conditions, under which the scan is asymptotically valid, are:

   1. Sample size squared large relative to population size times cM length threshold (n^2 = o(Nw))
   2. Scaled population size large relative to sample size (Nw = o(n))

   The Gaussian model is often reasonable whenever sample size and scaled population size are large, even if the above conditions don't hold.

   There is a generalization of the main Temple and Thompson CLT for flexible demographic scenarios, i.e., large recent effective population sizes. 

Potential errors
##############

* SLURM jobs may fail at the Beagle or ibd-ends steps because of RAM. Re-run with more resources.
* Make sure your VCF files are tab-indexed (``tabix -fp vcf [...]``)
* You are not in the ``mamba activate isweep`` Python environment. (Failure b/c you don't have some package.)
* Binary incompatibility with ``statsmodels`` in ``multiple-testing-analytica-*.py``. Run ``pip install --upgrade numpy statsmodels scipy pandas networkx matplotlib seaborn``.
* A locus fails at the ``rule first_rank`` in ``workflow/model-selection`` because no excess IBD sharing group exists
* Sometimes the `Browning Lab software <https://github.com/browning-lab/>`_ (JAR files) on GitHub gets corrupted. Ask Brian to recompile it, or recompile it yourself.
* Genetic maps have a header or are not tab-separated. Four column (PLINK style) genetics maps should be tab-separated and headerless. 
* Genetic maps and VCFs have different prefixes, e.g. chr7 versus 7 in CHROM column
* You are using Snakemake 8, or too old a version of Snakemake 7.
* "Error: Could not find or load main class ibdends.IbdEndsMain". Re-download the java files with ``rm -r software/`` and then ``bash get-software.sh``.
* Java class file is too old. The Flare Jan 24 version worked for JDK11, but I upgraded to JDK23 to get Flare Oct 24 version to run.

Reproducing paper results
##############

The tag v1.0 is closest to the code used in our publications. The scripts in the tag to simulate data with msprime and capture the IBD segments with tskibd are used in the Temple and Browning (2025+) publication.

.. note::

   The simulation study in ``workflow/simulate`` was used in Temple, Waples, and Browning (2024). The scripts are older versions of this software. I will provide minimal/some support if one wants to replicate our results or use our SLiM simulation scripts.

.. note::

   The branch ``bring_clues_update`` has ``workflow/other-methods`` for the comparisons in Temple, Waples, and Browning (2024).
