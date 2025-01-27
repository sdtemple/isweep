Documentation
===================================

.. image:: images/isweep-icon.png
   :align: center
   :width: 600px

|
|

**isweep** is a Python package and a series of automated workflows to study natural selection with identity-by-descent (IBD) segments. The Python package simulates IBD segments around a locus and estimates selection coefficients. The automated workflows are:

* :ref:`selection-scan`: detects selected loci with rigorous multiple testing thresholds
* :ref:`modeling-hard-sweeps`: estimates location, allele frequency, and selection coefficient of sweep
* :ref:`case-control-scan`: detects loci where IBD rates differ between binary cases and controls
* :ref:`phasing-and-ancestry`: supports haplotype phasing, local ancestry, and kinship inference

Each automated workflow has a dedicated page under :doc:`usage`. The general way to run these methods is:

1. Navigate to the appropriate workflow directory
2. Modify parameters in YAML configuration files
3. Send the jobs to a cluster with ``nohup snakemake [...] &``

I made a Zenodo repository with some simulated data to test the workflows. See :ref:`testing-workflows`.

The source code is `here <https://github.com/sdtemple/isweep>`_

Installation
--------

.. code-block:: shell

    git clone https://github.com/sdtemple/isweep.git
    mamba env create -f isweep-environment.yml
    bash get-software.sh
    pip install isweep

.. note::

   
   This project is in a stable state. I am commited to providing quick support via GitHub Issues at least into 2026.


Data Requirements
--------

The main requirements are a tab-separated genetic recombination map and enough samples to detect more than 0 IBD segments at all positions.

* Phased haplotypes (VCF files)
* Samples with a similar ancestry
* No close relatives
* Recombining autosomes

If you don't already have your data phased or cohort selected, we support :ref:`phasing-and-ancestry` with another workflow.

In humans, more than 1000 is enough samples, but more than 3000 samples is recommended. 

At some point, there are no gains in statistical power with more samples. I do not recommend analyzing more than 100k samples in a biobank. See our Temple and Browning (2025+) publication.

The tree of life is messy. Email or make a `GitHub Issue <https://github.com/sdtemple/isweep/issues>`_ for analysis advice about the nuances of your sample population.


Vignette
--------

Outside of the running the workflows, the main functions are:

* ``isweep.coalescent.simulate_ibd_isweep``: generate long IBD segments around a locus (w/ selection)
* ``isweep.coalescent.chi2_isweep``: use a uniroot finder with this to estimate the selection coefficient
* ``isweep.utilities.read_Ne``: load in recent effective population sizes

For instance:

.. code-block:: python

   from isweep import *
   # parameter settings
   s = 0.03
   p=float(0.5) # allele freq
   Ne=read_Ne('constant-10k.ne') # demo history
   model='m'
   long_ibd=3.0
   ab=[long_ibd,np.inf]
   nsamples=200
   # calculate denominator
   ploidy=2
   msamples=ploidy*nsamples
   N=msamples*(msamples-1)/2-msamples
   # simulate data
   out=simulate_ibd_isweep(
    nsamples,
    s,
    p,
    Ne,
    long_ibd,
    long_ibd,
    one_step_model=model,
    ploidy=ploidy,
   )
   ibd=out[0][0]
   # estimating the selection coefficient
   se = minimize_scalar(
    chi2_isweep,
    args=(p,Ne,N,(ibd,),ab,model,0,-0.01,ploidy),
    bounds=(0,0.5),
    method='bounded'
   ).x
   print('true selection coefficient')
   print(s)
   print('estimate selection coefficient')
   print(se)

There are also some IPython notebooks as examples in ``vignettes/``.

Citations
--------

This software and its methods are the basis of six publications.

* `Modeling hard sweeps <https://www.sciencedirect.com/science/article/pii/S0002929724003331>`_
* `Simulating IBD segments <https://www.biorxiv.org/content/10.1101/2024.12.13.628449v2>`_
* `Central limit theorems <https://www.biorxiv.org/content/10.1101/2024.06.05.597656v2>`_
* `Thesis on recent positive selection <https://www.proquest.com/docview/3105584569>`_
* TBD (Multiple testing in selection scan)
* TBD (Multiple testing in case-control scan)

The software `Beagle <https://faculty.washington.edu/browning/beagle/beagle.html>`_, `ibd-ends <https://github.com/browning-lab/ibd-ends/>`_, `hap-ibd <https://github.com/browning-lab/hap-ibd>`_, and `flare <https://github.com/browning-lab/flare>`_ are also used and should be cited.

* ``workflow/scan-selection``: hap-ibd and ibd-ends
* ``workflow/scan-case-control``: hap-ibd and ibd-ends
* ``workflow/model-selection``: hap-ibd
* ``workflow/phasing-ancestry``: Beagle, flare, and hap-ibd


Contents
--------

.. toctree::

   usage

API
--------
.. toctree::
   
   api


Contact
--------

Seth Temple (sethtem@umich.edu) or `GitHub Issues <https://github.com/sdtemple/isweep/issues>`_