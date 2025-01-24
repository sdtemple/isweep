Welcome to isweeps's documentation!
===================================

**isweep** is a Python package and a series of automated workflows to study natural selection with identity-by-descent (IBD) segments.

The name is *i*ncomplete *S*elective sweep *W* *E*xtended haplotypes *E*stimation *P*rocedure.

The Python package simulates IBD segments around a locus and estimates selection coefficients.

The automated workflows are:
- `workflow/scan-selection`: detects selected loci with rigorous multiple testing thresholds
- `workflow/model-selection`: estimates the location, allele frequency, and selection coefficient of a sweep
- `workflow/scan-case-control`: detects loci where IBD rates differ between binary cases and controls
- `workflow/prepare`: supports haplotype phasing, local ancestry, and kinship inference

Check out the :doc:`usage` section for further information, including
how to :ref:`installation` the project.

Installation
--------

.. code-block:: bash

    git clone https://github.com/sdtemple/isweep.git
    mamba env create -f isweep-environment.yml
    bash get-software.sh

.. note::

   This project is under active development.

Contents
--------

.. toctree::

   usage
   api
