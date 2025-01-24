Welcome to isweeps's documentation!
===================================

**isweep** is a Python package and a series of automated workflows to study natural selection with identity-by-descent (IBD) segments.

The name is "incomplete **S**elective sweep **W** **E**xtended haplotypes **E**stimation **P**rocedure.

The Python package simulates IBD segments around a locus and estimates selection coefficients.

The automated workflows are:
- `workflow/scan-selection`: detects selected loci with rigorous multiple testing thresholds
- `workflow/model-selection`: estimates the location, allele frequency, and selection coefficient of a sweep
- `workflow/scan-case-control`: detects loci where IBD rates differ between binary cases and controls
- `workflow/prepare`: supports haplotype phasing, local ancestry, and kinship inference



**Lumache** (/lu'make/) is a Python library for cooks and food lovers
that creates recipes mixing random ingredients.
It pulls data from the `Open Food Facts database <https://world.openfoodfacts.org/>`_
and offers a *simple* and *intuitive* API.

Check out the :doc:`usage` section for further information, including
how to :ref:`installation` the project.

.. note::

   This project is under active development.

Contents
--------

.. toctree::

   usage
   api
