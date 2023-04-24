# isweep <img src="isweep-icon.png" align="right" width="300px"/>

incomplete Selective sweep With Extended haplotypes Estimation Procedure

This software presents a statistical methodology for inferring the selection coefficient and the allele frequency of a single adaptive allele introduced to a population in the recent 50-600 generations. It uses identity-by-descent segments, i.e., long haplotypes, to make these inferences. The methodology proposes a parametric statistical model that relates segment data to population genetics and coalescent trees. Please read the following paper for more context:

### Software Overview

This repository contains a Python package and some Snakemake bioinformatics pipelines to conduct simulation experiments and real data analysis. Functions in the Python package have been documented and some functions have been unit tested or compared against other methods (e.g. tskit). The main package contains the subpackage, cis (coalescent identity-by-descent segments), which is for studying and simulating identity-by-descent segments with respect to a coalescent process.

### Installation

iSWEEP Python package
1. git clone https://github.com/sdtemple/iSWEEP.git
2. python -c 'import site; print(site.getsitepackages())'
3. cp -r src/iSWEEP/ "LOCATION OF YOUR PYTHON PACKAGES PRINTED ABOVE"

iSWEEP dependencies
* Use conda or pip to install/update numpy, pandas, tskit, stdpopsim, networkx, matplotlib
* Install SLiM and put a slim executable in the folder `software/` (https://messerlab.org/slim/)
* Install or update Browning Lab (http://faculty.washington.edu/browning/) JAR files
* Some JAR files already in `software/`
* Install bcftools, tabix, gzip

PARAGRAPH TO DESCRIBE INSTALATION DETAILS

### User Guide

See the file scripts/vignette.ipynb.

Procedure 1 (study true IBD from coalescent)

A. Simulate IBD segments and IBD communities from structured coalescent \
B. Use IBD communities to infer allele frequency \
C. Use IBD segments to infer selection coefficient (conditional on allele frequency) \
D. Parametric boostrap given (s,p) \
E. Summarize results

Procedure 2 (study inferred IBD from simulated sequence data)

A. Create or modify files in env/*
B. Simulate sequence data with tskit and msprime \
C. Infer IBD segments and IBD communities with JAR programs \
D. Use IBD communities to infer allele frequency \
E. Use IBD segments to infer selection coefficient (conditional on allele frequency) \
F. Parametric boostrap given (s,p) \
G. Summarize results

### Developer Guide

See DETAILS.md.

- documents/* : written references
- env/* : specifications for Snakemake runs
- scripts/* : analysis files 
- src/* : location of Python packages
- snakescripts/* : modified scripts for Snakemake
- workflow/* : place for Snakemake file

src/iSWEEP
  - iSWEEP/ 
    - isweep.py 
      - Computing chi-squared statistics (estimate selection coefficient) 
      - Bootstrapping 
      - Performance metrics 
    - isweepUtilities.py 
      - Putting IBD segments from .ibd file into bins
    - diameter.py
      - Greedy algorithm to infer 2K diameter communities
      - 3 or 5 sigma rule outlier detection  
    - cis/ 
      - coalescentIBD.py [MOST IMPORTANT FILE] 
        - Perform Wright Fisher simulations 
        - Simulate coalescent processes 
        - Simulate IBD segments from coalescent process 
        - Study and simulate quasi-geometric and Gamma(2,.) random variables 
      - cisUtilities.py 
        - Functions to create, modify .ne files 
        - Plotting functions 

### Statistical Methodology

PARAGRAPH TO DESCRIBE THE METHODOLOGIES AVAILABLE
