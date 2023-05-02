# isweep <img src="isweep-icon.png" align="right" width="200px"/>

incomplete Selective sweep With Extended haplotypes Estimation Procedure

This software presents a statistical methodology for inferring the selection coefficient and the allele frequency of a single adaptive allele introduced to a population in the recent 50-600 generations. It uses identity-by-descent segments, i.e., long haplotypes, to make these inferences. The methodology proposes a parametric statistical model that relates segment data to population genetics and coalescent trees. Please read the following paper for more context:

### Software Overview

This repository contains a Python package and some Snakemake bioinformatics pipelines to conduct simulation experiments and real data analysis. Functions in the Python package have been documented and some functions have been unit tested or compared against other methods (e.g. tskit). The main package contains the subpackage, cis (coalescent identity-by-descent segments), which is for studying and simulating identity-by-descent segments with respect to a coalescent process.

### Pipeline

1. Set up a conda environment with the environment.yml file in workflow/
  * conda create env -f environment.yml
2. Install all required programs (and bcftools, tabix)
  * See documentation below
3. Initial steps
  * Run maps.sh
    * First arg, MAPS, is the folder where your genetic maps are
    * Second arg, MAPPREFIX, is the prefix to your map files
    * Third arg, MAPSUFFIX, is the suffix to your map files
    * Fourth arg, WHERE, should be where you will do your analysis
    * Fifth arg, CHRLOW, is the min chromosome number
    * Sixth arg, CHRHIGH, is the max chromosome number
  * Check relatedness, population structure and write and direct to subsample of interest
    * Use PCA or ADMIXTURE-like to study population structure
    * There are various tools to infer relatedness in genetic samples
    * (Not advised) An analysis with mainly closely related genetic samples
4. Navigate to workflow/scan/ and run the Snakemake
  * Modify the .yaml file for your analysis (the CHANGE parts)
    * Note directing to where .vcf.gz are
  * Check: snakemake all -c1 -n --configfile .yaml
  * Activate conda environment: conda activate isweep
  * On cluster: nohup snakemake all -c1 --jobs 100 --latency-wait 120 --keep-going --cluster "" &
    * nohup snakemake all -c1 --jobs 100 --latency-wait 120 --keep-going --cluster "qsub -q b-students.q -N {rule} -e ~/logs/{rule}.e -o ~/logs/{rule}.o -V -m e -M sdtemple@uw.edu -pe local 16" &
      * -pe local 16 will assign 16 cores to cluster jobs
      * -V is important to pass your conda environment in
      * -N , -e , -o , -m e -M sdtemple@uw.edu are not necessary
  * If necessary, modify memory resources in rules/scan.smk
5. Look at roi.tsv , plots/ , stats/
  * Rename and modify roi.tsv as desired
  * This is you deciding which regions of interest
  * Add a few Mb to left and right of BPLEFT and BPRIGHT !!!
6. Navigate to workflow/roi/ and run the Snakemake 
  * Same instructions as step 4
7. You can run manhattan-rates.R in scan/terminalscripts/ to make Manhattan plots
  * head -n 20 manhattan-rates.R
  * Rscript --vanilla manhattan-rates.R args

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
