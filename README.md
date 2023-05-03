# isweep <img src="isweep-icon.png" align="right" width="200px"/>

incomplete Selective sweep With Extended haplotypes Estimation Procedure

This software presents a statistical method for inferring the selection coefficient and the allele frequency of a single adaptive allele introduced to a population in the recent 50-600 generations. It uses identity-by-descent segments, i.e., pairwise matching long haplotypes, to make these inferences. The methodology proposes a parametric statistical model that relates segment data to population genetics and coalescent trees. Please read the following paper for more context:

## Computer Overview

This repository contains a Python package and some Snakemake bioinformatics pipelines to conduct simulation studies and real data analysis. Functions in the Python package have been documented and tested. The main package contains the subpackage, cis (coalescent identity-by-descent segments), which is for studying and simulating identity-by-descent segments with respect to a coalescent process.

### Installation

isweep real data analysis
1. git close https://github.com/sdtemple/isweep.git
2. pip install the package
3. Check installation
  * python -c 'import site; print(site.getsitepackages())'
4. Download other software
5. Point to this software in your .yaml

isweep simulation study
1. git clone https://github.com/sdtemple/isweep.git
2. conda create env -f environment.yml
  * In simulate/
3. conda activate isweep-simulate
4. Check installation
  * python -c 'import site; print(site.getsitepackages())'
  * python -c 'import tskit; import msprime; import pyslim'
5. Download other software
6. Point to this software in your .yaml

Software dependencies
* Install SLiM and put a slim executable in the folder `software/` (https://messerlab.org/slim/)
* Install or update Browning Lab (http://faculty.washington.edu/browning/) JAR files
* Install bcftools, tabix, gzip

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
  * Check relatedness, population structure and write and direct to subsample of interest (see considerations section)
    * Use PCA or ADMIXTURE-like to study population structure
    * There are various tools to infer relatedness in genetic samples
    * (Not advised) An analysis with mainly closely related genetic samples
  * If you have very few samples relative to a meta biobank vcf, consider using bcftools view -S
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
5. Run IBDNe software (http://faculty.washington.edu/browning/)
 * zcat chr*.ibd.gz | java -Xmx100g -jar ibdne.jar map= out=
 * Do this on a cluster !!!
 * -Xmx100g is the heap memory
 * map= is a genetic map where you concatenate all chromosomes
 * out= is prefix for your .ne file
 * Point to your .ne file in step 7
 * You may need to downsample for this
   * isweep only requires enough samples to have uniformly non-zero ibd segment counts
   * Make an issue if you run out of memory
6. Look at roi.tsv , plots/ , stats/
  * Rename and modify roi.tsv as desired
  * This is you deciding which regions of interest
  * Add a few Mb to left and right of BPLEFT and BPRIGHT !!!
  * See considerations section
7. Navigate to workflow/roi/ and run the Snakemake
  * Same instructions as step 4
8. You can run manhattan-rates.R in scan/terminalscripts/ to make Manhattan plots
  * head -n 20 manhattan-rates.R
  * Rscript --vanilla manhattan-rates.R args
  * Assuming you have a version of R programming

### considerations

* Something about selecting the finalroi.tsv
* Something about beagle phasing
* Something about PCA
* Something about IBDNe

## Statistical Methods

The premise of isweep is to estimate the selection coefficient and the current adaptive allele frequency p(0). Together these imply a trajectory back in time of the adaptive allele. The Wright-Fisherian trajectory back assumes de novo sweep, ongoing sweep, known Mendelian inheritance, and known Ne(t) for past 150 or so generations. The method of moments estimator is also conditional, namely s | p(0), Ne(t), standing variation, sweep end, Mendelian inheritance. The method of moments matches expected IBD segment counts to observed IBD segment counts. We estimate p(0) with a weighted allele frequency estimator where weights come from a novel community detection algorithm (see our paper). You can estimate Ne(t) with IBDNe. You can change default settings to allow sweeps from standing variation or sweeps that recently ended, if this is supported by literature. Below we provide a high level look at three procedures in workflow/ to study isweep.

Procedure 1 (study true IBD from coalescent; coalescent/ workflow)

A. Simulate IBD segments and IBD communities from structured coalescent \
B. Use IBD communities to infer allele frequency \
C. Use IBD segments to infer selection coefficient (conditional on allele frequency) \
D. Parametric boostrap given (s,p) \
E. Summarize results

Procedure 2 (study inferred IBD from simulated sequence data; simulate/ workflow)

A. Simulate sequence data with SLiM, tskit, and msprime \
B. Manage .vcf.gz with bcftools \
C. Infer IBD segments and IBD communities with JAR programs \
D. Use IBD communities to infer allele frequency \
E. Use IBD segments to infer selection coefficient (conditional on allele frequency) \
F. Parametric boostrap given (s,p) \
G. Summarize results
  * For comparison, this may involve iSAFE and CLUES
  * We find that iSAFE downgrades bcftools, so do separately

## References

## Developer Guide

See DETAILS.md.

- documents/* : written references
- env/* : specifications for Snakemake runs
- scripts/* : analysis files
- src/* : location of Python packages
- snakescripts/* : modified scripts for Snakemake
- workflow/* : place for Snakemake file

src/
  - isweep/
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
