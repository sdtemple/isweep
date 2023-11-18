# isweep 

CAUTION: THIS PACKAGE IS ACTIVELY BEING REORGANIZED, CLEANED FOR USER-FRIENDLINESS. 
CHECK BACK IN EARLY DECEMBER 2023.

incomplete Selective sweep With Extended haplotypes Estimation Procedure

<img src="isweep-icon.png" align="center" width="400px"/>

This software presents statistical methods to study very recent, very strong positive selection.
- By very recent, we mean within the last 500 generations
- By very strong, we mean selection coefficient s >= 0.015 (1.5%)

#### The methods relate the lengths of IBD tracts to a coalescent model under selection. We assume 1 selected allele at a biallelic locus.
#### There are ~ six major ideas.

1. A genome-wide selection scan that looks for anomalously large IBD rates
2. Inferring anomalously large IBD clusters
3. Ranking alleles based on evidence that they are selected or in strong correlation w/ a selected allele.
4. Computing a measure of cluster agglomeration (an IBD information entropy)
5. Estimating an allele frequency of a sweeping allele, w/o knowing its identity (it also doesn't have to be genotyped)
6. Estimating selection coefficients (w/ nice statistical properties)
7. Estimating a confidence interval (w/ nice statistical properties)

(I plan to create some isolated bash scripts for these.
Currently, they are built in snakemake pipelines.)

#### The input data is:

- Whole genome sequences
  - Probably at least > 500 diploids
  - Phased data
    - Use Beagle (https://faculty.washington.edu/browning/beagle/beagle.html)
  - No apparent population structure (i.e., admixture, diverged groups in same dataset)
  - No apparent close relatedness
    - Use IBDkin (https://github.com/YingZhou001/IBDkin)
    - Or something else like plink (https://zzz.bwh.harvard.edu/plink/)
- Access to cluster computing
  - You should have at least 25 Gb of RAM and 6 CPUs on a node
    - More for larger datasets
  - Have not extended to cloud computing yet

#### Running the isweep procedure:

Phase data w/ Beagle or Shapeit beforehand.

1. Edit yaml files in the different workflow directories
2. Make pointers to large vcf files
3. nohup snakemake -s Snakefile-scan.smk -c1 --cluster "[options]" --jobs X --configfile *.yaml &
4. Run ibdne scripts
5. Make manhattan plot
6. Checkout the roi.tsv table
  - Edit it with locus names if you want
7. nohup snakemake -s Snakefile-roi.smk -c1 --cluster "[options]" --jobs X --configfile *.yaml &
  - Evaluate your analysis

Refer to the README.md in "cd workflow/*"

#### Computer Overview

This repository contains a Python package and some Snakemake bioinformatics pipelines.
- The package ---> src/
- The pipelines ---> workflow/

We have made README.md files in most subfolders.

#### Installation

1. git clone https://github.com/sdtemple/isweep.git
2. conda env create -f workflow/environment.yml
  - conda activate isweep
  - python -c 'import site; print(site.getsitepackages())'
3. Download other software
  - See below
  - Remember where you put these
  - Put them in the same place
  - You need to cite these software

#### Other software

- Beagle
- hap-ibd.jar
- ibd-ends.jar
- filter-lines.jar
- ibdne.jar
- bcftools
- tabix
- gzip

bcftools, tabix, gzip are in isweep conda environment.

If you perform simulation study, you require more software.
- SLiM
- add-uniform-err.jar

If you want to compare against other methods (using our pipelines), you require more software.
See workflow/other-methods/ folder.


### Citation

Please cite if you use this package
Temple, S.D., Waples, R.K., Browning, S.R. (2023) "Modeling recent positive selection in Americans of European ancestry"
https://www.biorxiv.org/content/10.1101/2023.11.13.566947v1
