# isweep 

incomplete Selective sweep With Extended haplotypes Estimation Procedure

<img src="isweep-icon.png" align="center" width="400px"/>

This software presents statistical methods to study very recent, very strong positive selection.
- By very recent, we mean within the last 500 generations
- By very strong, we mean selection coefficient s >= 0.015 (1.5%)

The methods relate the lengths of IBD tracts to a coalescent model under selection. 

We assume 1 selected allele at a biallelic locus.

There are ~ seven major ideas.

1. A genome-wide selection scan that looks for anomalously large IBD rates
2. Inferring anomalously large IBD clusters
3. Ranking alleles based on evidence for selection
4. Computing a measure of cluster agglomeration (an IBD information entropy)
5. Estimating allele frequency of unknown sweeping allele
6. Estimating selection coefficients (w/ nice statistical properties)
7. Estimating a confidence interval (w/ nice statistical properties)

These steps are implemented automatically in a `snakemake` pipeline.

## The input data is:

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

## Installation

1. ` git clone https://github.com/sdtemple/isweep.git `
2. ` conda env create -f isweep-environment.yml `
  - ` conda activate isweep `
  - ` python -c 'import site; print(site.getsitepackages())' `
3. Download other software
  - See below
  - Remember where you put these
  - Put them in the same place
  - You need to cite these software

## Computer Overview

This repository contains a Python package and some Snakemake bioinformatics pipelines.
- The package ---> src/
- The pipelines ---> workflow/

You should run all `snakemake` pipelines in their `workflow/some-pipeline/`.

You should be in the `conda activate isweep` environment for analyses.

You should run the analyses using cluster jobs.

We have made README.md files in most subfolders.

## Running the isweep procedure:

Phase data w/ Beagle or Shapeit beforehand.

1. Make pointers to large (phased) vcf files
2. Edit yaml files in the different workflow directories
3. Run the selection scan
- ` nohup snakemake -s Snakefile-scan.smk -c1 --cluster "[options]" --jobs X --configfile *.yaml & `
4. Estimate recent effective sizes
- ` workflow/scan/terminalscripts/run-ibdne.sh `
5. Make the Manhattan plot
- ` workflow/scan/terminalscripts/manhattan.py `
6. Checkout the roi.tsv file
  - ` cat roi.tsv `
  - Edit it with locus names if you want
7. Run the region of interest analysis
  - ` nohup snakemake -s Snakefile-roi.smk -c1 --cluster "[options]" --jobs X --configfile *.yaml & `

Refer to the README.md in ` cd workflow/some-pipeline ` for more instructions.

Tip: define variables for file, folder names, e.g., `variable1=1224 ` then `echo $variable1 `

## Other software

Put these in a common folder!

- Beagle (https://faculty.washington.edu/browning/beagle/beagle.html)
- hap-ibd.jar (https://github.com/browning-lab/hap-ibd)
- ibd-ends.jar (https://github.com/browning-lab/ibd-ends)
- filter-lines.jar (https://faculty.washington.edu/browning/beagle_utilities/utilities.html)
- ibdne.jar (https://faculty.washington.edu/browning/ibdne.html)
- IBDkin (https://github.com/YingZhou001/IBDkin)
- bcftools
- tabix
- gzip

bcftools, tabix, gzip are in the isweep conda environment.

If you perform simulation study, you require more software.
- SLiM (https://messerlab.org/slim/)
- add-uniform-err.jar (https://faculty.washington.edu/browning/beagle_utilities/utilities.html)

If you want to compare against other methods (using our pipelines), you require more software.
See `workflow/other-methods/` folder.

## Citation

Please cite if you use this package.

Temple, S.D., Waples, R.K., Browning, S.R. (2023) "Modeling recent positive selection in Americans of European ancestry"
https://www.biorxiv.org/content/10.1101/2023.11.13.566947v1

## Development things to do

- Provide some scripts to summarize analyses
- Give example commands, scripts for pre-processing steps
- Change all snakescripts/ to terminalscripts/
- Severely simplify the yaml files
  - So many parameters for methods development, not user-friendly
