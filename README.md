<img src="isweep-icon.png" align="center" width="600px"/>

See `announcements.md` for high-level updates on this repo.

See `usage.md` to evaluate if this methodology fits your study.

## Citation

Please cite if you use this package.

Temple, S.D., Waples, R.K., Browning, S.R. (2023) "Modeling recent positive selection in Americans of European ancestry"
https://www.biorxiv.org/content/10.1101/2023.11.13.566947v1

## Methodology

Acronym: Incomplete Selective sweep With Extended haplotypes Estimation Procedure

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

## Installation

1. ` git clone https://github.com/sdtemple/isweep.git `
2. ` conda env create -f isweep-environment.yml `
  - ` conda activate isweep `
  - ` python -c 'import site; print(site.getsitepackages())' `
3. Download software
  - You need to cite these software
  - ` bash get-software.sh software `
    - Puts these in a folder called `software/`
    - Requires wget
  - For simulation study, download SLiM yourself
    - Put in `software/`
    - https://messerlab.org/slim/

Use IBDkin (https://github.com/YingZhou001/IBDkin) to remove close relatives.

Use PCA or ADMIXTURE to subset based on ancestry.

If you want to compare against other methods (using our pipelines), you require more software.

See `workflow/other-methods/` folder.
  
## Overview

This repository contains a Python package and some Snakemake bioinformatics pipelines.
- The package ---> src/
- The pipelines ---> workflow/

You should run all `snakemake` pipelines in their `workflow/some-pipeline/`.

You should be in the `conda activate isweep` environment for analyses.

You should run the analyses using cluster jobs.

We have made README.md files in most subfolders.

## The input data is:

- Whole genome sequences
  - Probably at least > 500 diploids
  - Phased vcf data 0|1
  - No apparent population structure
  - No apparent close relatedness
  - A genetic map (bp ---> cM)
  - Recombining diploid chromosomes
    - Not extended to human X chromosome (yet?)
- Access to cluster computing
  - You should have at least 25 Gb of RAM and 6 CPUs on a node
    - More for larger datasets
  - Have not extended to cloud computing (yet?)

## Running the procedure:

This is the overall procedure. You will see more details for each step in workflow/some-pipeline/README.md files.

Phase data w/ Beagle or Shapeit beforehand.

1. Make pointers to large (phased) vcf files
2. Edit yaml files in the different workflow directories
3. Run the selection scan (workflow/scan)
- ` nohup snakemake -s Snakefile-scan.smk -c1 --cluster "[options]" --jobs X --configfile *.yaml & `
- Recommendation: Do a test run with your 2 smallest chromosomes.
- Check the *.log files from ibd-ends. If it recommends an estimated err, change the error rate in yaml.
- Then, run with all your chromosomes.
4. Estimate recent effective sizes (workflow/scan)
- ` workflow/scan/scripts/run-ibdne.sh `
5. Make the Manhattan plot (workflow/scan)
- ` workflow/scan/scripts/manhattan.py `
6. Checkout the roi.tsv file
  - ` cat roi.tsv `
  - Edit it with locus names if you want
7. Run the region of interest analysis (workflow/roi)
  - ` nohup snakemake -s Snakefile-roi.smk -c1 --cluster "[options]" --jobs X --configfile *.yaml & `

Tip: define variables for file, folder names, e.g., `variable1=1224 ` then `echo $variable1 `

## Considerations

- Positive selection
  - s >= 0 for estimation (use bounds in scipy.optimize.minimize_scalar) 

## Development things to do

- Provide option for model selection, standing variation, etc. in `roi.tsv` 
- Provide some scripts to summarize analyses
- Test performance in array data, less dense sequence data
- Test performance in dominance selection (sequence data)
- Time varying selection coefficients?
- Not designed for ploidy != 2 (yet)

## Current bugs
