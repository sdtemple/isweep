<img src="isweep-icon.png" align="center" width="400px"/>

CAUTION: I did some re-organizing. Testing them out, debugging in late November.
- Need to check the workflow/roi/ and workflow/scan/

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
- Check the *.log files from ibd-ends. If it recommends an estimated err, change the error rate in yaml.
4. Estimate recent effective sizes (workflow/scan)
- ` workflow/scan/terminalscripts/run-ibdne.sh `
5. Make the Manhattan plot (workflow/scan)
- ` workflow/scan/terminalscripts/manhattan.py `
6. Checkout the roi.tsv file
  - ` cat roi.tsv `
  - Edit it with locus names if you want
7. Run the region of interest analysis (workflow/roi)
  - ` nohup snakemake -s Snakefile-roi.smk -c1 --cluster "[options]" --jobs X --configfile *.yaml & `

Tip: define variables for file, folder names, e.g., `variable1=1224 ` then `echo $variable1 `

## Do these methods apply to my sample size?

It'd be nice if you have:
- Ancestral Ne >= 3000
- No population bottleneck with Ne < 3000 in last 500 generations
- Sample from 1 generation
  - I.e., use other methods for flies, mosquitoes, worms

#### You could test if there is enough IBD:

- In python, use the `simulate_ibd` and `read_Ne` functions in `isweep`
  - Ne file is tab-separated w/ two columns: generation, Ne
  - You should make at least 500 generations
  - It is okay to have a guess of Ne size before inferring with IBDNe
    - Use the different *_Ne() functions to make such an *.ne file
- Do you have on average # of 3.0 cM ibd tracts >= 200?
  - Then, yes, you probably have enough samples to run this analysis.

#### Some species from `stdpopsim` that may apply:

- Great apes
  - Pongo abelii
  - Pan troglodytes
  - Papio anubis
- Anas platyrhynchos (duck)
- Apis mellifera (bee)
- Canis familiaris (dog)
  - Address population structure, admixture

#### Should I use these methods on biobank datasets > 100k individuals?

- No, unless you want to spend a lot on computing costs.
  - The methods should scale, but you'd wait awhile and spend a lot.
- If you downsample to ~10k, you should have more than enough data.
 - (If IBD **were** binomially distributed, downsampling keeps same distributional properties.)

## Development things to do

- Provide some scripts to summarize analyses
- Test performance in array data, less dense sequence data
- Test performance when phase is inferred with beagle
- Test performance with no gene conversion
- Not designed for ploidy != 2 (yet)
- Debug population structure simulation script
  - Appears to fail when coming together is later than *.ne file
- Debug when there is no excess IBD
