# Statistical inference using long IBD segments

<img src="isweep-icon.png" align="center" width="600px"/>

New features are actively under construction (Fall 2024).
- cM, not bp, windowing
- Integrate multiple testing correction into pipeline
- Extension to IBD mapping

Contact sethtem@umich.edu or Github issues for troubleshooting.

See `misc/usage.md` to evaluate if this methodology fits your study.

See `misc/cluster-options.md` for some suggested cluster options to use in pipelines.

See on GitHub "Issues/Closed" for some comments **I/Seth** left about the pipeline. 

## Citation

Please cite if you use this package.

#### Methods to model selection:

Temple, S.D., Waples, R.K., Browning, S.R. (2024). Modeling recent positive selection using identity-by-descent segments. The American Journal of Human Genetics. <a href="https://doi.org/10.1016/j.ajhg.2024.08.023">https://doi.org/10.1016/j.ajhg.2024.08.023</a>. 

#### Methods to simulate IBD segments and our central limit theorems:

Temple, S.D., Thompson, E.A. (2024). Identity-by-descent in large samples. Preprint at bioRxiv, 2024.06.05.597656. <a href="https://www.biorxiv.org/content/10.1101/2024.06.05.597656v1">https://www.biorxiv.org/content/10.1101/2024.06.05.597656v1</a>.

#### Multiple testing correction for selection scan

Temple, S.D. (2024). "Statistical Inference using Identity-by-Descent Segments: Perspectives on Recent Positive Selection. PhD thesis (University of Washington). <a href="https://www.proquest.com/docview/3105584569?sourcetype=Dissertations%20&%20Theses">https://www.proquest.com/docview/3105584569?sourcetype=Dissertations%20&%20Theses</a>.


## Methodology

Acronym: *i*ncomplete *S*elective sweep *W*ith *E*xtended haplotypes *E*stimation *P*rocedure

This software presents methods to study recent, strong positive selection.
- By recent, we mean within the last 500 generations
- By strong, we mean selection coefficient s >= 0.015 (1.5%) 

In modeling a sweep, we assume 1 selected allele at a locus.

### Automated analysis pipeline(s):

1. A genome-wide selection scan for anomalously large IBD rates
 - With multiple testing correction
2. Inferring anomalously large IBD clusters
3. Ranking alleles based on evidence for selection
4. Computing a measure of cluster agglomeration (Gini impurity index)
5. Estimating frequency and location of unknown sweeping allele
6. Estimating a selection coefficient
7. Estimating a confidence interval

Step 1 may be standalone, depending on the analysis. (You may not care to model putative sweeps (Steps 2-7).)

### The input data is:

See `misc/usage.md`.

- Whole genome sequences
  - Probably at least > 500 diploids
  - Phased vcf data 0|1
  - No apparent population structure
  - No apparent close relatedness
  - Tab-separated genetic map (bp ---> cM) 
  - Recombining diploid autosomes
- Access to cluster computing
  - Not extended to cloud computing

Chromosome numbers in genetic maps should match chromosome numbers in VCFs.

## Repository overview

This repository contains a Python package and some Snakemake bioinformatics pipelines.
- The package ---> `src/`
- The pipelines ---> `workflow/`

You should run all `snakemake` pipelines in their `workflow/some-pipeline/`.

You should be in the ```mamba activate isweep``` environment for analyses.

You should run the analyses using cluster jobs.

## Installation

See `misc/installing-mamba.md` to get a Python package manager.

1. Clone the repository
``` 
git clone https://github.com/sdtemple/isweep.git 
```
2. Get the Python package
``` 
mamba env create -f isweep-environment.yml
```
```
mamba activate isweep
```
```
python -c 'import site; print(site.getsitepackages())'
```
3. Download software.
``` 
bash get-software.sh 
```
  - Requires `wget`.
  - You need to cite these software.

## Pre-processing

Phase data w/ Beagle or Shapeit beforehand.
Subset data in light of global ancestry and close relatedness.
Example scripts are in `scripts/pre-processing/`.
- Here is a pipeline we built for these purposes: `https://github.com/sdtemple/flare-pipeline`
- You could use IBDkin to detect close relatedness: `https://github.com/YingZhou001/IBDkin`
- You could use PCA, ADMIXTURE, or FLARE to determine global ancestry. 

## Main analysis

You will see more details for each step in `workflow/some-pipeline/README.md` files.

### For all workflows
1. Make pointers to large (phased) vcf files.
2. Edit YAML files in the different workflow directories.

### Detecting recent selection

Run the selection scan (`workflow/scan-selection`).
``` 
nohup snakemake -s Snakefile-scan.smk -c1 --cluster "[options]" --jobs X --configfile *.yaml & 
```
  - See the file `misc/cluster-options.md` for support.
- Recommendation: do a test run with your 2 smallest chromosomes.
- Check `*.log` files from `ibd-ends`. If it recommends an estimated err, change error rate in YAML file.
- Then, run with all your chromosomes.

Make the Manhattan plot: ` workflow/scan-selection/scripts/manhattan.py `.

### Modeling putative sweeps

1. Estimate recent effective sizes :` workflow/scan-selection/scripts/run-ibdne.sh `.
2. Checkout the `roi.tsv` file.
  - Edit with locus names if you want.
  - Edit to change defaults: additive model and 95% confidence intervals.
3. Run the region of interest analysis (`workflow/model-selection`).
``` 
nohup snakemake -s Snakefile-roi.smk -c1 --cluster "[options]" --jobs X --configfile *.yaml & 
``` 

The script to estimate recent Ne can be replaced with any method to estimate recent Ne, as it happens before the `snakemake` command. This method [HapNe](https://palamaralab.github.io/software/hapne/) is one such option.

## Picture of selection scan workflow

The flow chart below shows the steps ("rules") in the selection scan pipeline.

Diverting paths "mle" versus "scan" refer to different detection thresholds (3.0 and 2.0 cM).

See `dag-roi.png` for the steps in the sweep modeling pipeline.

<img src="dag-scan.png" align="center" width="600px"/>

<!-- ## Picture of selection modeling workflow

<img src="dag-roi.png" align="center" width="600px"/> -->
