
## What to do:

- Make an experiments file
  - Use experiments/setup-snakemake-experiment.py
  - python setup-snakemake-experiment.py
- Modify text files yamls/
  - See README in this subfolder

## How to perform simulation study:

1. `conda activate isweep`

2. `nohup snakemake -c1 -s Snakefile-simulate.smk --cluster "[options]" --jobs X --configfile *.yaml & `
 - I like to use the following snakemake options
    - `-n` to do a dry-run (what will be done)
    - `--keep-going`
    - `--rerun-triggers mtime` (if things didn't finish)
    - `--rerun-incomplete` (if things didn't finish)
    - `--latency-wait 300`
    - `nohup ... & ` let's you run pipeline in background

3. ` nohup snakemake -c1 -s Snakefile-ibd.smk --cluster "[options]" --jobs X --configfile *.yaml & `
  - For the qsub cluster manager, I pass `-V` into `--cluster` options and my active environment is isweep

4. Write some python/R scripts to summarize the data and make plots

The `Snakefile` solo does both `Snakefile-simulate.smk` + `Snakefile-ibd.smk`.

## How to clean up files:

- You should only do this once you're done / finalized with your work.
  - This includes running any comparisons you care about (isafe, clues, etc.)
- Linux command `rm`  
  - Use wildcards *
  - Don't remove these main data files:
    - large.chr1.vcf.gz
    - slimulation.trees.tsz
    - slimulation.vcf.gz or slimulation.bcf or slimulation.bcf.gz
    - results.hap.tsv
    - results.snp.tsv
    - ibd.entropy.tsv
  - MACRO in your yamls/ is the main folder directory (a string)
  - MICRO in your experiment file defines subfolders (a string)
  - REP in your experiment file defines subsubfolders (a number)

### Special for method developers

You can comment out `rule genotyping_error` under true phase, and comment in the three rules under infer phase. This change makes your simulation study based on genotypes phased by the Beagle software. `bash get-software.sh` to get `beagle.jar` and the `remove-phase.jar` programs.