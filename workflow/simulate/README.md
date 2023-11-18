
### What to do

- Download the required software (if you haven't already) and put it in a folder location
  - See bottom of page
- Make an experiments file
  - Use experiments/setup-snakemake-experiment.py
  - python setup-snakemake-experiment.py
- Modify text files yamls/
  - See README in this subfolder

### How to perform simulation study

1. nohup snakemake -c1 -s Snakefile-simulate.smk --cluster "[options]" --jobs X --configfile *.yaml &
  - I like to use the following snakemake options
    - "-n" to do a dry-run (what will be done)
    - "--keep-going"
    - "--rerun-triggers mtime"
    - "--rerun-incomplete"
    - "--latency-wait 300"
    - "nohup ... &" let's you run pipeline in background
2. nohup snakemake -c1 -s Snakefile-ibd.smk --cluster "[options]" --jobs X --configfile *.yaml &
  - For the qsub cluster manager, I pass "-V" into --cluster options and my activate environment is isweep
3. Write some python/R scripts to summarize the data and make plots

### How to clean up files

- You should only do this once you're done / finalized with your work.
  - This includes running any comparisons you care about (isafe, clues, etc.)
- Linux command "rm"  
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

### Required software
- https://messerlab.org/slim/
- https://github.com/browning-lab/hap-ibd
- https://github.com/browning-lab/ibd-ends
- http://faculty.washington.edu/browning/add-uniform-err.jar
- http://faculty.washington.edu/browning/filter-lines.jar