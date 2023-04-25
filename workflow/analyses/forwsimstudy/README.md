What to do
  * Download Browning Lab software and put in same place as yaml definition
  * Copy and modify simstudies.reference.yaml
  * Or, use simstudies.*.yaml
  * Run setup-snakemake-experiment.py to get define replicate studies
    * Modify your yaml file accordingly

Terminal commands
  * snakemake -c1 -n
    * Dry run of bioinformatics pipeline
  * snakemake -c2
    * Run with 2 cores
  * yada
    * Version for a computing cluster

How to clean up
  * Only do this after you have all desired results
    * This means you may still want to compare to CLUES, iSAFE
  * MACRO is the macro folder
  * MICRO is the micro folder
  * REP is a number
  * rm MACRO/MICRO/REP/
    * Any vcf.gz
    * Any ibd.gz
    * Any hbd.gz
    * Any .tbi
    * Any .csi
    * Any .bgz
    * Any firstpass.

What to keep
  * .slim files
  * .log files
  * .tsz files
    * Use this to conduct array-like studies
  * .tsv files
  * .tsv.gz files
  * arguments.yaml
