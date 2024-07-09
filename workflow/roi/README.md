# Second pipeline in data analysis:

## The procedure

1. Edit the YAML file.
2. ` mamba activate isweep `
3. Dry-run of the pipeline. Check what rules will run.
```
snakemake -c1 -s Snakefile-roi.smk --configfile *.yaml 
```
- ` -n ` option is the dry run
4. Run the pipeline for real.
```
nohup snakemake -c1 -s Snakefile-roi.smk --cluster "[options]" --jobs X --configfile *.yaml & 
```
- See `misc/cluster-options.md` to choose SLURM or QSUB options.
- Make sure to pass in the isweep environment to cluster.

<!-- ## Example

I use the qsub system to manage cluster jobs.

1. `conda activate isweep`

2. `snakemake -c1 -n -s Snakefile-roi.smk --configfile eur.yaml`

3. `nohup snakemake -c1 -s Snakefile-roi.smk --keep-going --latency-wait 300 --cluster "qsub -q some-queue.q -N {rule} -m e -M your.email.@university.edu -pe local 12 " --configfile eur.yaml --jobs 10 & `
    - `-pe local 12` allocates 12 CPUs on same node
    - `-m e -M your.email@university.edu` may send you 100s of emails

4. `head eur/*/results.hap.tsv` ; `head eur/*/results.snp.tsv` ; `head eur/*/ibd.gini.tsv` -->

## Other considerations

Review these files:
- An IBDNE file from Step 6 in `workflow/scan/README.md`
- `roi.tsv` from your selection scan
    - You may wish to copy/rename this file and make some of the following suggestions.
    - You can rename the hit* to loci names if you know some.
    - You can set values under column `ALPHA` to choose your one minus the confidence level.
    - You can use a different sweep model (column `MODEL`):
        - Default 'a' for additive
        - 'm' for multiplicative
        - 'd' for dominance
        - 'r' for recessive
- `plots/`
    - Check if there are lots of recombination hotspot
    - Check if there is low marker density

<!-- Make sure that you have an IBDNE file from .
- Using `workflow/scan/scripts/run-ibdne.sh` -->

The pipeline uses `hap-ibd.jar`. You can use different IBD segment detection methods if you wish.

Refer to the scripts in `rules/` for how to run the commands with different IBD segment methods.

<!-- We have done limited data analysis on array data. You can use the `workflow/simulate` pipeline, adjust simulation studies for your data, and assess the accuracy of methods for difference configurations. -->