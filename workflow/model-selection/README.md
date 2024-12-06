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

The pipeline uses `hap-ibd.jar`. You can use different IBD segment detection methods if you wish.

Refer to the scripts in `rules/` for how to run the commands with different IBD segment methods.