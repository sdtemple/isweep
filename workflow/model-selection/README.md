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

You can do an analysis with percentile-based confidence intervals.
- Go into `Snakefile-roi.smk` for rule `all` and uncomment list with percentile output
- Making percentile-based intervals takes much longer
- Modify the `NBOOTPERC` parameter as desired.