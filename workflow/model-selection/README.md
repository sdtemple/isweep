# Second pipeline in data analysis:

## The procedure

1. Run `ibdne.jar` to get recent effective sizes. Or another method estimate recent Ne.
    - Read `cat scripts/run-ibdne.sh ` or `head -n 20 scripts/run-ibdne.sh ` or 'vim scripts/run-ibdne.sh `.
    - `bash scripts/run-ibdne.sh [options] `
        - $1 : location of `ibdne.jar` applet
        - $2 : memory limit in Gb
        - $3 : folder name for the study
        - $4 : subfolder path for ibd segments
        - $5 : lowest chromosome number (e.g., 1)
        - $6 : highest chromosome number (e.g., 22)
        - $7 : prefix for outfiles 
    - You should do this with cluster computing. Send in a job!
      - E.g., `sbatch [cluster-options] scripts/run-ibdne.sh [script-options]`

2. Edit the YAML file.

3. ` mamba activate isweep `

4. Dry-run of the pipeline. Check what rules will run.
```
snakemake -c1 -s Snakefile-roi.smk --configfile *.yaml 
```
- ` -n ` option is the dry run

5. Run the pipeline for real.
```
nohup snakemake -c1 -s Snakefile-roi.smk --cluster "[options]" --jobs X --configfile *.yaml & 
```
- See `misc/cluster-options.md` to choose SLURM or QSUB options.
- Make sure to pass in the isweep environment to cluster.


Outputs:
- `summary.hap.norm.tsv` are estimated selection coefficients, and other estimates, for regions of interest.
  - Read Temple, Waples, and Browning (AJHG, 2024) to learn about the estimates.
  - Confidence intervals assume IBD rates are (asymptotically) normally distributed.
  - Frequency estimate is based on the best differentiated SNP **subset**.
  - Models are 'a' additive, 'm' multiplicative, 'd' dominance, and 'r' recessive.
- Other types of confidence intervals.
  - 'perc' wildcard means percentile-based confidence intervals.
  - 'snp' wildcard means that frequency estimate is based on best differentiated SNP.

## Other considerations

Review these files:
- An IBDNE file from Step 1 in `workflow/scan/README.md`. The Ne estimates shouldn't go haywire.
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