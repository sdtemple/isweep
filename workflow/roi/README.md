## Second pipeline in data analysis:

#### Run this pipeline after the `workflow/scan/` pipeline.

You should review the files:
- `roi.tsv` from your selection scan
    - You can rename the hit* to loci names if you know some
- `plots/`
    - Check if there are lots of recombination hotspot
    - Check if there is low marker density
        - ibd-ends.jar should address this
        - Read the paper about how sparsity could affect results

The pipeline uses hap-ibd.jar. You can use different IBD segment detection methods if you wish.

In this case, refer to the scripts in `rules/` for how to run the commands with different IBD segment methods.

We have done limited data analysis on array data. You can use the `workflow/simulate` pipeline, adjust simulation studies for your data, and assess the accuracy of methods for difference configurations.

## The procedure

1. Edit the *.yaml file
2. ` conda activate isweep `
3. ` snakemake -c1 -s Snakefile-roi.smk --configfile *.yaml `
    - ` -n ` option is the dry run
4. ` nohup snakemake -c1 -s Snakefile-roi.smk --cluster "[options] -V " --jobs X --configfile *.yaml & `
    - See `misc/cluster-options.md` to choose SLURM or QSUB options.
    - Make sure to pass in the isweep environment to cluster

## Example

I use the qsub system to manage cluster jobs.

1. `conda activate isweep`

2. `snakemake -c1 -n -s Snakefile-roi.smk --configfile eur.yaml`

3. `nohup snakemake -c1 -s Snakefile-roi.smk --keep-going --latency-wait 300 --cluster "qsub -q some-queue.q -N {rule} -m e -M your.email.@university.edu -pe local 12 " --configfile eur.yaml --jobs 10 & `
    - `-pe local 12` allocates 12 CPUs on same node
    - `-m e -M your.email@university.edu` may send you 100s of emails ---> use an email filter

4. `head eur/*/results.hap.tsv` ; `head eur/*/results.snp.tsv` ; `head eur/*/ibd.entropy.tsv`

