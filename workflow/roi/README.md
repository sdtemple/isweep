## This is the second pipeline in a real data analysis:

#### You run this pipeline after the workflow/scan/ pipeline.

You should review the files:
- `roi.tsv` from your selection scan
    - You can rename the hit* to loci names if you know some
- `plots/`
    - Check if there are lots of recombination hotspot
    - Check if there is low marker density
        - ibd-ends.jar should address this
        - Read the paper about how sparsity could affect results

## The procedure

1. Edit the *.yaml file
2. ` conda activate isweep `
3. ` snakemake -c1 -s Snakefile-roi.smk --configfile *.yaml `
    - ` -n ` option is the dry run
4. ` nohup snakemake -c1 -s Snakefile-roi.smk --cluster "[options] -V " --jobs X --configfile *.yaml & `
    - Make sure to pass in the isweep environment to cluster

I will write some scripts (soon-ish) to concatenate results for the loci.

## Example

I use the qsub system to manage cluster jobs.

1. `conda activate isweep`

2. `snakemake -c1 -n -s Snakefile-roi.smk --configfile eur.yaml`

3. `nohup snakemake -c1 -s Snakefile-roi.smk --keep-going --latency-wait 300 --cluster "qsub -q some-queue.q -N {rule} -m e -M your.email.@university.edu -pe local 12 " --configfile eur.yaml --jobs 10 & `
    - `-pe local 12` allocates 12 CPUs on same node
    - `-m e -M your.email@university.edu` may send you 100s of emails ---> use an email filter

4. `head eur/*/results.hap.tsv` ; `head eur/*/results.snp.tsv` ; `head eur/*/ibd.entropy.tsv`

