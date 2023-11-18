## You run this pipeline after the workflow/scan/ pipeline.

You should review the files:
- roi.tsv from your selection scan
    - You can rename the hit* to loci names if you know some
- plots/
    - Check if there are lots of recombination hotspot
    - Check if there is low marker density
        - ibd-ends should address this
        - Read the paper about how sparsity could affect results

## The procedure

1. Edit the *.yaml file
2. conda activate isweep
3. snakemake -c1 -s Snakefile-roi.smk --configfile *.yaml
    - "-n" option is the dry run
4. nohup snakemake -c1 -s Snakefile-roi.smk --cluster "[options]" --jobs X --configfile *.yaml &
    - Make sure to pass in the isweep environment to cluster

I will write some scripts (soon-ish) to concatenate results for the loci.
