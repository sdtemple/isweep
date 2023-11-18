## This is the first pipeline you run for real data analysis

You must phase your data beforehand (Beagle, Shapeit, etc.)

You must remove close relatives beforehand (IBDkin, plink, etc.)

You must subset to different "ancestry" groups (ADMIXTURE, principal components, etc.)

If you are studying a subset of a larger (biobank) dataset, create a subsample text file.
- Link to it in your yaml file
- If not, create an empty subsample.txt file
    - "touch subsample.txt"

## The procedure

1. Edit the *.yaml file
2. conda activate isweep
3. snakemake -c1 -s Snakefile-scan.smk --configfile *.yaml
    - "-n" option is the dry run
4. nohup snakemake -c1 -s Snakefile-scan.smk --cluster "[options]" --jobs X --configfile *.yaml &
    - Make sure to pass in the isweep environment to cluster
    - Search for other snakemake options if you so wish
5. Make a manhattan plot?
    - put a path to my plot script
    - python manhattan.py
6. Run ibdne.jar to get recent effective sizes
    - `cat run-ibdne.sh ` or `head -n 20 run-ibdne.sh``
    - `bash run-ibdne.sh [options] `

Steps 5. and 6. are if you want to estimate selection coefficients for loci.