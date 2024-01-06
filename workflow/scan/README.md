## This is the first pipeline you run for real data analysis:

You must phase your data beforehand (Beagle, Shapeit, etc.)
- For Beagle, something like ` java -jar beagle.jar gt=chr.vcf.gz map=chr.map out=chr.phased `
    - Send this as a cluster job !

You should remove close relatives beforehand (IBDkin, plink, etc.)
- Recommendation: at least remove relations w/ kinship > 0.125

You should subset to different "ancestry" groups (ADMIXTURE, principal components, etc.)
- Recommendation: use subsets w/ at least 85% the same global ancestry inference

If you are studying a subset of a larger (biobank) dataset, create a subsample text file.
- Link to it in your yaml file
If you are using the entire dataset,
- ` bcftools query -l your-smallest-chr.vcf.gz > subsample.txt `
- Link this in your yaml file

You should do an initial run to set the inferred error rate in your *.yaml.
- Study one of your smallest chromosomes.
- Use hap-ibd.jar and ibd-ends.jar
- Look at the `err` result in the log file from ibd-ends.jar.

## The procedure:

1. Edit the *.yaml file
2. ` conda activate isweep `

3. ` snakemake -c1 -s Snakefile-scan.smk --configfile *.yaml `
    - "-n" option is the dry run

4. `nohup snakemake -c1 -s Snakefile-scan.smk --cluster "[options]" --jobs X --configfile *.yaml & `
    - Make sure to pass in the isweep environment to cluster
    - Search for other snakemake options if you so wish

5. Make a manhattan plot?
    - put a path to my plot script
    - ` python scripts/manhattan.py [options]`
        - $1 : tsv file name for ibd at each window
        - $2 : prefix for output image files
        - $3 : lowest chromosome number (e.g, 1)
        - $4 : highest chromosome number (e.g., 22)
        - $5 : excess IBD rate threshold (e.g., 4)
        - $6 : telomere, odd regions threshold (e.g., 3)
            - This futzes w/ regions where ibd calling is difficult (anomalously low)
        - $7 : title of plot
        - $8 : width of plot
        - $9 : height of plot

6. Run ibdne.jar to get recent effective sizes
    - `cat scripts/run-ibdne.sh ` or `head -n 20 scripts/run-ibdne.sh `
    - `bash scripts/run-ibdne.sh [options] `
        - $1 : location of ibdne.jar applet
        - $2 : memory limit in Gb
        - $3 : folder name for the study
        - $4 : subfolder path for ibd segments
        - $5 : lowest chromosome number (e.g., 1)
        - $6 : highest chromosome number (e.g., 22)
        - $7 : prefix for outfiles 
    - You should do this with cluster computing. Send in a job!

Step 5 is if you want to make a GWAS plot.
Steps 6 is if you want to estimate selection coefficients for loci.

## Example

I use the qsub system to manage cluster jobs.

1. `conda activate isweep`
2. `snakemake -c1 -n -s Snakefile-scan.smk --configfile eur.yaml`

3. `nohup snakemake -c1 -s Snakefile-scan.smk --keep-going --latency-wait 300 --cluster "qsub -q some-queue.q -N {rule} -m e -M your.email.@university.edu -pe local 12 " --configfile eur.yaml --jobs 22 & `
    - `-pe local 12` allocates 12 CPUs on same node
    - `-m e -M your.email@university.edu` may send you 100s of emails ---> use an email filter

4. `python scripts/manhattan.py scan.ibd.tsv eur.manhattan 1 22 4. 3. "TOPMed European Americans" 12 3` 

5. ` echo "bash scripts/run-ibdne.sh ibdne.jar 100 eur ibdsegs/ibdends/modified/scan 1 22 ibdne " | qsub -q some-queue.q -N ibdne -m e -M your.email@university.edu `

You don't have to use my `manhattan.R` or `manhattan.py` files for plotting. Feel free to take inspiration from the scripts.