## This is the first pipeline you run for real data analysis

You must phase your data beforehand (Beagle, Shapeit, etc.)
- For Beagle, something like ` java -jar beagle.jar gt=chr.vcf.gz map=chr.map out=chr.phased `
-- Send this as a cluster job !

You must remove close relatives beforehand (IBDkin, plink, etc.)
- Recommendation: at least remove relations w/ kinship > 0.125

You must subset to different "ancestry" groups (ADMIXTURE, principal components, etc.)
- Recommendation: use subsets w/ at least 85% the same global ancestry inference

If you are studying a subset of a larger (biobank) dataset, create a subsample text file.
- Link to it in your yaml file
- If not, create an empty subsample.txt file
    - ` touch subsample.txt `

## The procedure

1. Edit the *.yaml file
2. ` conda activate isweep `
3. `bash maps.sh maps-folder map-prefix map-suffix study-name chrlow chrhigh`
    - $1 : folder where genetic maps are
    - $2 : the prefix before the chromosome number
    - $3 : the suffix after the chromosome number
    - $4 : define your study folder name (match to .yaml file)
    - $5 : lowest chromosome number (e.g., 1)
    - $6 : highest chromosome number (e.g., 22)

3. ` snakemake -c1 -s Snakefile-scan.smk --configfile *.yaml `
    - "-n" option is the dry run

4. `nohup snakemake -c1 -s Snakefile-scan.smk --cluster "[options]" --jobs X --configfile *.yaml & `
    - Make sure to pass in the isweep environment to cluster
    - Search for other snakemake options if you so wish

5. Make a manhattan plot?
    - put a path to my plot script
    - ` python terminalscripts/manhattan.py [options]`
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
    - `cat terminalscripts/run-ibdne.sh ` or `head -n 20 terminalscripts/run-ibdne.sh `
    - `bash terminalscripts/run-ibdne.sh [options] `
        - $1 : location of ibdne.jar applet
        - $2 : memory limit in Gb
        - $3 : folder name for the study
        - $4 : subfolder path for ibd segments
        - $5 : lowest chromosome number (e.g., 1)
        - $6 : highest chromosome number (e.g., 22)
        - $7 : prefix for outfiles 

Steps 5. and 6. are if you want to estimate selection coefficients for loci.

## Example

I use the qsub system to manage cluster jobs.

1. `conda activate isweep`
2. ` bash maps.sh path-to-maps-folder/ decode.chr .map eur 1 22 `
3. `snakemake -c1 -n -s Snakefile-scan.smk --configfile eur.yaml`

4. `nohup snakemake -c1 -s Snakefile-scan.smk --keep-going --latency-wait 300 --cluster "qsub -q some-queue.q -N {rule} -m e -M your.email.@university.edu -pe local 12 " --configfile eur.yaml --jobs 22 & `
    - `-pe local 12` allocates 12 CPUs on same node
     - `-m e -M your.email@university.edu` may send you 100s of emails ---> use an email filter

5. `python terminalscripts/manhattan.py scan.ibd.tsv eur.manhattan 1 22 4. 3. "TOPMed European Americans" 12 3` 


6. ` echo "bash terminalscripts/run-ibdne.sh ibdne.jar 100 eur ibdsegs/ibdends/modified/scan 1 22 ibdne " | qsub -q some-queue.q -N ibdne -m e -M your.email@university.edu `

You don't have to use my `manhattan.R` or `manhattan.py` files for plotting. Feel free to take inspiration from the scripts.