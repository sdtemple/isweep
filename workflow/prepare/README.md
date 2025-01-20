# Pipeline for local ancestry & IBD inference 

Perform phasing and local ancestry inference for an analysis of target admixed samples. 

Also, call IBD segments.

<img src="dag-240418.png" align="center" width="600px"/>

<!-- You can look at `dag-*` files to see what the pipeline looks like / looked like. -->

## Installation

1. Clone the repository.
```
git clone https://github.com/sdtemple/flare-pipeline
```
2. Install `java` so as to be able to run `java -jar` on terminal.
3. Download additional software (requires `wget`).
```
bash get-software.sh software
```
4. Make the Python enviroment.
```
mamba env create -f conda-env.yml
```
- I recommend using `mamba` of some sort (https://mamba.readthedocs.io/en/latest/index.html)
- In which case, `mamba env create -f conda-env.yml`
- See `misc/installing-mamba.md` for further support

## Requirements

- GDS files for each chromosome
    - Reference samples
    - Target admixed samples
- Map between reference sample IDs to their reference panel
    - Point to this file in your YAML settings

To subset or manipulate VCFs and GDSs:
- `scripts/vcf-to-gds.R`
- `scripts/get-sample-ids-gds.R`
- `scripts/subset-gds.R`: implements subsetting, and minor allele count and missingness filters
    - Use MAC 1 and missingness 0.00 if you only want to subset samples

## Run the pipeline 

1. `mamba activate flare24`
2. Modify the `your.analysis.arguments.yaml` file
    - See the `change:` settings
    - You need to choose a reference sample!
3. Do a dry-run of the pipeline. Check which rules will run.
```
snakemake -c1 -n --configfile *.yaml
```
4. Run the pipeline for real.
```
nohup snakemake -c1 --latency-wait 300 --keep-going --cluster "[options]" --configfile *.yaml --jobs XXX &
```
- Other useful `snakemake` commands
    - `--rerun-incomplete`
    - `--rerun-triggers mtime`
    - `--force-all` ???
- Commands for `qsub`
    - `--cluster "qsub -q your-queue.q -m e -M your.email@uni.edu -pe local XXX -l h_vmem=XXXG -V" `
    - "-V" is important to pass in your conda environment!
    - "-pe local XXX" is how many threads you will use
    - You don't have to send emails to yourself if you don't want to
- Commands for `slurm`
    - `--cluster "sbatch [options]" `
        - "-e ~/your-logs/{rule}.e" and "-o ~/your-logs/{rule}.o" will control where stderr, stdout go
        - "--cpus-per-task=XX" says how many cpus per job
        - "--nodes=XX" says how many nodes per job
        - "--partition=SOMENAME" says which partitions to use
        - "--mem=XX" says how much memory in MB
        - "--mail-type=ALL" and "--mail-user=your.email@university.edu" sends mail to you when a job finishes
        - "--job-name={rule}"
    - You can sign out of cluster. `no hup ... &` will keep this as an ongoing process until complete
5. Your LAI results in a `lai/` folder
6. Your IBD results in a `ibdsegs/` folder
    - By default, do not detect IBD segments

For reproducibility, the `arguments.yaml` in the main folder says what you ran. Don't change it ever!

For robustness, you can create different `*.yaml` settings and see how results change. 

Make sure to change the `your-analysis-folder` setting. 

## Other notes

### Possible bugs and errors

- Running out of memory: increase cluster-resouces:xmxmem in yaml and terminal pass into `--cluster`
- Names are different in CHROM column between reference and target samples
- JAR file is corrupted: download a fresh version
    - Changed the `remove-phase.jar` to `remove-phase.py`

### Phasing strategies

There are two phasing strategies:
- Use the reference to phase the target sample (could introduce imputed values)
    - This may result in more markers in the target sample data
    - You can also use imputed markers if you change the flag in Beagle
- Rephase target and reference targets altogether

### A pilot study on small chromosomes

Your chromosome files should be named numerically, e.g., chr1 all the way to chr22. If you want to study some sex chromosome or otherwise, use symbolic links to say call it chr23.

In your YAML configuration file, the chromosome files between `chr-low` to `chr-high` will be analyzed. Use a subset of the smallest chromosomes to test the pipeline, e.g., `chr-low: "21"` and `chr-high: "22"`. Then, to run the entire analysis, use `chr-low: "1"` and `chr-high: "22"`.

<!-- ### Initiating with *.vcf.gz files instead of *.gds files

The pipeline starts with GDS files for each chromosome. It converts these to VCF files.

If you have VCF files already, don't waste your time converting VCF to GDS.

Use `ln -s initial/folder/location/chr*.vcf.gz your/analysis/folder/gtdata/adxpop/chr*.vcf.gz` to create a symbolic link.
- Do this in `your-analysis-folder/gtdata/adxpop` and `your-analysis-folder/gtdata/refpop/`
    - You may need to
        1. `mkdir your-analysis-folder`
        2. `mkdir your-analysis-folder/gtdata`
        3. `mkdir your-analysis-folder/gtdata/adxpop`
        4. `mkdir your-analysis-folder/gtdata/refpop`
- If you have VCFs for references and GDS for admixed, do the symbolic link for the references only.
    - And vice versa

The minor allele count and minor allele frequency filters are applied redundantly:
- In the GDS to VCF conversion
- And shrinking the VCF if you start from VCF -->

### keep-samples and exclude-samples arguments

The `keep-samples` in the rule `gds_to_vcf_adx` controls memory and time if you have a huge admixed target sample, but you only want to study a subset of them.
- If you want to study all of them, make a one column text file with sample IDs for all samples

The `exclude-samples` is an option in the Browning Lab software (beagle, flare, etc). This is a one column text file with sample IDs for the samples you do not want to study. Most of the time this should be the empty file `excludesamples.txt` in this repo. If you use it with a non-empty file, it will subset the VCFs in the phasing step.

### rename-chrs-map argument

Sometimes the recombination map and the VCF/GDS data have different names for the chromosome column, e.g., chr22 versus 22. Choose between one of four text files in the `rename-chrs/` from this repo.
- The first column is the old name, the name of chromosome in the VCF file.
- The second column is the new name, the name of the chromosome in the map file.
- There are separate YAML options for the reference and admixed target samples. 

### Calling IBD segments

You can call detect IBD segments by removing the comments in `record_yaml` rule.
- This may be useful if you want to do:
    - IBD mapping in conjunction w/ local ancestry analyses
    - Relatedness inference using `IBDkin` (look for github repo)

### Limited cluster resources

- See the Beagle paper about manipulating window size
    - Impacts memory and runtime
    - The greatest concern is memory in large samples
    - https://www.cell.com/ajhg/fulltext/S0002-9297(21)00304-9
- Subset the target samples and implement reference phasing
    1. Make subsets of the target samples
    2. Make separate configuration YAML files w/ different `keep-samples` text files
    3. Comment out this: `[macro+'/lai/chr'+str(i)+'.rephased.flare.anc.vcf.gz' for i in range(low,high+1)],`
    4. Keep this: `[macro+'/lai/chr'+str(i)+'.referencephased.flare.anc.vcf.gz' for i in range(low,high+1)],`

## Development things to do

- This repo currently uses snakemake 7.25.2
    - May extend to version 8 as I develop familiarity in other repos
<!-- - Implement other local ancestry inference software
    - For example, MOSAIC from Salter-Townshend and Myers
- Impute other phasing software
    - For example, SHAPEIT -->
- Initial data can be `*.gds` OR `*.vcf`

## Citation

It takes time to make these pipelines. I am a early career researcher and appreciate credit. 

(smiley face) Please:
- At time of publication, check this repo to see if we have a paper introducing the pipeline.
    - (If we publish this pipeline somewhere, I will point out the paper.)
- Acknowledge me in publication.

You should be citing (see `get-software.sh`):
- `beagle`
- `flare`

<!-- ## Contact

Seth D. Temple

sdtemple.github.io

sdtemple@uw.edu -->