# First pipeline in data analysis:

## Pre-processing data

You must phase your data beforehand.
- For Beagle, something like ` java -jar beagle.jar gt=chr.vcf.gz map=chr.map out=chr.phased `
    - Send this as a cluster job !
- Or, you can get haplotype phase using our pipeline `https://github.com/sdtemple/flare-pipeline`.

You should remove close relatives beforehand.
- You can use IBDkin w/ IBD >= 2.0 cM from our pipeline `https://github.com/sdtemple/flare-pipeline`.
- Recommendation: at least remove relations w/ kinship > 0.125.

You should subset to different "ancestry" groups
- You can get global ancestry estimates from our pipeline `https://github.com/sdtemple/flare-pipeline`.
- Recommendation: use subsets w/ at least 85% the same global ancestry inference.

You should do an initial run to set the inferred error rate in your `*.yaml`.
- Study one of your smallest chromosomes.
- Use `hap-ibd.jar` and `ibd-ends.jar`
- Look at the `err` result in the `*.log` file from `ibd-ends.jar`.

## The procedure:

1. Edit the YAML file.
- The `array.yaml` file gives default settings for SNP array data.
- The `sequence.yaml` file gives default settings for whole genome sequence data.

2. ` mamba activate isweep `

3. Dry-run of the pipeline.
```
snakemake -c1 -s Snakefile-scan.smk --configfile *.yaml
```
- "-n" option is the dry run

4. Run the pipeline for real.
```
nohup snakemake -c1 -s Snakefile-scan.smk --cluster "[options]" --jobs X --configfile *.yaml &
```
- See `misc/cluster-options.md` to choose SLURM or QSUB options.
- Make sure to pass in the `isweep` environment to cluster.
- Search for other `snakemake` options if you so wish.

5. Make a customized IBD rates scan plot.
    - `python scripts/plotting/plot-scan.py [options]`

Outputs:
- `scan.modified.ibd.tsv` should have all the data for the scanning statistics and thresholds.
  - 'Z' variables are standardized/normalized.
  - 'RAW' are counts.
  - p values assume that IBD rates are (asymptotically) normally distributed.
- `roi.tsv` are your significant regions.
- `autocovariance.png` is autocovariance by cM distance. The black line is a fitted exponential curve.
- `zhistogram.png` is a default histogram for the IBD rates standardized. It should "look Gaussian".
- `scan.png` is a default plot for the selection scan.
- `fwer.analytical.txt` gives parameters and estimates for multiple-testing selection scan.

## Other considerations

#### CHROM column in VCF

If the CHROM column in your VCF files doesn not align with the chromosome column in your genetic maps:
- Adjust the chromosome column in your genetic maps to match as such.
- For example, 'chr1' in VCF file and '1' in genetic map. Change to 'chr1' in genetic map.

#### Cohort-specific analyses 

If you are studying a subset of a larger (biobank) dataset, create a subsample text file.
- Link to it in your YAML file.
If you are using the entire dataset,
- ` bcftools query -l your-smallest-chr.vcf.gz > subsample.txt `
- Link this in your YAML file.

#### Ignoring some chromosomes

To not analyze some chromosomes:
- If they are in a contiguous ordering, adjust the CHRLOW and CHRHIGH parameters in YAML file
  - For example, if I want to analyze 1,2,3 or chromosomes 1,2,3,4,5,6, I set CHRLOW=1 and CHRHIGH=3.
- If they are not in a contiguous ordering, make a column of which chromosomes to exclude.
  - For example, if I want to analyze 1,2,4,6, I set CHRLOW=1, CHRHIGH=6, and make a text file with 3 and 5 in a column.
  - The path of the text file is then the parameter CHREXCLUDE in the YAML file.

#### Ignoring small chromosomes

Chromosome sizes (in cM) smaller than the FIXED:ISWEEP:CHRSIZECUTOFF will not be analyzed.

We recommend against modifying this parameter.

#### Sequence data

In our simulation studies, we have evaluated our methods for sequence data. 

The recommended sequence parameters in the YAML file are:

```
  CANDHAPIBD: # hap-ibd.ar candidate segments for input to ibd-ends

    MINSEED: '0.5' # min-seed

    MINEXT: '0.2' # min-extend

    MINOUT: '1.0' # min-output

    MINMAF: '0.40' # minimum minor allele frequency

  IBDENDS: # refined ibd segments; ibd-ends.jar

    QUANTILES: '0.5' # quantiles from endpoint posterior

    NSAMPLES: '0' # num samples from endpoint posterior

    MINMAF: '0.001'

  ISWEEP: # scanning parameters

    SCANCUTOFF: '2.0' # ibd cM threshold
```

#### Array data

We have done limited real data analysis (UKBB) on array data. 

You could use the `workflow/simulate` pipeline, adjust settings for your data, and assess accuracy for different configurations. 

The recommended array settings in the YAML file are:

```
  CANDHAPIBD: # hap-ibd.ar candidate segments for input to ibd-ends

    MINSEED: '1.8' # min-seed
    
    MINEXT: '0.5' # min-extend
    
    MINOUT: '1.8' # min-output

    MINMAF: '0.001' # minimum minor allele frequency

  IBDENDS: # refined ibd segments; ibd-ends.har

    QUANTILES: '0.5' # quantiles from endpoint posterior

    NSAMPLES: '0' # num samples from endpoint posterior

    MINMAF: '0.001'

  ISWEEP: # scanning parameters

    SCANCUTOFF: '2.0' # ibd cM threshold
```
