# Optional third pipeline in data analysis:

If you have binary/categorical phenotypes, you can evaluate difference in IBD rates for these groups.

The common application is IBD mapping where you have binary case and control phenotypes.\

You run this pipeline after running the selection scan pipeline.
1. Strong selection can confound the test.
- You want to make sure "significant" regions here are not confounded by "significant" regions in the selection scan.
2. The selection scan automatically will detect IBD segments for you and organize the folder.

Methodologically, this procedure is very similar to the selection scan. We just apply it to the difference in IBD rates.
- Under the same Temple and Thompson asymptotic conditions, this test should be reasonable.
- p values come from the Gaussian distribution, when the Temple and Thompson conditions and the null model approximately hold.

## The procedure:
---

#### Scanning for anomalous differences in IBD rates
---

1. Edit the YAML file.
- The `case.yaml` file gives default settings.

2. ` mamba activate isweep `

3. Dry-run of the scanning pipeline.
```
snakemake -c1 -s Snakefile-case.smk --configfile *.yaml
```
- "-n" option is the dry run

4. Run the scanning pipeline for real.
```
nohup snakemake -c1 -s Snakefile-scan.smk --cluster "[options]" --jobs X --configfile *.yaml &
```
- See `misc/cluster-options.md` to choose SLURM or QSUB options.
- Make sure to pass in the `isweep` environment to cluster.
- Search for other `snakemake` options if you so wish.

5. Make a customized IBD rates scan plot.
    - `python scripts/plotting/plot-scan-case-pipeline.py [options]`

Outputs:
- `scan.case.ibd.tsv` should have all the data for the scanning statistics and thresholds.
  - 'Z1' and 'Z0' variables are standardized/normalized IBD rates for cases and controls.
  - 'ZDIFF' is the difference of 'Z1' and 'Z0'.
  - 'ZDIFFZ' is the standardized/normalized version of 'ZDIFF'.
  - The 'RAW' thresholds apply to 'ZDIFF'. The 'Z' thresholds correspond to 'ZDIFFZ'.
  - p values assume that IBD rates are (asymptotically) normally distributed.
- `roi.case.tsv` are your significant regions.
- `autocovariance*.png` is autocovariance by cM distance. The black line is a fitted exponential curve.
- `zhistogram*.png` is a default histogram for the IBD rates standardized. It should "look Gaussian".
    - There are ones for cases, controls, and their differences in rates.
- `scan.case.png` is a default plot for the IBD difference scan.
- `fwer.analytical.case.txt` gives parameters and estimates for multiple-testing difference scan.

#### Evaluating phenotypes of clusters with excess IBD sharing
---

1. Edit the YAML file.
- The `case.roi.yaml` file gives default settings.

2. ` mamba activate isweep `

3. Dry-run of the putlier pipeline.
```
snakemake -c1 -s Snakefile-outlier.smk --configfile *.yaml
```
- "-n" option is the dry run

4. Run the outlier pipeline for real.
```
nohup snakemake -c1 -s Snakefile-outlier.smk --cluster "[options]" --jobs X --configfile *.yaml &
```

Outputs:
- Look into the subfolders for `outlier*.phenotype.txt` files.
  - First column is the haplotype ID, with '_1' and '_2' being the haplotype and the string before the individual.
  - Second column is the case or control status, as defined in your phenotypes file.
- Also, look into the `matrix.outlier.phenotypes.tsv` files.
  - These have categorical phenotypes (binary for each outlier group) in a linear model X matrix.
  - And a binary column 'OUTLIER_ANY' if the haplotype is in an outlier group period.


## Broad use of IBD outlier group detection
--

**Technically, this `Snakefile-outlier.smk` workflow should work for any numeric phenotype.**

The important script to make the design matrix (X in linear regression) is `scripts/utilities/make-design-matrix.py`. Example output is here in `design.sorted.tsv` and `design.tsv`.

### Phenotype file
---

This is a tab-separated file with two columns. First column is a sample ID. Second column is numeric. An example is `phenotypes.txt`

### Regions of interest file
---
This is a tab-separated file with a header and the following columns (at least):
- NAME (example: 'LCT')
- CHROM (example: 2)
- BPCENTER (example: 136000000) (IBD segments detected overlapping this)
- BPLEFTCENTER (example: 134000000) (defines left end of region)
- BPRIGHTCENTER (example: 138000000) (defines right end of region)