## Which are which pipelines

There are multiple snakemake pipelines in these subfolders.

- `scan/` implements the genome-wide selection scan for real data

- `roi/` implements the loci-specific analyses for real data

- `simulate/` implements simulation study on genetic data
  - `Snakefile-simulate.smk` will use SLiM to get `.trees` and `.vcf.gz`
  - `Snakefile-ibd.smk` will run scripts like in `scan/` and `roi/`

- `other-methods/` implements other methods to compare in simulation study

- `coalescent/` has some python scripts for simulation study on coalescent data 

You will likely use `scan/` and then `roi/`.

There is documentation in each subfolder.

You will edit the .yaml files to adjust settings and direct paths to where data is.