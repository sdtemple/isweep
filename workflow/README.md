## Which are which pipelines

There are multiple snakemake pipelines in these subfolders.

- `scan-selection/` implements the genome-wide selection scan for real data.

- `model-selection/` implements the loci-specific sweep analyses for real data.

- `simulate/` implements selection simulation study on genetic data.
  - `Snakefile-simulate.smk` will use SLiM to get `.trees` and `.vcf.gz`
  - `Snakefile-ibd.smk` will run scripts like in `scan-selection/` and `model-selection/`

- `scan-case-control` implements a genome-wide scan for differences in case and control IBD rates.

There is documentation in each subfolder.

You will edit the `*.yaml` files to adjust settings and direct paths to where data is.
