## What to do

These are broadly not useful to applied work.

These are for reproducibility in our papers, and me being able to go back and see what I did.

We set up snakemake pipelines or a documented sequence of shell scripts for comparing methods.

You will need to download the software from these packages.

- The clues folder implements `https://github.com/standard-aaron/clues` or `https://github.com/avaughn271/CLUES2` with trees from `https://myersgroup.github.io/relate/`
    - 240529 folder refers to May 29, 2024 work
    - 240220 folder refers to February 20, 2024 work
- The ehh folder implements iSAFE `https://github.com/alek0991/iSAFE` and EHH methods from selscan `https://github.com/szpiech/selscan`
    - EHH methods are not normalized by neutral simulations.
    - I doubt this matters when ranking variants at small selected loci.
- The imagene folder gives some scripts, an a snakemake smk file for fitting and predicting a neural network
    - `https://github.com/mfumagalli/ImaGene`
    - I provide two different python scripts to fit Regression or Categorical neural networks.
- The tskibd folder is an early exploratory analysis using true ibd segments https://github.com/bguo068/tskibd
    - We think this method may have an issue with simulated gene conversion tracts
    - It also may not get ALL the long ibd segments if window size is not small enough.