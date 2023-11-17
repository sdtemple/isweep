# isweep-simstudy-ehh
Pipeline for simulation study applying methods based on extended haplotype homozygosity
- Assumes existing files from a forward simulation w/ selection
	- Mainly, a file slimulation.trees
	- An existing folder structure
	- Made by snakemake pipeline in isweep package
- Assumes a package isweep here (https://github.com/sdtemple/isweep)
- Assumes a package selscan (https://github.com/szpiech/selscan)
- Assumes a package isafe (https://github.com/alek0991/iSAFE)

Statistics computed (see their citations, manual):
- selscan
	- iHS
	- iHH12
	- nSL
- isafe

To run:
- Create an environment with some standard python packages
	- See isafe link
	- Bring in snakemake as well
- Use pre-compiled selscan bin/linux/selscan, or compile according to their instructions

Run the pipeline on a computing cluster w/ the command (qsub-based):
- `conda activate <your-isafe-environment>`
- `nohup snakemake all -s Snakefile-isafe-selscan -c1 -jobs XXX --cluster "qsub -q <name-of-queue.q> -V" --keep-going --latency-wait XXX --rerun-triggers mtime --rerun-incomplete --configfile compare.yaml &`
