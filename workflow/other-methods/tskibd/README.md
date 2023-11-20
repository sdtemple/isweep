# Results with true ibd segments

This pipeline is very preliminary. I explored getting true ibd tracts from a tree sequence.

This is mostly for people developing methods to call ibd segments.

It does not look like the software gets all of them based on our methods comparison.

This may be due to the gene conversion features we added to our simulation study.

Pipeline for simulation study based on true IBD from `tskibd` package
- Assumes existing files from a forward simulation w/ selection
	- Mainly, a file slimulation.trees
	- An existing folder structure
- Assumes a package tskibd (https://github.com/bguo068/tskibd)

To run:
- Create an environment with some standard python packages
	- numpy
	- pandas
	- scipy
	- etc.
- Compile tskibd according to their instructions

Run the pipeline on a computing cluster w/ the command (qsub-based):

`nohup snakemake all -s Snakefile-tskibd.smk -c1 -jobs XXX --cluster "qsub -q name-of-queue.q -V" --keep-going --latency-wait XXX --rerun-triggers mtime --rerun-incomplete --configfile trueibd.yaml &`

