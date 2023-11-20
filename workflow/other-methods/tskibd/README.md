# isweep-simstudy-tskibd
Pipeline for simulation study based on true IBD from tskibd package
- Assumes existing files from a forward simulation w/ selection
	- Mainly, a file slimulation.trees
	- An existing folder structure
- Assumes a package isweep here (https://github.com/sdtemple/isweep)
- Assumes a package tskibd (https://github.com/bguo068/tskibd)

To run:
- Create an environment with some standard python packages
- Compile tskibd according to their instructions

Run the pipeline on a computing cluster w/ the command (qsub-based):

`nohup snakemake all -s Snakefile-tskibd -c1 -jobs XXX --cluster "qsub -q name-of-queue.q -V" --keep-going --latency-wait XXX --rerun-triggers mtime --rerun-incomplete --configfile trueibd.yaml &`

