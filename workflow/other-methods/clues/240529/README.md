## My instructions to run CLUES version 2
## Seth Temple
## May 29, 2024

### Setup

You need to get their repository
`git clone https://github.com/avaughn271/CLUES2`

You need install the required packages. Follow their instructions on repo.
`conda env create -f clues-env-v2.yaml` or `mamba env create -f clues-env-v2.yaml`

You will also need to install Relate for these scripts. Keep track of where you put your Relate folder.
`https://myersgroup.github.io/relate/index.html`

You can use a different ARG reconstruction tool like SINGER. See CLUES repo.
You will need to modify the scripts in such a scenario.

Here I provide a python file `inference-Ne-indexing-fixed.py`.
This is credit to the CLUES authors, as a modification of their `inference.py` file.
I only modify how they bring in the Ne file, which I find is error-prone.
If running an analysis with constant Ne, you can use their `inference.py`.
The associated shell script is clues-run-constant-N.sh.


### Instructions

CLUES and CLUES2 and Relate are not easy to use with Snakemake.
They sometimes crash for unknown reasons, which can end the entire pipeline.

Check out CLUES or Relate if any of these crash.
Double check that input in shell scripts are right.

See arguments below. Meant to be run in order.
Submit these scripts solo to the cluster:
- coarse-coal.py
- `bcftools query -l $FOLDER/$VCFFILE | some-number > $SUBSAMPLE`
- tb.sh 
- - This runs tabix
- - `tb.sh $FOLDER $VCFFILE`
- subset-chr.sh
- - This runs bcftools
- - `subset-chr.sh $FOLDER $PREFIX $VCFFILE $LEFT $RIGHT $CHRNUM $SUBSAMPLE`
- relate-file-format.sh
- - `relate-file-format.sh $RELATEFOLDER $FOLDER $PREFIX`
- relate-run.sh
- `relate-run.sh $RELATEFOLDER $FOLDER $PREFIX $MU $DEMO $GENMAP`
- relate-branch-timeb.sh
- - This is to run CLUES1
- - `relate-branch-timeb.sh $RELATEFOLDER $FOLDER/$PREFIX $DEMO clues2 $MU $NUMMCMC $LOC`
- relate-branch-newick.sh
- - `get-daf.py slimulation.freq`
- - This is to run CLUES2
- - `relate-branch-newick.sh $RELATEFOLDER $FOLDER $PREFIX $DEMO clues2 $MU $NUMMCMC $LOC`
- `python $ALLELESCRIPT $FOLDER/$PREFIX.haps $FOLDER/derived.alleles.txt $LOC`
- relate-to-clues.sh
- - `relate-to-clues.sh output-file-from-relate-branch-newick.sh output-file-from-derived-allele-command-above $FOLDER clues2`
- clues-run.sh
- - `clues-run.sh $CLUESSCRIPT $DAFSCRIPT $FOLDER clues2 clues2_times.txt slimulation.freq  $DEMO $DOMINANCECOEF $DF $TCUTOFF`
- - slimulation.freq is text file with current allele frequency in last row, first column 
- clues-run-constant-N.sh
- - Same as clues-run.sh but replace $DEMO with the constant Ne
- - if you have constant Ne

Use ` sbatch --partition= --job-name= --mail-type= --mail-user= --mem= `.
They all run on one CPU, so no need to specify ` --cpus-per-task= ` 

Relate creates a lot of temporary files and folders. Delete these after.
CLUES2 wants to use a large Newick tree. Delete this after for memory purposes.
Supposedly CLUES2 is much faster than CLUES1, so use the `relate-branch-newich.sh` w/ `relate-to-clues.sh` approach.

You will need to specify the following parameters:
- FOLDER : where you have simulated data (vcf)
- LOC : the location of the causal allele
- LEFT : left endpoint for smaller vcf (bp)
- RIGHT : right endpoint for smaller vcf (bp)
- CHRNUM : `#CHROM` in vcf file
- VCFFILE : name of a vcf file
- PREFIX : probably best to keep this prefix consistent throughout procedure
- DF : the number of frequency bins for CLUES
- NUMMCMC : the number of samples of tree branches
- MU : the mutation rate
- GENMAP : the genetic map
- DEMO : the demography file
- - Use my script ` coarse-coal.py ` to make this file
- - `python coarse-coal.py input-Nefile.ne header-file.py output-Ne `
- - input-Nefile.ne is the file format in IBDNe.jar
- - header-file.txt is one row with `group1` and a second row w/ times in generations
- - output-Ne is what you want to name it as your DEMO parameter
- TCUTOFF : the time in generation cutoff for CLUES
- DOMINANCECOEF : h in 1:(1+hs):(1+s)
- DAFSCRIPT : my script to get the population allele frequency
- - `get-daf.py`
- ALLELESCRIPT : my script to get the derived alleles text file
- - `derived-allele-file.py`
- CLUESSCRIPT : script to run CLUES
- - I modify the small section on Ne to work with my Ne file specification
- - `inference-Ne-indexing-fixed.py`

CLUES recommends for a higher DF and NUMMCMC:
- They use DF=600 in their 2023 pre-print
- They use NUMMCMC=3000 in their 2023 pre-print
- They use TCUTOFF=536 in their 2023 pre-print

You will want to use ` for k in $(seq k1 1 k2) ; do ; done ;` to submit jobs in mass. These may take a few hours to a couple days.

