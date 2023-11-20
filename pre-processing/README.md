## Pre-processing data

The whole point is to get a subsample from a larger dataset.

Input data
  - Phased data
    - Use Beagle (https://faculty.washington.edu/browning/beagle/beagle.html)
  - No apparent population structure (i.e., admixture, diverged groups in same dataset)
  - No apparent close relatedness
    - Use IBDkin (https://github.com/YingZhou001/IBDkin)
    - Or something else like plink (https://zzz.bwh.harvard.edu/plink/)

I will write put some pre-processing scripts in this folder, like
- Running a global ancestry inference
- Running a prinicipal components analysis
- Running IBDkin
- Phasing

## Caution:
- I don't remember when I last ran or edited these scripts.
- Use them as inspiration for pre-processing your data.
  - The point of these scripts is to address population structure and relatedness.
  - Use the linux command `head -n 10` to look at inputs
- For the Rscripts you probably need GENESIS and SNPRelate
  - See our key resources table in paper
  - https://bioconductor.org/packages/release/bioc/html/GENESIS.html
  - https://bioconductor.org/packages/release/bioc/html/SNPRelate.html
- You shouldn't do PCA or IBD calling or global ancestry w/ large samples on your personal laptop
  - Use a computing cluster

## Probably what some scripts do

- high-kinship.py filters the ibdkin output to only use 
- keep-one-family-member.py
  - I think this will form connected components that have inferred relation < 3 degrees
  - Then keep only 1 individual from each connected component
  - These connected components are like extended families
    - So if you have a family reunion with all your first and second cousins
- The "meta" python scripts probably writes R scripts for plotting w/ ggplot2
  - The sns.*.py probably does the same but w/ matplotlib
  - These probably assume you have a RACE and STUDY column
- The "run-*.sh" bash scripts probably calls your R or python scripts for cluster computing
  - run-ibdkin.sh probably assumes you have some existing file about the chromosme beginning and ends
- In hapibd.o you can see how I ran hap-ibd.jar and how long it ran for 37k samples
  - Ran this for ibdkin analysis
- I probably used some commands like
  - `bcftools query -l ... > file.txt`
  - `bcftools query -f "%CHROM\t%POS" ... > file.txt`
- I probably followed the instructions in admixture for global ancestry inference
  - http://dalexander.github.io/admixture/

See examples to run beagle.jar, hap-ibd.jar from Brian Browning

## Citation

You need to cite these software!