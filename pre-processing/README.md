# Pre-processing data

The whole point is to get a subsample from a larger dataset.

<!-- Input data
  - Phased data
  - No apparent population structure
  - No apparent close relatedness -->

I put some pre-processing scripts in this folder, like:
- running a global ancestry inference;
- running a prinicipal components analysis;
- running `IBDkin`;
- and, haplotype phasing.

Now, in an automated fashion, you could get:
- `https://github.com/sdtemple/flare-pipeline`
- Haplotype phase
- Local & global ancestry inference
- Long IBD segments

## What some scripts do

- `high-kinship.py` filters the ibdkin output to only use 
- `keep-one-family-member.py`
  - This forms connected components that have inferred relation < 3 degrees.
  - Then, keep only 1 individual from each connected component.
  - These connected components are like extended families
- The "meta" python scripts write R scripts for plotting w/ `ggplot2`.
  - The `sns.*.py` probably does the same but w/ `matplotlib`
  - These assume you have a RACE and STUDY column.
- The `run-*.sh` bash scripts call your R or python scripts for cluster computing.
  - `run-ibdkin.sh` assumes you have some existing file about the chromosme beginning and ends.
- In `hapibd.o` you can see how I ran `hap-ibd.jar` and how long it ran for 37k samples.
<!-- - I used some commands like
  - `bcftools query -l ... > file.txt`
  - `bcftools query -f "%CHROM\t%POS" ... > file.txt` -->
- I followed the instructions in ADMIXTURE for global ancestry inference.
  - http://dalexander.github.io/admixture/

See examples to run `beagle.jar`, `hap-ibd.jar` from Brian Browning.

## Citation

You need to cite these software!

## Caution:

I don't remember when I last ran or edited these scripts.

Use them as inspiration for pre-processing your data.
  <!-- - The point of these scripts is to address population structure and relatedness.
  - Use the linux command `head -n 10` to look at inputs. -->
<!-- - For the R scripts you need `GENESIS` and `SNPRelate`.
  - See our key resources table in paper.
    - https://bioconductor.org/packages/release/bioc/html/GENESIS.html
    - https://bioconductor.org/packages/release/bioc/html/SNPRelate.html -->
<!-- - You shouldn't do PCA or IBD calling or global ancestry w/ large samples on your personal laptop
  - Use a computing cluster -->