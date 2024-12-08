# Pre-processing data

The whole point is to get a subsample from a larger dataset. The subsample should have no apparent population structure or relatedness. I put some scripts here to address close familial relatedness.

As well, in an automated fashion, you could get:
- `https://github.com/sdtemple/flare-pipeline`
- Haplotype phase
- Local & global ancestry inference
- Long IBD segments

## What scripts do

- `high-kinship.py` filters the ibdkin output to only use 
- `keep-one-family-member.py`
  - This forms connected components that have inferred relation < 3 degrees.
  - Then, keep only 1 individual from each connected component.
  - These connected components are like extended families
- `run-ibdkin.sh` assumes you have some existing file about the chromosome beginning and ends. It runs the [IBDkin](https://github.com/YingZhou001/IBDkin) program.
  - You will need to modify the hard-coded variables in the script.

## Citation

You need to cite any software you used for pre-processing!