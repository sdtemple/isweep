## By and large, do not analyze super small chromosomes

In chromosomes <= 10.0 cM, many IBD tracts could be **truncated** at the left and right ends of the genetic map.

For instance, if an IBD segment starts 1.5 centiMorgans (cM) from the right end of the chromosome, it can only be maximum 1.5 cM long, because of **truncation**.

In general, `hap-ibd.jar` and `ibd-ends.jar` are less accurate at IBD segment detection near the telomere.

We recommend you do not analyze chromosomes <= 10 cM.

The automated pipeline requires contiguous ordering of chromosomes. You should use the `ln -s [current-vcf-path] [new-vcf-path]` and `ln -s [current-map-path] [new-map-path]` to create pseudo-chromosome orderings.

You will need to adjust the CHANGE:ISWEEP:GENOMESIZE parameter, which is the total cM length of sequences analyzed, to reflect the removal of small chromosomes from the analysis. The parameter controls the effective number of hypothesis tests.

### For example,
---

Sizes of chromosome
- Chromosome 1 is 20 cM
- Chromosome 2 is 5 cM
- Chromosome 3 is 15 cM
- Chromosome 4 is 25 cM

Adjust names of vcf files.
- `ln -s chr1.vcf.gz [new-vcf-folder]/chr1.vcf.gz`
- `ln -s chr3.vcf.gz [new-vcf-folder]/chr2.vcf.gz`
- `ln -s chr4.vcf.gz [new-vcf-folder]/chr3.vcf.gz`

Adjust names of genetic maps
- `ln -s chr1.map [new-map-folder]/chr1.map`
- `ln -s chr3.map [new-map-folder]/chr2.map`
- `ln -s chr4.map [new-map-folder]/chr3.map`

### Modeling sweeps
---

When you create these pointer files, you may create an issue with `workflow/model-selection`. This issue occurs because of region subsetting with `bcftools query`. For sweep estimation, you may need to:
- Rename the chromosome number in the ROI output file to be the true chromosome number
- Reset the VCF folders in the YAML file to by the folders with the true chromosome numbers
- Adjust the IBD data in 'ibdsegs/ibdends/mle' to reflect the true chromosome numbers. 


