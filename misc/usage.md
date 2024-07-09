## Do these methods apply to my study and dataset?

### Genetic maps

If you don't have a genetic maps, assume a uniform recombination rate.

If you don't have an estimate of the uniform recombination rate, the methods do not apply.

To make a uniform map, you need four columns without a header.
- Let 1e-8 be the the recombination rate. 
- Let there be 100 Mb of data on a chromosome.
- Let chr1 be the name of the chromosome
- The first row is `chr1  . 0 1`
- The second row is `chr1 . 100 100000000`

Make a uniform rate genetic map for each chromosome.

### Sample size

It'd be nice if you have:
- Ancestral Ne >= 3000
- No population bottleneck with Ne < 3000 in last 500 generations
- Sample from 1 generation
  - I.e., use other methods for flies, mosquitoes, worms

#### You could test if there is enough IBD:

- In python, use the `simulate_ibd` and `read_Ne` functions in `isweep`
  - Ne file is tab-separated w/ two columns: generation, Ne
  - You should make at least 500 generations
  - It is okay to have a guess of Ne size before inferring with IBDNe
    - Use the different `*_Ne()` functions to make such an *.ne file
- Do you have on average # of 3.0 cM ibd tracts >= 200?
  - Then, yes, you probably have enough samples to run this analysis.

#### Some species from `stdpopsim` that may apply:

- Great apes
  - Pongo abelii
  - Pan troglodytes
  - Papio anubis
- Anas platyrhynchos (duck)
- Apis mellifera (bee)
- Canis familiaris (dog)
  - Address population structure, admixture

Another species may be Plasmodium falciparum (https://www.nature.com/articles/s41467-024-46659-0)

### Should I use these methods on biobank datasets > 100k individuals?

- No, unless you want to spend a lot on computing costs.
  - The methods should scale, but you'd wait awhile and spend a lot.
- If you downsample to ~10k, you should have more than enough (human) data.

### Can I use these methods for array data?

- We have **not tested** on array data.
  - If you wish to try, we recommend using the `hap-ibd` array parameters (see their paper).
- You could increase CHANGE:SIMULATE:MSPMAF to 0.01 or 0.05 and do a simulation study `workflow/simulate/`.

### Can I use these methods for low quality genotyping?

- We have **not tested** for CHANGE:SIMULATE:GTERR > 0.001.
- You could change this genotyping error rate and do a simulation study `workflow/simulate/`.
- See `ibd-ends` paper for discussion.

### Can I use these methods for non-diploid species?

- We **do not** have plans to extend methodology in this direction.
- Math in `isweep` Python package would be fine.
 - You have to find another way to phase and infer IBD segments.
 - You have to find another way to identify candidate variants and estimate frequency.
 - See rule `third_hap_infer` in `rules/third.smk` for example. 

## Modeling choices

### Can I assume/fit a non-additive model?

You do not have to assume additivity. There are 4 models available.
- 'a' Additive 1:(1+s):(1+2s)
- 'm' Multiplicative 1:(1+s):(1+s)^2
- 'd' Dominant 1:(1+s):(1+s)
- 'r' Recessive 1:1:(1+s)

The additive model is the default.

You modify this in your ROI tab-separated file for `workflow/roi/`.

### Can I assume a sweep from de novo mutation?

Yes, review the standing variation parameter in the above `help()`.

You can define a frequency at which selection began.

We do not provide utility for this in the pipelines. 

You can copy and modify the script `workflow/roi/scripts/estimate.py`.

### Can I assume an ongoing sweep?

Yes, review the `tau0` parameter in the above `help()`.

You can define a recent generation at which selection ends (s = 0).

We do not provide utility for this in the pipelines. 

You can copy and modify the script `workflow/roi/scripts/estimate.py`.

If `tau0 >> 15`, you could bias your estimate. Or, there isn't much IBD to learn from.