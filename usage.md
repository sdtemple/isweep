## Do these methods apply to my sample size?

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
    - Use the different *_Ne() functions to make such an *.ne file
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

## Should I use these methods on biobank datasets > 100k individuals?

- No, unless you want to spend a lot on computing costs.
  - The methods should scale, but you'd wait awhile and spend a lot.
- If you downsample to ~10k, you should have more than enough data.
 - (If IBD **were** binomially distributed, downsampling keeps same distributional properties.)

## Should I use these methods for array data?

- We have not tested on array data.
- You could increase CHANGE:SIMULATE:MSPMAF to 0.01 or 0.05 and do a simulation study.
- It may that selection scan and coefficient estimation is okay.
- It may be that identifying candidate variants is not.

## Should I use these methods for low quality genotyping?

- We have not tested for CHANGE:SIMULATE:GTERR != 0.001,0.0002,0.0001.
- You could change this genotyping error rate and do a simulation study.
- See ibd-ends paper for discussion.

## Should I use these methods for non-diploid species?

- We do not have plans to extend methodology in this direction.
- Math in isweep python package would be fine.
 - You have to find another way to phase and infer IBD segments.
 - You have to find another way to identify candidate variants, estimate frequency.
 - See rule `third_hap_infer` in `rules/third.smk` for example. 