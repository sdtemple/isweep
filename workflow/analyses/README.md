These folders contain code to reproduce simulation studies and real data analyses for the paper.

coalsimstudy:
  Simulation studies on coalescent IBD data
  Uses simulation scripts in isweep paradigm
  These simulation studies base the methodology on robust statistical inference
  And/or explainable biases in model-based inference under model misspecification
  Python files developed with Copy and Paste, modifying accordingly
  Separate README.md for results of coalescent IBD studies
  
forwsimstudy
  * Workhorse simulation studies based on genetic data
  * Forward simulate 1 de novo sweep variant with SLiM
  * Recapitate and place mutations on coalescent tree with tskit, msprime
  * Perform IBD-based inference steps
  * Use isweep paradigm to conduct inference
    * Create a separate conda environment for iSAFE comparison
    * Create a separate conda environment for CLUES comparison
      * These methods are from 2018-2020 and may downgrade versions required for isweep
      
isweep inference on real genetic data
  realprepstudy
    Initial procedures to identify regions of interest
      Some data preprocessing
      Genome-wide IBD calls for IBDNe inference
      Genome-wide IBD calls for hap-ibd + ibd-ends
      Manhattan plot style inference
      Some plotting
  realpoststudy
    Estimation procedures in the isweep paradigm
      Infer p(0) for allele frequency of putatively causal mutant
      Method-of-moments and parametric bootstrap for s | p(0), Ne(t), ...
      Some tables
    Do this after the realprepstudy, identifying regions of interest!!!
