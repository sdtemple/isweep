#### This document is for the developer(s)

- msprime recapitation error in complex SLiM setups
    - Use addSubpopSplit() instead of addSubpop()
    - Included `+1` in the `extend_Ne` function
- Writing `plots/` and `stats/`
    - File name was the mode of last region of interest
    - Simple typo of "meBP" versus "moBP" variable name
- Double checked that dominance and recessive selection were implemented correctly.
- Incorporated opportunity for more MCMC, thinning, and burn-in for CLUES pipe.