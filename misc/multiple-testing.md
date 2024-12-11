## Updates with multiple-testing corrections

#### Asymptotic normality assumption
---

The scan conducts one sample Z tests. The hypothesis test is reasonable when the Temple and Thompson conditions hold.
- Large recent population sizes (Ne)
- Large detection threshold (>= 2.0 cM)
- Large sample size 
- Recent population sizes times detection threshold are much larger than sample size.
- Sample size squared is much larger than recent population sizes times detection threshold.

You should evaluate the `zhistogram.png` file to see if your results "appear normally distributed".
- Make sure it is symmetric.
- It is okay if there are very small lower tails and very large upper tails.
    - Putatively places where IBD segments are poorly detected.
    - Putatively genomic regions under selection (an alternative model).

#### Model assumption
---

The multiple-testing correction is based on modeling the IBD rates as a stochastic process. The assumptions are:
- IBD rates under null model are normally distributed. (Hence, verifying that the Temple and Thompson conditions approximately hold).
- Decay of IBD rates is well fit as exponental. Look at the plot `autocovariance.png`.
- There is approximate "stationary" in mean and variance.
- There is approximate "first-order Markov" property.
- If your estimate of theta is much smaller than 20, your test may have low power.
- The test may be slightly anti-conservative using the >= 2.0 cM threshold.
    - View genomic regions barely exceeding the GW significance level with skepticism.
- The test may be slightly conservative using the >= 3.0 cM threshold.
- p values are based on assumption that test statistic is normally distributed.

#### Deviation from mean
---

The scan is looking for deviations from the mean IBD rate.
- If some chromosomes have uniformly inflated means, you need to analyze these subsets separately.
- Make two *.yaml files with different `CHRLOW` and `CHRHIGH` versions.
- Adjust the total cM length parameter accordingly.

#### What these changes mean for the pipeline

- Regions of interest file is based on the analytical multiple-testing correction. You can review the heuristic four standard deviations and simulation-based multiple-testing correction as well in the file `fwer.analytical.txt` and `fwer.simulation.txt` (the simulated distribution of first hitting time).
- Significance must now last for 0.5 cM continguous instead of 1.0 cM.
- IBD rate is computed every fixed cM step size. (We interpolate the genetic map.)
    - The old heuristic calculated variance every fixed bp step size.
- Genome-wide significance level (p value) is reasonable when Temple and Thompson conditions hold.
