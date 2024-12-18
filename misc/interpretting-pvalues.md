## When the p-values make sense

The p-values are calculated assuming that the IBD rates are normally distributed. This assumption is only reasonable when the Temple and Thompson asymptotic conditions hold.

Please check the `zhistogram.png` and related files after running scans.
- If these histograms do not look Gaussian, do not use or interpret the output column `PVALUE`!
    - You can still investigate regions with high IBD rates. Just don't use the language "statistically significant". 
- These histograms may not look Gaussian because:
    - Truncations / large point mass at/near 0.
        - Sample size is too small for the scaled population size.
        - Or telomeric/centromeric effects.
    - Pervasive forms of selection genome-wide.
    - Chromosomes subsets with different average/median IBD rates.
        - See related Markdown file. 
    - Poor genetic recombination map. 