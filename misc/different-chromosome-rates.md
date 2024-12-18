## Subsets of chromosomes have vastly different average/median IBD rates

The hypothesis test has the null model of a common genome-wide mean IBD rate. You may observe that subsets of the chromosomes have visibly different average IBD rates (in the selection scan figure). (This may reflect different chromosome-specific recent Ne's.) 

Consider the following example:

- Chromosome 1 has average IBD rate 1e-4
- Chromosome 2 has average IBD rate 1e-4
- Chromosome 3 has average IBD rate 3e-4
- Chromosome 4 has average IBD rate 3e-4

You should analyze chromosomes 1 and 2 separately from chromosomes 3 and 4.

In this case, create to configuration files by adjusting the CHRLOW and CHRHIGH parameters.

The CHANGE:ISWEEP:GENOMESIZE parameter should still reflect the total sum of cM in chromosomes 1 through 4. The parameter controls the effective number of hypothesis tests.

