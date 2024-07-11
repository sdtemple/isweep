import msprime
import numpy as np
import sys
fileout,sample_size,recombination_rate,sequence_length,dftw_duration=sys.argv[1:]
nsamp = int(float(sample_size))
demography=25000
print(demography)
ts = msprime.sim_ancestry(
  samples = nsamp,
  ploidy=2,
  population_size = demography,
  model=[
    msprime.DiscreteTimeWrightFisher(duration=int(float(dftw_duration))),
    msprime.StandardCoalescent(),
  ],
  recombination_rate=float(recombination_rate),
  sequence_length=int(float(sequence_length)),
#   gene_conversion_rate=2e-8,
#   gene_conversion_tract_length=300,
#   random_seed = 1234
)
ts.dump(fileout)
