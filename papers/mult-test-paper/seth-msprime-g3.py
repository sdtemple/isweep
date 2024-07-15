import msprime
import numpy as np
import sys
fileout,sample_size,recombination_rate,sequence_length,dftw_duration=sys.argv[1:]
nsamp = int(float(sample_size))
# population configuration
N1 = 5000 * np.exp(0.01*240) * np.exp(0.07*50) * np.exp(0.15*10) # Population size at present
N2 = 5000 * np.exp(0.01*240) * np.exp(0.07*50)
N3 = 5000 * np.exp(0.01*240)
N4 = 5000
T1 = 0
T2 = 10
T3 = 60
T4 = 300
growth_rate1 = -np.log(N2/N1) / (T2-T1)
growth_rate2 = -np.log(N3/N2) / (T3-T2)
growth_rate3 = -np.log(N4/N3) / (T4-T3)
demography = msprime.Demography()
demography.add_population(initial_size=N1, growth_rate=growth_rate1)
demography.add_population_parameters_change(time=T2, growth_rate=growth_rate2)
demography.add_population_parameters_change(time=T3, growth_rate=growth_rate3)
demography.add_population_parameters_change(time=T4, growth_rate=0)
print(demography)
ts = msprime.sim_ancestry(
  samples = nsamp,
  ploidy=2,
  demography = demography,
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