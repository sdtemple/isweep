import sys
import pyslim, msprime
import tskit as ts

treesin, treesout, rho, ancNe = sys.argv[1:]
rho = float(rho)
Ne = int(float(ancNe))

# load
orig_ts = ts.load(treesin)
orig_ts = pyslim.update(orig_ts) # pyslim, tskit new versions

# recapitate
recap_ts = pyslim.recapitate(orig_ts, ancestral_Ne = Ne, recombination_rate = rho) # address warnings

# dump/save
recap_ts.dump(treesout)