"""
Put mutations on tree sequence

@author: sdtemple
"""

import sys, os, subprocess, pyslim, msprime, tskit
from isweep import *

mu=snakemake.config['FIXED']['MU']
rho=snakemake.config['FIXED']['RHO']
Ne=snakemake.config['FIXED']['tNe']
trees=snakemake.input.trees
maf=snakemake.config['FIXED']['MSPMAF']
n=snakemake.config['FIXED']['SAMPSIZE']
ploidy=snakemake.config['FIXED']['PLOIDY']
filename=snakemake.output.bcf

mu = float(mu)
rho = float(rho)
maf = float(maf)
n = int(float(n))
ploidy = int(float(n))
N = n * ploidy
mac = int(N * maf)
mac = min(2, mac)
Ne = read_Ne(Ne)
Ne = int(float(list(Ne.values())[-1]))

orig_ts = tskit.load(trees)
orig_ts = pyslim.update(orig_ts) # pyslim, tskit new versions

# random nucleotide(s)
orig_ts = pyslim.generate_nucleotides(orig_ts)
orig_ts = pyslim.convert_alleles(orig_ts)

# recapitate
recap_ts = pyslim.recapitate(orig_ts, ancestral_Ne = Ne, recombination_rate = rho) # address warnings

# root check
roots = max(tree.num_roots for tree in recap_ts.trees())
if roots > 1:
    raise ValueError('More than 1 root')

# add neutral mutations
mut_ts = msprime.sim_mutations(recap_ts, rate=mu, keep=True)

# remove rare variants
# courtesy of Ryan Waples (waplesr@uw.edu)
def strip_MAC(ts, MAC):
    """
    Removes sites with minor allele count <= MAC
    Returns a new tree sequence with sites removed
    """
    initial_sites = ts.num_sites
    samples = ts.samples()
    nsamp = len(samples)
    sites_to_remove = []
    for tree in ts.trees():
        for site in tree.sites():
            if len(site.mutations) == 1:
                mut = site.mutations[0]
                if (tree.num_samples(mut.node) <= MAC) or \
                    (tree.num_samples(mut.node) >= (nsamp - MAC)):
                    sites_to_remove.append(site.id)
            else:
                sites_to_remove.append(site.id)
    ts = ts.delete_sites(sites_to_remove)
    final_sites = ts.num_sites
    print(f"MAC filter (<={MAC}):")
    print(f"removed {initial_sites - final_sites} sites ({(initial_sites - final_sites) / (initial_sites):.0%}), {final_sites} sites remain")
    return ts
def strip_adjacent_sites(ts, dist=1.5):
    """Remove sites within dist bp of each other.
    Removes the right-most site in each pair.
    Returns a new tree-sequence.
    """
    initial_sites = ts.num_sites
    present_samples = np.intersect1d(
    	ts.samples(),
    	np.where(ts.tables.nodes.asdict()['time'] == 0)[0]
    )
    sites_to_remove = []
    prev_pos = 0
    for tree in ts.trees(tracked_samples=present_samples):
        for site in tree.sites():
            if len(site.mutations) == 1:
                pos = site.position
                if (pos - prev_pos) < dist:
                    sites_to_remove.append(site.id)
                    prev_pos = pos
            else:
                sites_to_remove.append(site.id)
    ts = ts.delete_sites(sites_to_remove)
    final_sites = ts.num_sites
    print('Adjacent sites filter:')
    print(f"removed {initial_sites-final_sites} sites ({(initial_sites-final_sites) / (initial_sites):.0%}), {final_sites} sites remain")
    return(ts)

# apply MAC, bp distance filter
ts = strip_MAC(mut_ts, mac)
ts = strip_adjacent_sites(ts)
# ts = strip_adjacent_sites(ts, 10)

# saving

read_fd, write_fd = os.pipe()
write_pipe = os.fdopen(write_fd, "w")
with open(filename, "w") as bcf_file:
    proc = subprocess.Popen(
        ["bcftools", "view", "-O", "b"], stdin=read_fd, stdout=bcf_file
    )
    ts.write_vcf(write_pipe)
    write_pipe.close()
    os.close(read_fd)
    proc.wait()
    if proc.returncode != 0:
        raise RuntimeError("bcftools failed with status:", proc.returncode)
