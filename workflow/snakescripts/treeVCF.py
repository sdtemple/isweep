"""
Put mutations on tree sequence

@author: sdtemple
"""

import sys, os, subprocess, pyslim, msprime, tskit
from iSWEEP import *

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
# print('UPDATE: recapitated tree sequence')

# root check

roots = max(tree.num_roots for tree in recap_ts.trees())
if roots > 1:
    raise ValueError('More than 1 root')

# add neutral mutations

mut_ts = msprime.sim_mutations(recap_ts, rate=mu, keep=True)
# print('UPDATE: added neutral mutations')

# remove rare variants
# courtesy of Ryan Waples (waplesr@uw.edu)

def simple_strip_MAC(ts, MAC):
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
    ts = ts.delete_sites(sites_to_remove)
    final_sites = ts.num_sites
    # print(f"MAC filter (<={MAC}):")
    # print(f"removed {initial_sites - final_sites} sites ({(initial_sites - final_sites) / (initial_sites):.0%}), {final_sites} sites remain")
    return ts

ts = simple_strip_MAC(mut_ts, mac)
# print('UPDATE: applied MAC filter')

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
