# Put mutations on tree sequence


import sys, os, subprocess, pyslim, msprime
import tskit as ts
import numpy as np
import pandas as pd

filename, trees, mu, rho, ancNe, maf, n, ploidy = sys.argv[1:]

mu = float(mu)
rho = float(rho)
maf = float(maf)
n = int(float(n))
ploidy = int(float(n))
N = n * ploidy
mac = int(N * maf)
mac = min(2, mac)
Ne = int(float(ancNe))

orig_ts = ts.load(trees)
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
def strip_MAC(treeseq, MAC):
    """
    Removes sites with minor allele count <= MAC
    Returns a new tree sequence with sites removed
    """
    initial_sites = treeseq.num_sites
    samples = treeseq.samples()
    nsamp = len(samples)
    sites_to_remove = []
    for tree in treeseq.trees():
        for site in tree.sites():
            if len(site.mutations) == 1:
                mut = site.mutations[0]
                if (tree.num_samples(mut.node) <= MAC) or \
                    (tree.num_samples(mut.node) >= (nsamp - MAC)):
                    sites_to_remove.append(site.id)
            else:
                sites_to_remove.append(site.id)
    treeseq = treeseq.delete_sites(sites_to_remove)
    final_sites = treeseq.num_sites
    print(f"MAC filter (<={MAC}):")
    print(f"removed {initial_sites - final_sites} sites ({(initial_sites - final_sites) / (initial_sites):.0%}), {final_sites} sites remain")
    return treeseq
def strip_adjacent_sites(treeseq, dist=1.5):
    """Remove sites within dist bp of each other.
    Removes the right-most site in each pair.
    Returns a new tree-sequence.
    """
    initial_sites = treeseq.num_sites
    present_samples = np.intersect1d(
    	treeseq.samples(),
    	np.where(treeseq.tables.nodes.asdict()['time'] == 0)[0]
    )
    sites_to_remove = []
    prev_pos = 0
    for tree in treeseq.trees(tracked_samples=present_samples):
        for site in tree.sites():
            if len(site.mutations) == 1:
                pos = site.position
                if (pos - prev_pos) < dist:
                    sites_to_remove.append(site.id)
                    prev_pos = pos
            else:
                sites_to_remove.append(site.id)
    treeseq = treeseq.delete_sites(sites_to_remove)
    final_sites = treeseq.num_sites
    print('Adjacent sites filter:')
    print(f"removed {initial_sites-final_sites} sites ({(initial_sites-final_sites) / (initial_sites):.0%}), {final_sites} sites remain")
    return(treeseq)

# apply MAC, bp distance filter
treeseq = strip_MAC(mut_ts, mac)
treeseq = strip_adjacent_sites(treeseq)

# saving

read_fd, write_fd = os.pipe()
write_pipe = os.fdopen(write_fd, "w")
with open(filename, "w") as bcf_file:
    proc = subprocess.Popen(
        ["bcftools", "view", "-O", "b"], stdin=read_fd, stdout=bcf_file
    )
    treeseq.write_vcf(write_pipe)
    write_pipe.close()
    os.close(read_fd)
    proc.wait()
    if proc.returncode != 0:
        raise RuntimeError("bcftools failed with status:", proc.returncode)
