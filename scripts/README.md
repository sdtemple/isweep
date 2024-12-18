## Documenting python scripts

You may want to individually run some of the scripts. Here, I describe the scripts I imagine a user may want to re-run individually.

You can use `python *.py --help` to see the inputs to all scripts. For bash scripts, read in `vim` or `more`.

You can customize your plots after the pipeline. There are motivating scripts in 'plotting/'.

Other subfolders are organized by workflow (`scan/` versus `model/` sweeps) or general utilities.

### Scripts
---

#### `plotting/`
- `plot-scan.py` : makes a plot of scanning statistic.
    - Used to visualize the selection scans.
- `plot-histogram.py` : makes a plot of some column statistic.
    - Used as a diagnostic check to see genome-wide IBD rates.
- `plot-autocovariance.py` : makes a plot for the log linear model of IBD rate decay.
    - Used as a diagnostic check to see if IBD rate process is approximately OU.
- `plot-sweep.py` : line plot for allele frequency over time.
    - Used when sure that loci is undergoing a sweep.

#### `utilities/`
- `filter-lines.py` : apply simple filter to tab-separated data.
- `count_ibd.py` : compute number of IBD segments overlapping loci along genome.
- `interpolate-map.py` : make genetic map with rows every fixed cM step size.
- `trim-telomere-map.py` : cut some cM from start and end of genetic map.
- `exclude-samples.py` : finds outer intersection of two files.
- `remove-phase.py` : create an unphased vcf file.
- `add-uniform-err.py` : introduce allelic errors in vcf file.
- `fake-phenotypes.py` : create file with fake/randomly simulated categorical phenotypes.

#### `scan/`
- `scan.py` : implement a genome-wide scan over some statistic.
- `multiple-testing-analytical.py` : compute the analytical threshold for multiple testing.
- `multiple-testing-simulation.py` : compute the simulation-based threshold for multiple testing.
    - Requires theta estimate from running `multiple-testing-analytical.py` beforehand.

#### `model/`
- `estimate-norm.py` : estimate selection coefficient with standard normal intervals.
- `estimate-perc.py` : estimate selection coefficient with percentile intervals.
- `rank.py` : score alleles based on how differentiated they are from IBD outgroup and rest of sample.
- `outliers.py` : detect clusters with excess IBD rate.




