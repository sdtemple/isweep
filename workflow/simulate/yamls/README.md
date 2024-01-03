#### Overall comments

Here I describe the *.yaml files in this directory.
These files operate like nested Python dictionaries.
There are two top levels: CHANGE and FIXED.
- You should make edits to CHANGE given your local setup and simstudy design.
- You should not make edits to FIXED w/o expertise or asking for help (git issues).

#### Description of CHANGE parameters

- FOLDERS (paths to key files)
    - MACRO (name of folder where simstudy results are)
    - MICRO (name of experiments files)
    - YAML (name of this file, copy + pasted for reference in MACRO folder)
    - TERMINALSCRIPTS (where some scripts are, should be in like this workflow/simulate/terminalscripts)
    - SOFTWARE (where required software is, you download these yourself)

- PROGRAMS
    - SLiM (where your slim executable is, https://messerlab.org/slim/)
    - HAPIBD (where your hap-ibd.jar executable is, https://github.com/browning-lab/hap-ibd)
    - IBDENDS (where your ibd-ends.jar executable is, https://github.com/browning-lab/ibd-ends)
    - GTERR: (where your add-uniform-err.jar executable, http://faculty.washington.edu/browning/)
    - FILTER: (where your filter-lines.jar executable is, http://faculty.washington.edu/browning/)

- CLUSTER
    - LARGEMEM (in gigabytes, depends on your compute nodes, due to slim, msprime, ibd-ends for some demography examples)

- SIMULATE
    - tNe (*.ne file formatted like ibdne.jar; at least 500 generations)
    - SAMPSIZE : 5000 (number of diploids)
    - MSPMAF : 0.001 (minimum allele frequency saved from msprime)
    - MIGRRATE : 0 (uniform migration rate)
    - NUMSUBPOP : 1 (number of subpopulations)
    - RHO : 1e-8 (recombination rate)
    - MU : 1e-8 (mutation rate)
    - GCLEN : 300 (mean length of gene conversion tract)
    - GCPROP : 2/3 (proportion of recombinations that are gene conversion)
    - CMLEN : 8 (region size in megabases)
    - LOC : 3,999,999 (location of adaptive mutation in basepair)
    - GTERR : 0.0002 (a rate of genotyping or other errors)
    - TSPLIT : 1000 (time backwards of population split)

Further comments
- MSPMAF will affect the allele frequencies in all vcf files
- You should only change MIGRRATE and NUMSUBPOP if you want to study cryptic population substructure
- You will throw errors if TSPLIT is not greater than the greatest time value in your experiments file
- Have only worked w/ uniform recombination maps

#### Description of FIXED parameters
- IBDENDS, CANDHAPIBD, HAPIBD
    - Names match those of options in ibd-ends.jar, hap-ibd.jar
    - CANDHAPIBD is the options for calling candidate IBD segments ---> ibd-ends.jar
        - You want to be liberal with this to not prematurely throw out pairs
        - ibd-ends.jar will adjust the endpoints

- SIMULATE
    - a,b : min,max allele frequency of adaptive mutation
        - isweep is designed for partial sweeps in [0.1,0.9]
    - Have not done much testing of NUMVAR (more than 1 mutation)
    - Have not done much testing  of SAMPLEONE, SAMPLEEQUAL (don't remember what they do)
    - PLOIDY should be 2;  haven't modified scripts to handle more

- ISWEEP
    - GENOMEEND (some large value in Mb to extend past the entire chromosome)
    - BY (computing IBD rates every BY kilobasepairs)
    - TELOCUT : 0.5 (cM at end of region (telomeres) to cutoff)
    - MLECUTOFF : 3.0 (cM threshold for selection coefficient estimation)
    - SCANCUTOFF: 2.0 (cM threshold for IBD rate selection scan)
    - DIAMETER : 6 (this number / 2 controls the maximum ibd path lengths)
    - GROUPCUTOFF : 3 (number of standard deviations to exceed in anomaly detection)
        - https://en.wikipedia.org/wiki/68%E2%80%9395%E2%80%9399.7_rule
    - WINSIZE : 250000 (size and step for location initliazing, in basepairs)
    - WINSTEP : 50000
    - HAPSIZE: 500000 (size and step for positions in haplotype rectangles)
    - HAPSTEP: 50000
    - FREQSIZE: 0.100 (size and step for frequencies in haplotype rectangles)
    - FREQSTEP: 0.025
    - NUMSNP: 5 (min count of snps in haplotype rectangle)
    - MINAAF: 0.10 (min allele frequencies considered)
    - MAXSPACING : 25000 (max positions w/o a snp)
    - QRANGE : 10 (num of percentiles to sum over in Qw statistic)


Further comments
- NUMSNP, MAXSPACING will probably depend on changes in mutation rate, minimum allele frequency
- SCANCUTOFF probably doesn't matter if you're doing a simstudy on selection coefficient estimation
