### What to do

You can modify the python scripts here for your own simulation study.
This is based on simulating true IBD tracts from a coalescent.
You need to be in a conda environment with the isweep package installed.
You probably want to send large simulations to a cluster.

This is the first simulation study in our paper.
    - Results in this main folder are those in our first table
    - Other results are in scripts from nuanced/ and misspecified/

1. "vim [filename].py" or "head -n [filename].py"
    - This is to see the script inputs
    - In vim, 
        - use "Esc" to leave insert "shift + I"
        - use ":wq" then "Enter" to save and exit
2. Write some python/R scripts to summarize the data
    - summarize-simstudy.py may do something like that