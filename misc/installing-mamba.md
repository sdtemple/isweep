
1. `wget https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-x86_64.sh` 
2. `bash Miniforge3-Linux-x86_64.sh`
3. `mamba` (check that mamba works)
    - You may need to `vim .bashrc`.
    - Put in a line `alias mamba="/path/to/miniforge3/bin/mamba"`.
    - `source .bashrc`
    - Sign out / exit terminal
    - Sign back in / begin terminal 
4. `python` (check that python works)
    - You may need to do the same for a python version under `/path/to/miniforge3/bin/`.
    - `exit()`
5. `mamba env create -f isweep-environment.yaml # make environment` (inside the git repo)
6. `mamba activate isweep`
7. `bcftools`
8. `tabix`
9. `snakemake -c1 -n --configfile YOUR.yaml`

