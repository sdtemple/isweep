
This is the workflow manager (slurm) that most people use in 2024.

- Commands for `slurm`
    - https://github.com/aws/aws-parallelcluster/wiki/Transition-from-SGE-to-SLURM
    - In the snakemake pipeline `--cluster "sbatch [options]" `
        - "-e ~/your-logs/{rule}.e" and "-o ~/your-logs/{rule}.o" will control where stderr, stdout go
        - "--cpus-per-task=XX" says how many cpus per job
        - "--nodes=XX" says how many nodes per job
        - "--partition=SOMENAME" says which partitions to use
        - "--mem=XX" says how much memory in MB
        - "--mail-type=ALL" and "--mail-user=your.email@university.edu" sends mail to you when a job finishes
        - "--job-name={rule}"

By and large, the Sun Grid systems (qsub) is not being supported anymore.

- Commands for `qsub`
    - In the snakemake pipeline `--cluster "qsub -q your-queue.q " `
    - "-V" is important to pass in your conda environment!
    - "-pe local XXX" is how many threads you will use
    - "-m e -M your.email@uni.edu"
        - You don't have to send emails to yourself if you don't want to
    - "-l h_vmem=XXXG" is how much hard RAM memory you will use
