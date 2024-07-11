#!/bin/bash
module load Java
argon=$1
prefix=$2
samplesize=$3
effectivesize=$4
chromosomelength=$5
morganlength=$6
blocksize=$7
recombinationrate=$8
mutationrate=$9
memory=${10}
java -Xmx${memory}g -jar $argon \
    -N $effectivesize \
    -pop 1 $samplesize \
    -size $chromosomelength \
    -rec $recombinationrate \
    -mut $mutationrate \
    -seq 1.0 \
    -IBD $morganlength \
    -len $blocksize \
    -gz \
    -out $prefix  
