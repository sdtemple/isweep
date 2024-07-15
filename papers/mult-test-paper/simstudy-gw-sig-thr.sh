#!/bin/bash
rscr=$1
fileout=$2
theta=$3
pval=$4
stepsize=$5
chrnum=$6
chrlen=$7
numtrain=$8
numtype1=$9
numcores=${10}
Rscript --vanilla $rscr $fileout $theta $pval $stepsize $chrnum $chrlen $numtrain $numtype1 $numcores