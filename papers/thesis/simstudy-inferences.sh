#!/bin/bash
pyscript=$1
fileout=$2
selcoef=$3
allelefreq=$4
morgan=$5
sampsize=$6
popsize=$7
numsim=$8
numboot=$9
python $pyscript $numsim $numboot $sampsize $allelefreq $popsize $fileout $selcoef $morgan
