#!/bin/bash
pyscript=$1
fileout=$2
prefix=$3
suffix=$4
st=$5
en=$6
slide=$7
columnname=$8
covariancelength=$9
chrlen=${10}
pval=${11}
python $pyscript \
    $fileout \
    $prefix \
    $suffix \
    $st \
    $en \
    $slide \
    $columnname \
    $covariancelength \
    $chrlen \
    $pval