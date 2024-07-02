#!/bin/bash
pyscript=$1
mainfolder=$2
prefix=$3
derfile=$4
out=$5
python $pyscript --RelateSamples $mainfolder/$prefix.newick --DerivedFile $mainfolder/$derfile --out $mainfolder/$out
