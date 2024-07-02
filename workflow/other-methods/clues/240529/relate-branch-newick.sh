#!/bin/bash
relatefolder=$1
mainfolder=$2
prefix=$3
coal=$4
outfix=$5
mu=$6
nummcmc=$7
loc=$8
$relatefolder/scripts/SampleBranchLengths/SampleBranchLengths.sh --input $mainfolder/$prefix --output $mainfolder/$outfix.temp -m $mu --coal $coal --format n --num_samples $nummcmc --first_bp $loc --last_bp $loc
mv $mainfolder/$outfix.temp.anc $mainfolder/$outfix.anc
mv $mainfolder/$outfix.temp.mut $mainfolder/$outfix.mut
mv $mainfolder/$outfix.temp.newick $mainfolder/$outfix.newick
rm $mainfolder/$outfix.temp.anc || true
rm $mainfolder/$outfix.temp.mut || true
rm $mainfolder/$outfix.temp.newick || true
rm $mainfolder/$outfix.temp.dist || true
# remove temp files, folders in directory you sent slurm from
