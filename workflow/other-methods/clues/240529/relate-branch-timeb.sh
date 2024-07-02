#!/bin/bash
relatefolder=$1
mainfolder=$2
prefix=$3
coal=$4
outfix=$5
mu=$6
nummcmc=$7
loc=$8
$relatefolder/scripts/SampleBranchLengths/SampleBranchLengths.sh --input $mainfolder/$prefix --output $mainfolder/$outfix.temp -m $mu --coal $coal --format b --num_samples $nummcmc -first_bp $loc --last_bp $loc
cp $mainfolder/$outfix.temp.anc $mainfolder/$outfix.anc
cp $mainfolder/$outfix.temp.mut $mainfolder/$outfix.mut
cp $mainfolder/$outfix.temp.timeb $mainfolder/$outfix.timeb
rm $mainfolder/$outfix.temp.anc || true
rm $mainfolder/$outfix.temp.mut || true
rm $mainfolder/$outfix.temp.timeb || true
rm $mainfolder/$outfix.temp.dist || true
# remove temp files, folders in directory you sent slurm from
