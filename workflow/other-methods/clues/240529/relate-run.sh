#!/bin/bash
relatefolder=$1
mainfolder=$2
prefix=$3
mu=$4
coal=$5
genmap=$6
$relatefolder/bin/Relate --mode All --haps $mainfolder/$prefix.haps --sample $mainfolder/$prefix.sample -m $mu --coal $coal -o $prefix.temp --map $genmap || true
cp $prefix.temp.anc $mainfolder/$prefix.anc
cp $prefix.temp.mut $mainfolder/$prefix.mut
rm $prefix.temp.anc || true
rm $prefix.temp.mut || true
