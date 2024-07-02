#!/bin/bash
mainfolder=$1
prefix=$2
vcffile=$3
left=$4
right=$5
chrnum=$6
subsample=$7
vcffile=$mainfolder/$vcffile
out=$mainfolder/$prefix.vcf
bcftools view $vcffile -r $chrnum:$left-$right -S $subsample -Ov -o $out 
