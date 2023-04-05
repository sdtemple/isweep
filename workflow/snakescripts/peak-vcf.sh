#!/bin/bash

# input
vcfin=$1
windowin=$2
folderin=$3
pm=$4
vcfout=$5
maf=$6
method=$7

# leftstream, rightstream
left=$(python ${folderin}/ibd-focus-min.py ${windowin} ${pm})
right=$(python ${folderin}/ibd-focus-max.py ${windowin} ${pm})

# changing to bgzip
gunzip -c ${vcfin}.vcf.gz | bgzip  > ${vcfin}.${method}.vcf.bgz
tabix -fp vcf ${vcfin}.${method}.vcf.bgz

# filtering
bcftools view ${vcfin}.${method}.vcf.bgz -r 1:${left}-${right} -q ${maf}:nonmajor -Oz > ${vcfout}
