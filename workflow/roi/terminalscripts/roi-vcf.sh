#!/bin/bash

# input
folder=$1
vcfin=$2
vcfout=$3
left=$4
right=$5
sampfile=$6
minaf=$7

# filtering
tabix -fp vcf ${vcfin}
bcftools view ${vcfin} -r 1:${left}-${right} -S ${sampfile} -q ${minaf} -Oz > ${vcfout}
