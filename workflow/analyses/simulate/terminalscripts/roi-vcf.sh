#!/bin/bash

# input
folder=$1
vcfin=$2
vcfout=$3
left=$4
right=$5
chrnum=$6
maf=$7

# changing to bgzip
gunzip -c ${vcfin} | bgzip  > ${folder}/chr${chrnum}.vcf.bgz
tabix -fp vcf ${folder}/chr${chrnum}.vcf.bgz

# filtering
bcftools view ${folder}/chr${chrnum}.vcf.bgz -r 1:${left}-${right} -q ${maf}:nonmajor -Oz > ${vcfout}
