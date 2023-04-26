#!/bin/bash

# First step: subsample phased vcf
# Seth D. Temple, sdtemple@uw.edu
# April 25, 2023

### fixed variables ###

# cluster
QUEUE=b-students.q
EMAIL=sdtemple@uw.edu
LOGS=~/logs/

# where
VCFS=/projects/browning/brwnlab/topmed_freeze8/merge_Mar21/
VCFPREFIX=merged_chr
VCFSUFFIX=_BAGS_BIOME_CCAF_FHS_HYPERGEN_JHS_MLOF_VTE_VUAF_WHI.PASS_polysnps.vcf.gz

# input variables
WHERE=$1
SOFTWARE=$2
CHRLOW=$3
CHRHIGH=$4
DIPLOIDSAMPLES=$5

mkdir -p ${LOGS}
mkdir -p $WHERE
mkdir -p ${WHERE}vcfs/

# subsample vcfs
for j in $(seq $CHRLOW 1 $CHRHIGH); do echo "bcftools view -Ob -S $DIPLOIDSAMPLES ${VCFS}${VCFPREFIX}${j}${VCFSUFFIX} | bcftools view -c 1:minor -Oz -o ${WHERE}vcfs/chr${j}.vcf.gz" | qsub -q ${QUEUE} -N subsamplevcf${j} -e ${LOGS} -o ${LOGS} -m e -M ${EMAIL}; done;
