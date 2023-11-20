#!/bin/bash

CHR=$1
MAF=$2
COHORT=$3
PRIVST=$4
INFOLDER=/home/users/sdtemple/brwn/seth/topmed/${COHORT}
QUEUE=b-students.q
EMAIL=sdtemple@uw.edu

echo "bcftools view -Oz --min-af ${MAF}:nonmajor ${INFOLDER}/${COHORT}.${PRIVST}.chr${CHR}.vcf.gz > ${INFOLDER}/${COHORT}.${PRIVST}.chr${CHR}maf${MAF}.vcf.gz" | qsub -q ${QUEUE} -N ${COHORT}MAF${CHR} -m e -M ${EMAIL}
