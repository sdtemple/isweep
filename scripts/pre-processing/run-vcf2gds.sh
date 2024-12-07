#!/bin/bash

VCF=$1
GDS=$2
COHORT=$3
VFOLDER=/home/users/sdtemple/brwn/seth/topmed/${COHORT}
SFOLDER=/home/users/sdtemple/brwn/seth/topmed
IGN=chr
MET=biallelic.only
QUEUE=b-students.q
EMAIL=sdtemple@uw.edu

echo "Rscript --vanilla ${SFOLDER}/vcf2gds.R ${VFOLDER}/${VCF} ${VFOLDER}/${GDS} ${IGN} ${MET}" | qsub -q ${QUEUE} -N vcf2gds -m e -M ${EMAIL}
