#!/bin/bash

# arguments
MINCHR=1
MAXCHR=22
COHORT=phase.paper.data
MAP=38
IBDKINEXEC=/home/users/sdtemple/brwn/software/IBDkin/src-v2.8.7.7/IBDkin
INDFILE=/home/users/sdtemple/brwn/seth/topmed/${COHORT}/topmed.all.samples.txt
SEGFOLDER=/home/users/sdtemple/brwn/seth/topmed/${COHORT}
MAPFOLDER=/projects/browning/maps/decode.2019.b${MAP}
KINFOLDER=/home/users/sdtemple/brwn/seth/topmed/${COHORT}/ibdkin
AUTOSOMES=/home/users/sdtemple/brwn/seth/topmed/decode.autosomes.txt
mkdir -p ${KINFOLDER}
EMAIL=sdtemple@uw.edu
QUEUE=b-students.q

# create --range input
sed 's/chr//' ${AUTOSOMES} > ${KINFOLDER}/range.autosomes.txt
awk '{if($1 >= a && $1 <= b) {print $0}}' a=$MINCHR b=$MAXCHR ${KINFOLDER}/range.autosomes.txt > ${KINFOLDER}/range.range.txt
sed -e 's/^/chr/' ${KINFOLDER}/range.range.txt > ${KINFOLDER}/range.txt
rm ${KINFOLDER}/range.autosomes.txt
rm ${KINFOLDER}/range.range.txt

# create --ibdfile inputs
for i in $(seq $MINCHR $MAXCHR);
do
    echo "${SEGFOLDER}/chr${i}.ibd.gz";
done > ${KINFOLDER}/ibdfile.txt

# create --map file
for i in $(seq $MINCHR $MAXCHR);
do
    cat "${MAPFOLDER}/decode2019.chrchr${i}.GRCh${MAP}.map";
done > ${KINFOLDER}/decode.2019.chr${MINCHR}-${MAXCHR}.b${MAP}.map

# ibdkin run
echo "${IBDKINEXEC} --ibdfile ${KINFOLDER}/ibdfile.txt --map ${KINFOLDER}/decode.2019.chr${MINCHR}-${MAXCHR}.b${MAP}.map --ind ${INDFILE} --range ${KINFOLDER}/range.txt --nthreads 12 --out ${KINFOLDER}/chr${MINCHR}-${MAXCHR} --outmask --outcoverage" | qsub -q ${QUEUE} -pe local 12 -N ibdkin -m e -M ${EMAIL} -cwd

