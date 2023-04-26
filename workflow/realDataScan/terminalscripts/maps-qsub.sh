#!/bin/bash

# First step: copy and paste genetic maps
# Seth D. Temple, sdtemple@uw.edu
# April 25, 2023

### fixed variables ###

# cluster
QUEUE=b-students.q
EMAIL=sdtemple@uw.edu
LOGS=~/logs/

# where
MAPS=/projects/browning/maps/decode.2019.b38/
MAPPREFIX=decode2019.chrchr
MAPSUFFIX=.GRCh38.map

# input variables
WHERE=$1
SOFTWARE=$2
CHRLOW=$3
CHRHIGH=$4

mkdir -p ${LOGS}
mkdir -p $WHERE
mkdir -p ${WHERE}maps/

# copy and paste maps
for j in $(seq $CHRLOW 1 $CHRHIGH); do echo "cp ${MAPS}${MAPPREFIX}${j}${MAPSUFFIX} ${WHERE}maps/chr${j}.map" | qsub -q $QUEUE -N map${j} -e ${LOGS} -o ${LOGS}; done;
