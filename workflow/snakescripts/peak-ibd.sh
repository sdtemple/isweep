#!/bin/bash

ibd=$1
window=$2
aboutfile=$3
cutoff=$4
jarprog=$5
index1=$6
index2=$7
index3=$8
genomeend=$(python -c "print(${9}*1000000)")
out=${10}

central=$(python ${aboutfile} ${window})
zcat ${ibd} | java -jar ${jarprog} 'I' -${index1} 0 ${cutoff} | java -jar ${jarprog} 'I' ${index2} 0 ${central} | java -jar ${jarprog} 'I' ${index3} ${central} ${genomeend} | gzip > ${out}
