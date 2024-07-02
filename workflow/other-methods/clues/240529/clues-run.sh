#!/bin/bash
clues=$1
dafscript=$2
folder=$3
prefix=$4
timesfile=$5
daffile=$6
coalfile=$7
dominancecoef=$8
df=$9
tCutoff=${10}
daffile=$folder/$daffile
allelefreq=$(python $dafscript $daffile)
echo $allelefreq
timesfile=$folder/$timesfile
python $clues --times $timesfile --noAlleleTraj --coal $coalfile --popFreq $allelefreq --out $folder/$prefix --h $dominancecoef --df $df --tCutoff ${tCutoff}
