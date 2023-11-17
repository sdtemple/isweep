
ibdfile=$2
freqfile=$3
cm=$4
n=$5
ploidy=2
effdemo=$6
fileout=$1

ibdest=$(zcat $ibdfile | wc -l)
freqest=$(python lines.py $freqfile 2 2)
python estimate2.py \
    $fileout \
    $ibdest \
    $freqest \
    $cm \
    $n \
    $effdemo \
    $ploidy