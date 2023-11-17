ibdne=$1
xmx=$2
folder=$3
prefix=$4
chrlow=$5
chrhigh=$6
outfile=$7

for j in $(seq ${chrlow} 1 ${chrhigh}); do cat ${folder}/maps/chr${j}.map >> ${folder}/maps/chr${chrlow}-${chrhigh}.map ; done;
for j in $(seq ${chrlow} 1 ${chrhigh}); do zcat ${prefix}chr${j}.ibd.gz >> ${prefix}chrall.ibd ; done;
cat ${prefix}chrall.ibd | java -Xmx${xmx}g -jar ${ibdne} map=${folder}/maps/chr${chrlow}-${chrhigh}.map out=${outfile} 
rm ${prefix}chrall.ibd
rm ${folder}/maps/chr${chrlow}-${chrhigh}.map