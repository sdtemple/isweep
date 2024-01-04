# run the IBDNe program on a set of chromosomes

# inputs

ibdne=$1 # jar file
xmx=$2 # memory in gb
study=$3 
ibd=$4 # study/ibd make directory for chr.ibd.gz files
chrlow=$5 # lowest chromosome number
chrhigh=$6 # highest chromosome number
outfile=$7 # prefix for ibdne files
prefix=${study}/${ibd}

# merge map files
for j in $(seq ${chrlow} 1 ${chrhigh}); do cat ${study}/maps/chr${j}.map >> ${study}/maps/chr${chrlow}-${chrhigh}.map ; done;

# merge ibd files
for j in $(seq ${chrlow} 1 ${chrhigh}); do zcat ${prefix}/chr${j}.ibd.gz >> ${prefix}/chrall.ibd ; done;

# run ibdne
cat ${prefix}/chrall.ibd | java -Xmx${xmx}g -jar ${ibdne} map=${study}/maps/chr${chrlow}-${chrhigh}.map out=${study}/${outfile}

# remove temp files 
rm ${prefix}/chrall.ibd
rm ${study}/maps/chr${chrlow}-${chrhigh}.map
