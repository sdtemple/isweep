# First step: copy and paste genetic maps

### fixed variables ###

# where
MAPS=$1
MAPPREFIX=$2
MAPSUFFIX=$3

# input variables
WHERE=$4
CHRLOW=$5
CHRHIGH=$6

mkdir -p ${WHERE}/
mkdir -p ${WHERE}/maps/
touch ${WHERE}/maps/chr${CHRLOW}-${CHRLOW}.map
rm ${WHERE}/maps/chr${CHRLOW}-${CHRLOW}.map

# copy and paste maps
for j in $(seq $CHRLOW 1 $CHRHIGH); do cp ${MAPS}/${MAPPREFIX}${j}${MAPSUFFIX} ${WHERE}/maps/chr${j}.map; done;
for j in $(seq $CHRLOW 1 $CHRHIGH); do cat $WHERE/maps/chr${j}.map >> $WHERE/maps/chr${CHRLOW}-${CHRHIGH}.map ; done;
