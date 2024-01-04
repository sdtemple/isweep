### First input is the directory where the software will be installed

# make folder
mkdir -p $1
cd $1

# get software from the internet
wget https://faculty.washington.edu/browning/beagle/beagle.22Jul22.46e.jar # look for the latest version
wget https://faculty.washington.edu/browning/hap-ibd.jar
wget https://faculty.washington.edu/browning/ibd-ends.jar # build from source if corrupted
wget https://faculty.washington.edu/browning/ibdne/ibdne.23Apr20.ae9.jar
wget https://faculty.washington.edu/browning/beagle_utilities/gtstats.jar

# get auxillary files
wget https://github.com/sdtemple/isweep/blob/main/misc/filter-lines.jar
wget https://github.com/sdtemple/isweep/blob/main/misc/add-uniform-err.jar
wget https://github.com/sdtemple/isweep/blob/main/misc/remove-phase.jar

# renaming some jar files
mv ibdne.23Apr20.ae9.jar ibdne.jar
mv beagle.22Jul22.46e.jar beagle.jar

chmod 755 *.jar

# go back to the original directory
cd ..