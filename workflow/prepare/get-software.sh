# make folder
mkdir -p $1
cd $1

# get software from the internet
wget https://faculty.washington.edu/browning/beagle/beagle.22Jul22.46e.jar # look for the latest version
wget https://faculty.washington.edu/browning/hap-ibd.jar
wget https://faculty.washington.edu/browning/flare.jar
wget https://github.com/sdtemple/isweep/blob/main/misc/remove-phase.jar

# renaming the beagle
mv beagle.22Jul22.46e.jar beagle.jar

chmod 755 *.jar

# go back to the original directory
cd ..