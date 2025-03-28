### First input is the directory where the software will be installed

# make folder
mkdir -p software
cd software

# get software from the internet
wget https://faculty.washington.edu/browning/beagle/beagle.17Dec24.224.jar # look for the latest version
wget https://faculty.washington.edu/browning/hap-ibd.jar
wget https://faculty.washington.edu/browning/ibd-ends.jar # build from source if corrupted
wget https://faculty.washington.edu/browning/ibdne/ibdne.23Apr20.ae9.jar
wget https://faculty.washington.edu/browning/flare.jar

# renaming some jar filesd .
mv ibdne.23Apr20.ae9.jar ibdne.jar
mv beagle.17Dec24.224.jar beagle.jar

chmod 755 *.jar

# go back to the original directory
cd ..
