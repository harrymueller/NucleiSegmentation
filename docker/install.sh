#!/bin/bash
##############################
# Init setup
##############################
apt-get install git build-essential vim wget apt-utils systemctl -y

# dir to store tools files
mkdir ~/tools && cd ~/tools 

##############################
# 00 ST_BarcodeMap v1.0: https://github.com/BGIResearch/ST_BarcodeMap
##############################
apt-get install libboost-thread-dev zlibc libhdf5-serial-dev -y # dependencies
export CPATH="/usr/include/hdf5/serial/" # add the hdf5 lib to the C compile path
git clone https://github.com/BGIResearch/ST_BarcodeMap

# symlink hdf5 libs
cd /usr/lib/x86_64-linux-gnu # version numbers might change
ln -s libhdf5_serial.so.103.0.0 libhdf5.so
ln -s libhdf5_serial_hl.so.103.0.0 libhdf5_hl.so

# add -lboost_serialization -lhdf5 to libs
cd ~/tools/ST_BarcodeMap
sed -i 's/LIBS := -lz -lpthread/LIBS := -lz -lpthread -lboost_serialization -lhdf5/g' Makefile

make

##############################
# 01.1 FastP: https://github.com/OpenGene/fastp
##############################
cd ~/tools

wget http://opengene.org/fastp/fastp
chmod a+x ./fastp

##############################
# 01.2 STAR
##############################
wget https://github.com/alexdobin/STAR/archive/refs/tags/2.7.2b.tar.gz
tar -xzf 2.7.2b.tar.gz
rm 2.7.2b.tar.gz
cd STAR-2.7.2b/source
make STAR

##############################
# 02? Bam2Gem
##############################
cd ~/tools
source ~/.bashrc

##############################
# 02? Lasso
##############################

##############################
# 03.1 R
##############################
if false; then
    apt-get install gfortran libbz2-dev liblzma-dev libpcre2-dev libcurl4-openssl-dev default-jdk -y

    wget https://cran.r-project.org/src/base/R-4/R-4.1.2.tar.gz
    tar -xzvf R-4.1.2.tar.gz
    rm R-4.1.2.tar.gz

    cd R-4.1.2
    ./configure --with-readline=no --with-x=no --enable-R-shlib; make; make install
fi

##############################
# 03 04 05 07
#   Seurat, monocle3, SingleR, SciBet, clusterProfiler dependencies
##############################
apt-get install libssl-dev libxml2-dev libpng-dev -y

##############################
# 06 CellPhoneDB
##############################
apt-get install python3.8 pipenv python3-venv -y

cd ~/tools
pip install cellphonedb

##############################
# Aliases
##############################
cat aliases >> ~/.bashrc
cat aliases >> /home/rstudio/.bashrc

##############################
# Finished
##############################
echo ""
echo "##############################"
echo "Finished installing all packages."
echo "##############################"
