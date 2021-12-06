#!/bin/bash
##############################
# Init setup
##############################
apt-get install git build-essential vim wget -y # git
echo "alias l='ls -lah'" >> ~/.bashrc
echo "alias py='python3'" >> ~/.bashrc

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

# add cli command 
echo "alias ST_BarcodeMap='~/tools/ST_BarcodeMap/./ST_BarcodeMap-0.0.1'" >> ~/.bashrc

##############################
# 01.1 FastP: https://github.com/OpenGene/fastp
##############################
cd ~/tools

wget http://opengene.org/fastp/fastp
chmod a+x ./fastp

echo "alias fastp='~/tools/./fastp'" >> ~/.bashrc

##############################
# 01.2 STAR
##############################
wget https://github.com/alexdobin/STAR/archive/refs/tags/2.7.2b.tar.gz
tar -xzf 2.7.2b.tar.gz
rm 2.7.2b.tar.gz
cd STAR-2.7.2b/source
make STAR

echo "alias STAR='~/tools/STAR-2.7.2b/source/./STAR'" >> ~/.bashrc

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
apt-get install gfortran libbz2-dev liblzma-dev libpcre2-dev libcurl4-openssl-dev default-jdk -y

wget https://cran.r-project.org/src/base/R-4/R-4.1.2.tar.gz
tar -xzvf R-4.1.2.tar.gz
rm R-4.1.2.tar.gz

cd R-4.1.2
./configure --with-readline=no --with-x=no --enable-R-shlib; make; make install

##############################
# 03 04 05 07
#   Seurat, monocle3, SingleR, SciBet, clusterProfiler
##############################
apt-get install libssl-dev libxml2-dev libpng-dev -y
Rscript r_pkgs.R

##############################
# 06 CellPhoneDB
##############################
apt-get install python3.8 pipenv python3-venv -y

cd ~/tools
python -m venv cpdb # docs recommend installing in venv
source cpdb/bin/activate
pip install cellphonedb
deactivate
##############################
# rStudio Server
##############################
apt-get install gdebi-core -y
wget https://download2.rstudio.org/server/bionic/amd64/rstudio-server-2021.09.1-372-amd64.deb
gdebi rstudio-server-2021.09.1-372-amd64.deb -y

##############################
# Finished
##############################
echo ""
echo "##############################"
echo "Finished installing all packages."
echo "##############################"