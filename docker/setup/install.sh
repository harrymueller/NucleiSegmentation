#!/bin/bash
##############################
# Init setup
##############################
apt-get install git build-essential vim wget apt-utils systemctl tmux htop sysstat -y

# dir to store tools files
mkdir ~/tools && cd ~/tools 

##############################
# Aliases
##############################
cat aliases >> /root/.bashrc
cat aliases >> /home/rstudio/.bashrc
source ~/.bashrc

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
# 03.1 Rstudio Working Directory config
##############################
RS_PATH=/mnt/local/rstudio
echo "session-default-working-dir=$RS_PATH" >> /etc/rstudio/rsession.conf
echo "session-default-new-project-dir=$RS_PATH" >> /etc/rstudio/rsession.conf

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
pip install CellPhoneDB

##############################
# x1 Singularity
##############################
apt-get install uuid-dev libgpgme-dev squashfs-tools libseccomp-dev \
    pkg-config cryptsetup-bin -y

# install GO
export VERSION=1.14.12 OS=linux ARCH=amd64 && \
    wget https://dl.google.com/go/go$VERSION.$OS-$ARCH.tar.gz && \
    sudo tar -C /usr/local -xzvf go$VERSION.$OS-$ARCH.tar.gz && \
    rm go$VERSION.$OS-$ARCH.tar.gz
echo 'export GOPATH=${HOME}/go' >> ~/.bashrc && \
    echo 'export PATH=/usr/local/go/bin:${PATH}:${GOPATH}/bin' >> ~/.bashrc && \
    source ~/.bashrc

wget https://github.com/sylabs/singularity/releases/download/v3.9.1/singularity-ce_3.9.1+6-g38b50cbc5-focal_amd64.deb
dpkg -i singularity-ce_3.9.1+6-g38b50cbc5-focal_amd64.deb

##############################
# x2 SAW
##############################
singularity build SAW_v1.0.0.sif docker://stomics/saw:01.0.0 
cd ~/tools
git clone https://github.com/BGIResearch/SAW.git

##############################
# others
##############################
# R raster
apt-get install libudunits2-dev libgdal-dev libgeos-dev libproj-dev libxt-dev -y

# imagemagick
apt-get install imagemagick -y

##############################
# Finished
##############################
echo ""
echo "##############################"
echo "Finished installing all packages."
echo "##############################"


export LD_LIBRARY_PATH=/lib:/usr/lib:/usr/local/lib