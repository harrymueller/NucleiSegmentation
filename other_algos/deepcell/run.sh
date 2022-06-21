#!/bin/bash
#export DATA_DIR=/mnt/perkinsdata/tongue_STOmics/discovery/deepcell/tongue-4
export DATA_DIR=/mnt/perkinsdata/tongue_STOmics/discovery/nuclei_segregation/gbg/deepcell
export MOUNT_DIR=/data
export APPLICATION=mesmer
#export NUCLEAR_FILE=01_cropping/cropped.png
export NUCLEAR_FILE=FP200000495BR_E5_2000_1000.png

docker run -it \
  -v $DATA_DIR:$MOUNT_DIR \
  vanvalenlab/deepcell-applications:0.3.1 \
  $APPLICATION \
  --nuclear-image $MOUNT_DIR/$NUCLEAR_FILE \
  --output-directory $MOUNT_DIR/02_deepcell \
  --output-name mask.tif \
  --compartment nuclear