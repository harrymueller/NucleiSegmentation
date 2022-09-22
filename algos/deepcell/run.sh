#!/bin/bash
#export DATA_DIR=/mnt/stomics/deepcell/tongue-5
export DATA_DIR=/mnt/stomics/other_stains/jen_c/deepcell
export MOUNT_DIR=/data
export APPLICATION=mesmer
#export NUCLEAR_FILE=FP200000495BR_E5.png
export NUCLEAR_FILE=Sydney_1_A6_v2.png

docker run -it \
  -v $DATA_DIR:$MOUNT_DIR \
  vanvalenlab/deepcell-applications:0.3.1 \
  $APPLICATION \
  --nuclear-image $MOUNT_DIR/$NUCLEAR_FILE \
  --output-directory $MOUNT_DIR/02_deepcell \
  --output-name mask.tif \
  --compartment nuclear
