#!/bin/bash
export DATA_DIR=/mnt/stomics/deepcell/tongue-5
export MOUNT_DIR=/data
export APPLICATION=mesmer
export NUCLEAR_FILE=FP200000495BR_E5.png

docker run -it \
  -v $DATA_DIR:$MOUNT_DIR \
  vanvalenlab/deepcell-applications:0.3.1 \
  $APPLICATION \
  --nuclear-image $MOUNT_DIR/$NUCLEAR_FILE \
  --output-directory $MOUNT_DIR/02_deepcell \
  --output-name mask.tif \
  --compartment nuclear