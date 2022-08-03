#!/bin/bash
DATA_DIR=$1
MOUNT_DIR=/data
APPLICATION=mesmer
NUCLEAR_FILE=$2
MASK_FILE="${NUCLEAR_FILE/.png/_mask.tif}"

docker run \
  -v $DATA_DIR:$MOUNT_DIR \
  vanvalenlab/deepcell-applications:0.3.1 \
  $APPLICATION \
  --nuclear-image $MOUNT_DIR/$NUCLEAR_FILE \
  --output-directory $MOUNT_DIR/02_deepcell \
  --output-name $MASK_FILE \
  --compartment nuclear