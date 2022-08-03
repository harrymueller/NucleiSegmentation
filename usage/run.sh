#!/bin/bash

if [ $1=cellpose ]; then
    echo "####################"
    echo "# Cellpose"
    echo "####################"

    DIR=/mnt/stomics/benchmarking/cellpose

    # run
    CMD="python -m cellpose --dir $DIR \
        --pretrained_model nuclei \
        --diameter 0. \
        --save_png"

    nohup ./record.sh $CMD > log.out &
    tail -f log.out
elif [ $1=deepcell ]; then
    echo "####################"
    echo "# Deepcell"
    echo "####################"

    export DATA_DIR=/mnt/stomics/benchmarking/deepcell
    export MOUNT_DIR=/data
    export APPLICATION=mesmer

    for f in $DATA_DIR/*.png; do
        export NUCLEAR_FILE=FP200000495BR_E5_2000_1000.png

        CMD="docker run -it \
                -v $DATA_DIR:$MOUNT_DIR \
                vanvalenlab/deepcell-applications:0.3.1 \
                $APPLICATION \
                --nuclear-image $MOUNT_DIR/$NUCLEAR_FILE \
                --output-directory $MOUNT_DIR/02_deepcell \
                --output-name mask.tif \
                --compartment nuclear"
    done

    DIR=/mnt/stomics/benchmarking/cellpose

    # run
    CMD="python -m cellpose --dir $DIR \
        --pretrained_model nuclei \
        --diameter 0. \
        --save_png"

    nohup ./record.sh $CMD > log.out &
    tail -f log.out
fi