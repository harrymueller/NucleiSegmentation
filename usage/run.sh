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

    nohup ./record $CMD > log.out &
fi