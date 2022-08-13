#!/bin/bash
# Ensure in conda env
DIR=/mnt/stomics/benchmarking/cellpose

# run
python -m cellpose --dir $DIR \
    --pretrained_model nuclei \
    --diameter 0. \
    --save_png