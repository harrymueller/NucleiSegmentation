#!/bin/bash
DIR=/mnt/perkinsdata/tongue_STOmics/benchmarking/25_2k_4_thresholded/cellpose

mkdir $DIR/orig $DIR/pickles $DIR/masks

mv $DIR/*npy $DIR/pickles
mv $DIR/*masks* $DIR/masks
mv $DIR/*png $DIR/orig

python3 algos/cellpose/convert.py $DIR
