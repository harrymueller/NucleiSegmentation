#!/bin/bash
DIR=/mnt/perkinsdata/tongue_STOmics/benchmarking/25_2k_4/cellpose

mkdir $DIR/orig $DIR/pickles $DIR/masks

mv $DIR/*npy $DIR/pickles
mv $DIR/*masks* $DIR/masks
mv $DIR/*png $DIR/orig

python3 other_algos/cellpose/convert.py