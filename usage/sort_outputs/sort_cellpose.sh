#!/bin/bash
DIR=/mnt/perkinsdata/tongue_STOmics/benchmarking/cellpose

mkdir $DIR/orig $DIR/pickles $DIR/masks

mv $DIR/*npy $DIR/pickles
mv $DIR/*masks* $DIR/masks
mv $DIR/*png $DIR/orig

py other_alogs/cellpose/convert.py