#!/bin/bash
DIR=/mnt/perkinsdata/tongue_STOmics/benchmarking/25_2k_4/deepcell

mkdir -p $DIR/orig $DIR/masks $DIR/results/segments $DIR/results/outlines
mv $DIR/*png $DIR/orig
mv $DIR/*.tif $DIR/masks

python3 algos/deepcell/convert.py