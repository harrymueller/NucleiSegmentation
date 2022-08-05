#!/bin/bash
DIR=/mnt/perkinsdata/tongue_STOmics/benchmarking/deepcell

mkdir -p $DIR/orig $DIR/results/segments $DIR/results/outlines
mv $DIR/*png $DIR/orig

py other_alogs/deepcell/convert.py