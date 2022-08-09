#!/bin/bash
DIR=/mnt/perkinsdata/tongue_STOmics/benchmarking/25_2k_2

# orig
mkdir $DIR/orig
mv $DIR/*png $DIR/orig

mkdir -p $DIR/results/outlines $DIR/results/segments $DIR/results/other

for FOLDER in $DIR/nuclei*; do
    SAMPLE=`basename "$FOLDER"`
    echo $FOLDER

    mv $DIR/$SAMPLE/4\ segmented.png $DIR/results/outlines/$SAMPLE.png
    mv $DIR/$SAMPLE/segments.csv $DIR/results/segments/$SAMPLE.csv
    mv $DIR/$SAMPLE $DIR/results/other
done