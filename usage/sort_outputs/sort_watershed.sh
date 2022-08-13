#!/bin/bash
DIR=/mnt/perkinsdata/tongue_STOmics/benchmarking/25_2k_3/watershed

# orig
mkdir $DIR/orig
mv $DIR/*png $DIR/orig

mkdir -p $DIR/results/outlines $DIR/results/segments $DIR/results/other

for FOLDER in $DIR/nuclei*; do
    SAMPLE=`basename "$FOLDER"`
    SEGMENT_NAME="${SAMPLE/nuclei/segments}"

    mv $DIR/$SAMPLE/4\ segmented.png $DIR/results/outlines/$SAMPLE.png
    mv $DIR/$SAMPLE/segments.csv $DIR/results/segments/$SEGMENT_NAME.csv
    mv $DIR/$SAMPLE $DIR/results/other
done