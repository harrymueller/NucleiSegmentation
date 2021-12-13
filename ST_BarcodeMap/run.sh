#!/bin/bash
INPUT=/mnt/tongue_data/tongue/original/tongue-4
OUTPUT=/mnt/data/TongueSTOmics/tongue-4

rec ~/tools/ST_BarcodeMap/./ST_BarcodeMap-0.0.1 -V \
    --in $INPUT/00.fq/tongue-4.h5 \
    --in1 $INPUT/00.fq/DP8400023552BR_L01_read_1.fq.gz \
    --in2 $INPUT/00.fq/DP8400023552BR_L01_read_2.fq.gz \
    --out $OUTPUT/combine_read.fq.gz \
    --report $OUTPUT/log4.out \
    --mismatch 0 --thread 2