#!/bin/bash
# script to run gem2rds.R multiple times with several bin sizes

SCRIPT_PATH=/mnt/local/scripts/rscripts/gem2rds_cli.R
OUTPUT_DIR=/mnt/data/tongue_STOmics/discovery/gemRDS
INPUT_DIRS=("/mnt/data/tongue_STOmics/original/tongue-4/tongue-4.gem" "/mnt/data/tongue_STOmics/original/tongue-5/tongue-5.gem")
BINS=(1 10 20 50 100)

for INPUT in "${INPUT_DIRS[@]}"; do
    for BIN in "${BINS[@]}"; do
        Rscript $SCRIPT_PATH \
            --infile $INPUT \
            --outdir $OUTPUT_DIR \
            --binsize $BIN \
            --image TRUE
        echo "Completed $BIN $INPUT"
    done
done