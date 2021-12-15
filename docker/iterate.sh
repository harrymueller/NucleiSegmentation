#!/bin/bash
# iterates through tongue ids and binsizes executing the Rscript parsed via $1

SCRIPT_PATH=$1
TONGUE_IDS=("tongue-4" "tongue-5")
#BINS=(1 10 20 50 100)
BINS=(100 50 20 10 1)

for ID in "${TONGUE_IDS[@]}"; do
    for BIN in "${BINS[@]}"; do
        echo "Starting $ID $BIN"
	    Rscript $SCRIPT_PATH \
            --binsize $BIN \
            --id $ID
        echo "Completed $ID $BIN"
    done
done
