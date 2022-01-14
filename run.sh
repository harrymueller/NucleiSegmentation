#!/bin/bash

if [ -z $1 ]; then
    echo "Please pass the name of a file to run"
else
    ID=("tongue-5" "tongue-5" "tongue-5" "tongue-5" "tongue-5" "tongue-5" "tongue-4" "tongue-4" "tongue-4" "tongue-4" "tongue-4")
    NORM_METHOD=("SCT" "SCT" "SCT" "SCT" "LN" "LN" "SCT" "SCT" "SCT" "LN" "LN")
    BINSIZE=("20" "30" "50" "100" "50" "100" "30" "50" "100" "50" "100")
    DIAMETER=("40" "25" "0" "0" "0" "0" "25" "0" "0" "0" "0")
    RESOLUTION=("0.3" "0.3" "0.4" "0.5" "0.5" "0.5" "0.3" "0.4" "0.5" "0.5" "0.5")

    for (( i=0; i<${#ID[@]}; i++ )); do
        Rscript $1 --id ${ID[$i]} \
            --method ${NORM_METHOD[$i]} \
            --binsize ${BINSIZE[$i]} \
            --diameter ${DIAMETER[$i]} \
            --resolution ${RESOLUTION[$i]}
    done

fi