#!/bin/bash
SCRIPT=watershed
DIR=/mnt/stomics/benchmarking/${SCRIPT}
NUCLEAR_SEG_EXE=/data/tongue/ssDNA_nuclei_segmentation/NucleiSegregation

rm $DIR/log.out $DIR/usage.txt

# run program via nohup
for f in $DIR/*.png; do
    filename=`basename "$f"`
    if [ $SCRIPT = deepcell ]; then
        nohup bash usage/run_deepcell.sh $DIR $filename >> $DIR/log.out 2>&1 &
    elif [ $SCRIPT = watershed ]; then
        newdir=$DIR/"${filename/.png/}"
        mkdir newdir
        nohup bash $NUCLEAR_SEG_EXE -i $DIR/$filename -o newdir >> $DIR/log.out 2>&1 &
    fi
    PID=$!

    # start pidstat
    nohup pidstat -h -r -u -p $PID 1 --human >> $DIR/usage.txt 2>&1 &
    PIDSTAT=$!

    echo "#########################################"
    echo "# PID = $PID | PIDSTAT = $PIDSTAT"
    echo "#########################################"

    # tail logs of program
    tail -f -n +0 $DIR/log.out &
    wait $PID

    # when finished -> analyse data
    kill $PIDSTAT # kill pid stat
    wait $PIDSTAT 2>/dev/null
done 

echo "#########################################"
echo "# Finished main execution, analysing... "
echo "#########################################"

# print results
# python3 /data/tongue/scripts/usage/extract_stats.py "$@" $DIR/usage.txt $DIR/stats.txt