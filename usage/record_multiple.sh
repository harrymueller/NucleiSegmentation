#!/bin/bash
DIR=/mnt/stomics/benchmarking/deepcell
rm $DIR/log.out $DIR/usage.txt

if [ -z $1 ]; then
    exit 0
fi

CURRENT_DIR=`pwd`

# run program via nohup
cd $DIR
for f in *.png; do
    nohup ./CURRENT_DIR/usage/run_deepcell.sh $DIR $f >> $DIR/log.out 2>&1 &
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


cd CURRENT_DIR
# print results
# python3 /data/tongue/scripts/usage/extract_stats.py "$@" $DIR/usage.txt $DIR/stats.txt