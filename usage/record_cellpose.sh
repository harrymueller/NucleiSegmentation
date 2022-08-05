#!/bin/bash
DIR=/mnt/stomics/benchmarking/cellpose
rm $DIR/log.out $DIR/usage.txt

if [ -z $1 ]; then
    exit 0
fi

CMD="python -m cellpose --dir $DIR \
        --pretrained_model nuclei \
        --diameter 0. \
        --save_png"

# run program via nohup
nohup $CMD > $DIR/log.out 2>&1 &
PID=$!

# start pidstat
nohup pidstat -h -r -u -p $PID 1 --human > $DIR/usage.txt 2>&1 &
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

echo "#########################################"
echo "# Finished main execution, analysing... "
echo "#########################################"

# print results
python3 /data/tongue/scripts/usage/extract_stats.py "$@" $DIR/usage.txt $DIR/stats.txt