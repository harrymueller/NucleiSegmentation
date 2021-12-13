#!/bin/bash
DIR=`pwd`/st_barcodemap_4
# TODO if files exist, delete
# TODO detect earlier stopping and kill processes
rm $DIR/log.out $DIR/usage.txt

# run program via nohup
nohup $@ > $DIR/log.out 2>&1 &
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
python3 /mnt/local/TongueSTOmics/usage/extract_stats.py "$@" $DIR/usage.txt $DIR/stats.txt