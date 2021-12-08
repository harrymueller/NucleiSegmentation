#!/bin/bash
DIR=`pwd`/logs
# TODO if files exist, delete
# TODO detect earlier stopping and kill processes
rm $DIR/$1_log.out $DIR/$1_usage.txt

# run program via nohup
nohup $@ > $DIR/$1_log.out 2>&1 &
PID=$!

# start pidstat
nohup pidstat -h -r -u -p $PID 1 --human > $DIR/$1_usage.txt 2>&1 &
PIDSTAT=$!

echo "#########################################"
echo "# PID = $PID | PIDSTAT = $PIDSTAT"
echo "#########################################"

# tail logs of program
tail -f -n +0 $DIR/$1_log.out &
wait $PID

# when finished -> analyse data
kill $PIDSTAT # kill pid stat
wait $PIDSTAT 2>/dev/null

echo "#########################################"
echo "# Finished main execution, analysing... "
echo "#########################################"

# print results
python extract_stats.py "$@" $DIR/$1_usage.txt $DIR/$1_stats.txt