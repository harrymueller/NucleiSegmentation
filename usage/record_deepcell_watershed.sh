#!/bin/bash
SCRIPT=deepcell
DIR=/mnt/stomics/benchmarking/5x5/${SCRIPT}
NUCLEAR_SEG_EXE=/data/tongue/ssDNA_nuclei_segmentation/NucleiSegregation

rm $DIR/log.out $DIR/usage.txt

# run program via nohup
#if [ $SCRIPT = watershed ]; then
    # start pidstat

 #   nohup pidstat -h -r -u -G NucleiSegregation 1 --human >> $DIR/usage.txt 2>&1 &
  #  PIDSTAT=$!
#fi

for f in $DIR/*.png; do
    filename=`basename "$f"`
    if [ $SCRIPT = deepcell ]; then
        nohup bash usage/run_deepcell.sh $DIR $filename >> $DIR/log.out 2>&1 &
        PID=$!

        # start pidstat
        nohup pidstat -h -r -u -G python 1 --human >> $DIR/usage.txt 2>&1 &
        PIDSTAT=$!

	wait $PID
    elif [ $SCRIPT = watershed ]; then
        nohup $NUCLEAR_SEG_EXE watershed -i $f -o $DIR >> $DIR/log.out 2>&1 &
        PID=$!
	
	nohup pidstat -h -r -u -p $PID 1 --human >> $DIR/usage.txt 2>&1 &
	PIDSTAT=$!

	wait $PID
	mv $DIR/02_watershed "${f/.png}"
    fi

    echo "#########################################"
    echo "# PID = $PID | PIDSTAT = $PIDSTAT"
    echo "#########################################"

    # tail logs of program

    # when finished -> analyse data
    #if [ $SCRIPT != watershed ]; then
        kill $PIDSTAT # kill pid stat
        wait $PIDSTAT 2>/dev/null
    #fi
done 
#if [ $SCRIPT = watershed ]; then
#    kill $PIDSTAT # kill pid stat
#    wait $PIDSTAT 2>/dev/null
#fi
echo "#########################################"
echo "# Finished main execution, analysing... "
echo "#########################################"

# print results
# python3 /data/tongue/scripts/usage/extract_stats.py "$@" $DIR/usage.txt $DIR/stats.txt
