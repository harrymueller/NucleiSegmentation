#!/bin/bash
if [ -z $2 ]; then
	echo "-p pawsey"
	echo "-l laptop"
elif [ $2 == -p ]; then
VOLUMES="-v /mnt/stomics:/mnt/data \
    -v /data/tongue/scripts:/mnt/data/scripts \
    -v /data/tongue/rstudio:/mnt/data/rstudio"
elif [ $2 == -l ]; then 
	VOLUMES="-v /mnt/perkinsdata/tongue_STOmics/discovery:/mnt/data \
    -v /data/Perkins/TongueSTOmics/scripts:/mnt/data/scripts \
    -v /data/Perkins/TongueSTOmics/rstudio:/mnt/data/rstudio"
fi

if [ -z $1 ]; then
    opt="-h"
else
    opt=$1
fi

if [ $opt == -t ]; then
    docker run --rm $VOLUMES \
        -p 8787:8787 \
        -ti tonguediscovery /bin/bash
elif [ $opt == -s ]; then
    docker run --rm $VOLUMES \
        -ti stomics/saw:02.1.0 /bin/bash
else
    echo "-t    Tongue Discovery Container"
    echo "-s    SAW container"
fi
