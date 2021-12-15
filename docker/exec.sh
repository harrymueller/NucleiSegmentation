#!/bin/bash
# iterates through tongue ids and binsizes executing the Rscript parsed via $1 in disconnected docker container
# $1 relative to rscripts dir
nohup docker run -v /data/tongue:/mnt/data \
    tonguediscovery \
    bash /mnt/data/scripts/docker/iterate.sh /mnt/data/scripts/rscripts/$1 \
    > log.out && tail -f log.out