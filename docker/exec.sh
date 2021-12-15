#!/bin/bash
# run bash script ($1) in disconnected docker container
# $1 relative to top dir of git repo
nohup docker run -rm -v /data/tongue:/mnt/data tonguediscovery bash /mnt/data/scripts/$1 >> log.out &
