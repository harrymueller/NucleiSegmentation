#!/bin/bash
nohup docker run -v /data/tongue:/mnt/data tonguediscovery bash $1 >> log.out &
