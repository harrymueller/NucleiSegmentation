#!/bin/bash
ulimit -n 100000
DIR='/mnt/stomics/alignment'

bash SAW_OLD/bin/stereoRun_singleLane_v4.0.0.sh \
	-m $DIR/mask/tongue-5.h5 \
	-1 $DIR/reads/tongue-5_read_1.fq.gz \
	-2 $DIR/reads/tongue-5_read_2.fq.gz \
	-g $DIR/reference/STAR_SJ100 \
	-a $DIR/reference/genes.gtf \
	-o $DIR/output \
	-i $DIR/images \
	-s /data/tongue/SAW_v4.0.0.sif \
	-t 16

