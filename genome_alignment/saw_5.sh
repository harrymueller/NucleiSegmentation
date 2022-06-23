#!/bin/bash
ulimit -n 100000
DIR='/mnt/stomics/alignment'
SN_ID="FP200000495BR_E5"
export SINGULARITY_BIND="$DIR,/mnt/stomics/mouse_references/STAR_SJ100" #,/mnt/stomics/mouse_references/STAR_SJ100:$DIR/reference/STAR_SJ100"
bash stereoRun_singleLane_v4.1.0.sh \
	-m $DIR/mask/$SN_ID.h5 \
	-1 $DIR/reads/${SN_ID}_read_1.fq.gz \
	-2 $DIR/reads/${SN_ID}_read_2.fq.gz \
	-g /mnt/stomics/mouse_references/STAR_SJ100 \ #$DIR/reference/STAR_SJ100 \
	-a $DIR/reference/genes.gtf \
	-o $DIR/output \ #-i $DIR/images \
	-s /mnt/stomics/SAW_v4.1.0.sif \
	-t 16 \
	-c 2.7

