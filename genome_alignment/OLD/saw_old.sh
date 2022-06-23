#!/bin/bash
ulimit -n 100000

bash SAW/bin/stereoRun_singleLane_v4.0.0.sh \
	-m /mnt/stomics/tongue_stomics/tongue-5/00.fq/tongue-5.h5 \
	-1 /mnt/stomics/tongue_stomics/tongue-5/00.fq/DP8400024518TL_L01_read_1.fq.gz \
	-2 /mnt/stomics/tongue_stomics/tongue-5/00.fq/DP8400024518TL_L01_read_2.fq.gz \
	-g /mnt/stomics/mouse_references/indexed \
	-a /mnt/stomics/mouse_references/gtf/GCA_000001635.9_GRCm39_genomic.gtf \
	-o /mnt/stomics/testing_alignment \ #-i /mnt/stomics/images/tongue-5 \
	-s /data/tongue/SAW_v4.0.0.sif \
	-t 16

