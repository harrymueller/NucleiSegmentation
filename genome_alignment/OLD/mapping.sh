#!/bin/bash

export SINGULARITY_BIND="/mnt/stomics"
uname -n 100000

DATA_DIR='/mnt/stomics/mouse_references'

singularity exec /data/tongue/SAW_v4.0.0.sif mapping \
	--outSAMattributes spatial \
	--outSAMtype BAM SortedByCoordinate \
	--genomeDir $DATA_DIR/STAR_SJ100 \
	--runThreadN 16 \
	--outFileNamePrefix /data/02.alignment/tongue5. \
	--sysShell /bin/bash \
	--stParaFile /data/tongue5.bcPara \
	--readNameSeparator \" \" \
	--limitBAMsortRAM 38582880124 \
	--limitOutSJcollapsed 10000000 \
	--limitIObufferSize=280000000 \
	> /mnt/stomics/testing_alignment_manual/00.fq/tongue5_barcodeMap.stat
