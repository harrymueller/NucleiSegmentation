#!/bin/bash

DATA_DIR='/mnt/stomics/mouse_references'
export SINGULARITY_BIND=$DATA_DIR

singularity exec /mnt/stomics/SAW_v4.1.0.sif mapping \
	--runMode genomeGenerate \
	--genomeDir $DATA_DIR/STAR_SJ100_GENCODE \
	--genomeFastaFiles $DATA_DIR/gencode/GRCm39.primary_assembly.genome.fa \
	--sjdbGTFfile $DATA_DIR/gencode/gencode.vM29.primary_assembly.annotation.gtf \
	--sjdbOverhang 99 \
	--runThreadN 16
