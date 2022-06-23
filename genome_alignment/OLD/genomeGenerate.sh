#!/bin/bash

DATA_DIR='/mnt/stomics/mouse_references'
export SINGULARITY_BIND=$DATA_DIR

singularity exec /mnt/stomics/SAW_v4.1.0.sif \ #/data/tongue/SAW_v4.0.0.sif mapping \
	--runMode genomeGenerate \
	--genomeDir $DATA_DIR/STAR_SJ100_410 \
	--genomeFastaFiles $DATA_DIR/gencode/GRCm39.primary_assembly.genome.fa \ #$DATA_DIR/genome_assemblies_fna/GCF_000001635.27_GRCm39_genomic.fna \
	--sjdbGTFfile $DATA_DIR/gencode/gencode.vM29.primary_assembly.annotation.gtf \ #genome_assemblies_gtf/GCF_000001635.27_GRCm39_genomic.gtf \ #$DATA_DIR/gtf/GCA_000001635.9_GRCm39_genomic.gtf \
	--sjdbOverhang 99 \
	--runThreadN 16
