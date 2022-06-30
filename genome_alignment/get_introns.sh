#!/bin/bash
dir=/mnt/stomics/alignment/output
bedtools intersect -a $dir/00.mapping/FP200000495BR_E5_read_1.Aligned.sortedByCoord.out.bam \
	-b /mnt/stomics/alignment/premrna_ref/ucsc_introns_no_exons.bed \
	-wa -s -u -f 0.5 -split > $dir/introns_no_exons_u_50.bam

