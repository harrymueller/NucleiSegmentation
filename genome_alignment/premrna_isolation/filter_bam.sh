#!/bin/bash
# files
INPUT=/mnt/stomics/alignment/output/00.mapping/FP200000495BR_E5_read_1.Aligned.sortedByCoord.out.bam
BED=/mnt/stomics/alignment/premrna_ref/ucsc_introns_no_exons.bed
OUTPUT=/mnt/stomics/alignment/output/introns_no_exons_u_50.bam

# intersect
bedtools intersect -a $INPUT -b $BED \
	-wa -s -u -f 0.5 -split > $OUTPUT

