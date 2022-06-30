#!/bin/bash
#bedtools intersect -a outputGenCode/introns_50.bam -b premrna_ref/ucsc_exons.bed -f 0.05 -wa -split > outputGenCode/introns_50_exons.bam
#bedtools intersect -a /mnt/stomics/alignment/outputGenCode/00.mapping/FP200000495BR_E5_read_1.Aligned.sortedByCoord.out.bam -b /mnt/stomics/alignment/premrna_ref/ucsc_exons.bed -f 0.5 -wa -split > /mnt/stomics/alignment/outputGenCode/axolotl/exons.bam
bedtools subtract -a /mnt/stomics/alignment/outputGenCode/axolotl/introns_s_exons.bam -b /mnt/stomics/alignment/premrna_ref/ucsc_exons.bed -f 0.5 -A -s > /mnt/stomics/alignment/outputGenCode/axolotl/introns_s_exons_s.bam
#bedtools intersect -a /mnt/stomics/alignment/outputGenCode/axolotl/introns_s.bam -b /mnt/stomics/alignment/premrna_ref/ucsc_exons.bed -wa -split > /mnt/stomics/alignment/outputGenCode/axolotl/introns_s_exons.bam
