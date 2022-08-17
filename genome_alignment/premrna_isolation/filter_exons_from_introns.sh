#!/bin/bash
DIR=/mnt/stomics/alignment/premrna_ref

bedtools subtract -a $DIR/ucsc_introns.bed 
                  -b $DIR/ucsc_exons.bed -A 
                  > $DIR/ucsc_introns_no_exons.bed
