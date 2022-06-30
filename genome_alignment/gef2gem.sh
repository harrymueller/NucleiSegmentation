#!/bin/bash
export SINGULARITY_BIND="/mnt/stomics"
DIR="/mnt/stomics/alignment/output/02.count_introns_no_exons_u_50"
singularity exec /mnt/stomics/SAW_v4.1.0.sif gefTools view \
	-i $DIR/FP200000495BR_E5.raw.gef \
	-o $DIR/FP200000495BR_E5.gem \
	-b 1
