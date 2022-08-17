#!/bin/bash
export SINGULARITY_BIND="/mnt/stomics"
name=introns_no_exons_u_50
mkdir -p /mnt/stomics/alignment/output/02.count_${name}
singularity exec /mnt/stomics/SAW_v4.1.0.sif count \
   -i /mnt/stomics/alignment/output/${name}.bam \
   -o /mnt/stomics/alignment/output/02.count_${name}/FP200000495BR_E5.Aligned.sortedByCoord.out.merge.q10.dedup.target.bam \
   -a /mnt/stomics/mouse_references/gencode/gencode.vM29.primary_assembly.annotation.gtf \
   -s /mnt/stomics/alignment/output/02.count_${name}/FP200000495BR_E5.Aligned.sortedByCoord.out.merge.q10.dedup.target.bam.summary.stat \
   -e /mnt/stomics/alignment/output/02.count_${name}/FP200000495BR_E5.raw.gef \
   --sat_file /mnt/stomics/alignment/output/02.count_${name}/FP200000495BR_E5_raw_barcode_gene_exp.txt \
   --umi_on \
   --sn FP200000495BR_E5 \
   -c 32 \
   -m 128
