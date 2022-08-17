#!/bin/bash
export SINGULARITY_BIND="/mnt/stomics"

imageDir=/mnt/stomics/alignment/images
outputDir=/mnt/stomics/alignment/output
sif=/mnt/stomics/SAW_v4.1.0.sif

sn=FP200000495BR_E5

echo singularity exec ${sif} register \
	-i ${imageDir}/FP200000495BR_E5.tar.gz \
	-c ${imageDir}/FP200000495BR_E5.json \
	-v ${outputDir}/02.count/FP200000495BR_E5.gem \
	-o ${outputDir}/03.register

echo ">>> REGISTER DONE"
singularity exec ${sif} tissueCut \
	--dnbfile ${outputDir}/01.merge/${sn}.barcodeReadsCount.txt \
	-i ${outputDir}/02.count/FP200000495BR_E5.raw.gef \
	-o ${outputDir}/04.tissuecut \
	-s ${outputDir}/03.register/7_result \
	-t tissue \
	--platform T10 \
	--snId ${sn}
