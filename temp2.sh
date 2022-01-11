#!/bin/bash
for id in "tongue-4" "tongue-5"; do
	echo $id
	for method in "SCT" "LN"; do
		Rscript rscripts/cli/dim_reduction.R --binsize 10 --id $id --diameter 60 --method $method
		Rscript rscripts/cli/dim_reduction.R --binsize 20 --id $id --diameter 40 --method $method
		Rscript rscripts/cli/dim_reduction.R --binsize 30 --id $id --diameter 25 --method $method
	done
done
