#!/bin/bash
for id in "tongue-4" "tongue-5"; do
	echo $id
	Rscript rscripts/cli/spatial_plotting.R --binsize 10 --id $id --diameter 60
	Rscript rscripts/cli/spatial_plotting.R --binsize 20 --id $id --diameter 40
	Rscript rscripts/cli/spatial_plotting.R --binsize 30 --id $id --diameter 25
done
