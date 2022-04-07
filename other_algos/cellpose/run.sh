#!/bin/bash

python -m cellpose --dir /mnt/perkinsdata/tongue_STOmics/discovery/cellpose/images \
 --pretrained_model nuclei --diameter 0. --save_png
nohup python -m cellpose --dir /data/tongue/cellpose/images --pretrained_model nuclei --diameter 0. --save_png &