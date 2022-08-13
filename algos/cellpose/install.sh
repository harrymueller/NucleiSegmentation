#!/bin/bash

# create and activate conda env
conda create -n cellpose python=3.8
conda activate cellpose

# install packages
python -m pip install 'cellpose[all]'

# for GPUs
# pip uninstall torch
# conda install pytorch cudatoolkit=10.2 -c pytorch