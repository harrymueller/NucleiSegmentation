#!/bin/bash
DIR="/mnt/perkinsdata/tongue_STOmics/benchmarking/25_2k_4"

echo "Watershed..."
python3 get_A_B.py $DIR watershed

echo "Cellpose..."
python3 get_A_B.py $DIR cellpose

echo "Deepcell..."
python3 get_A_B.py $DIR deepcell

echo "Measures..."
python3 get_measures.py $DIR

echo "Summarise..."
python3 summarise.py $DIR

echo "Format Summary..."
python3 format_summary.py $DIR