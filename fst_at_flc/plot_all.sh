#!/bin/bash

outdir=$1
echo "Outdir: $outdir"

for csv in `ls $outdir/*.csv` ; do
    echo "======================================="
    echo "$csv"
    echo
    
    ./plot_all.py $csv
    echo
done
