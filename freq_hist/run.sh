#!/bin/bash

#ml cmake gnu7

GRENEDALF=/Carnegie/DPB/Homes/Users/lczech/grenedalf/bin/grenedalf
DATA="/lustre/scratch/lczech/grenepipe-runs/ath-evo-francois-rerun/genotyped/all.vcf.gz"

$GRENEDALF frequency --vcf-file $DATA --omit-invariants > grenedalf.log

mv frequency.csv frequency-francois-vcf.csv
