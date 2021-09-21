#!/bin/bash

#ml cmake gnu7

GRENEDALF=/Carnegie/DPB/Homes/Users/lczech/grenedalf/bin/grenedalf
DATA="/lustre/scratch/lczech/grenepipe-runs/ath-evo-francois-rerun/mpileup/all-merged-units.mpileup.gz"
NAMES="/lustre/scratch/lczech/grenepipe-runs/ath-evo-francois-rerun/mpileup/all-merged-units.names.txt"

$GRENEDALF frequency --pileup-file $DATA --omit-invariants --sample-name-list ${NAMES} > grenedalf.log

