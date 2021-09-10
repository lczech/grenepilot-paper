#!/bin/bash

#ml cmake gnu7

GRENEDALF=/Carnegie/DPB/Homes/Users/lczech/grenedalf/bin/grenedalf
DATA="/lustre/scratch/lczech/grenepipe-runs/2021-08-26-ath-grenenet-release-10/mpileup/all-merged-units.mpileup.gz"
NAMES="/lustre/scratch/lczech/grenepipe-runs/2021-08-26-ath-grenenet-release-10/mpileup/all-merged-units.names-short.txt"

$GRENEDALF frequency --pileup-file $DATA --omit-invariants --sample-name-list ${NAMES} > grenedalf.log

