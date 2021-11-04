#!/bin/bash

#ml cmake gnu7

GRENEDALF=/Carnegie/DPB/Homes/Users/lczech/grenedalf/bin/grenedalf

DATA_ALL="/lustre/scratch/lczech/grenepipe-runs/2021-08-26-ath-grenenet-release-10/mpileup/all-merged-units.mpileup.gz"
DATA_MAPQ60="/lustre/scratch/lczech/grenepipe-runs/2021-08-26-ath-grenenet-release-10-mapq-60/all-merged-units.mpileup.gz"
NAMES="/lustre/scratch/lczech/grenepipe-runs/2021-08-26-ath-grenenet-release-10/mpileup/all-merged-units.names-short.txt"

#$GRENEDALF frequency --pileup-file $DATA_ALL --omit-invariants --sample-name-list ${NAMES} --file-suffix "-all" > grenedalf.log
$GRENEDALF frequency --pileup-file $DATA_MAPQ60 --omit-invariants --sample-name-list ${NAMES} --file-suffix "-mapq60" >> grenedalf.log

