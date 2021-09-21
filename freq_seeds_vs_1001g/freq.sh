#!/bin/bash

#ml cmake gnu7

GRENEDALF=/Carnegie/DPB/Homes/Users/lczech/grenedalf/bin/grenedalf
DATA="/lustre/scratch/lczech/grenepipe-runs/ath-evo-seeds-pileups/mpileup/all-merged-units.mpileup.gz"
NAMES="/lustre/scratch/lczech/grenepipe-runs/ath-evo-seeds-pileups/mpileup/all-merged-units.names.txt"

rm -f frequency.csv

$GRENEDALF frequency --pileup-file $DATA --omit-invariants --type total --sample-name-list ${NAMES} > grenedalf.log

# add a column CHROM_POS that combines both, so that we can easily merge them in R.
cat frequency.csv | sed "s/^\([^\t]*\)\t\([^\t]*\)\t/\1\t\2\t\1_\2\t/g" > frequency-extra.csv

rm frequency.csv
mv frequency-extra.csv frequency.csv
