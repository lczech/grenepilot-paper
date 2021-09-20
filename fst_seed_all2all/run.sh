#!/bin/bash

# Need modules. Load from outside the script though.
# ml cmake gnu7

GRENEDALF=/Carnegie/DPB/Homes/Users/lczech/grenedalf/bin/grenedalf
DATA="/lustre/scratch/lczech/grenepipe-runs/ath-evo-seeds-pileups/mpileup/all-merged-units.mpileup.gz"

$GRENEDALF fst --pileup-file $DATA --window-width 1 --omit-na-windows --pool-sizes 2500 --sample-name-prefix S > grenedalf-fst.log
