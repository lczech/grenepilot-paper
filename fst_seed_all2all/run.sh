#!/bin/bash

# Need modules. Load from outside the script though.
# ml cmake gnu7

GRENEDALF=/Carnegie/DPB/Homes/Users/lczech/grenedalf/bin/grenedalf
DATA="/lustre/scratch/lczech/grenepipe-runs/ath-evo-seeds-pileups/mpileup/all-merged-units.mpileup.gz"

#$GRENEDALF fst --pileup-file $DATA --window-width 1 --omit-na-windows --pool-sizes 2500 --sample-name-prefix S --file-suffix "-width-1" > grenedalf-fst-width-1.log

$GRENEDALF fst --pileup-file $DATA --window-width 10000 --omit-na-windows --pool-sizes 2500 --sample-name-prefix S --file-suffix "-width-10000" > grenedalf-fst-width-10000.log
