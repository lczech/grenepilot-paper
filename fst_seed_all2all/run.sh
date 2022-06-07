#!/bin/bash

# Need modules. Load from outside the script though.
# ml cmake gnu7

# calc cluster
#GRENEDALF=/Carnegie/DPB/Homes/Users/lczech/grenedalf/bin/grenedalf
#DATA="/lustre/scratch/lczech/grenepipe-runs/ath-evo-seeds-pileups/mpileup/all-merged-units.mpileup.gz"

# caltech cluster
GRENEDALF=/central/home/lczech/grenedalf/bin/grenedalf
DATA="/central/groups/carnegie_poc/lczech/lustre_scratch/grenepipe-runs/ath-evo-seeds-pileups/mpileup/all-merged-units.mpileup.gz"

# test
#$GRENEDALF fst --pileup-file $DATA --window-width 1 --omit-na-windows --pool-sizes 2500 --sample-name-prefix S --file-suffix "-width-1" > grenedalf-fst-width-1.log

# Old run with grenedalf before unbiased fst was implemented. This hence uses Kofler instead.
#$GRENEDALF fst --pileup-file $DATA --window-width 10000 --omit-na-windows --pool-sizes 2500 --sample-name-prefix S --file-suffix "-width-10000" > grenedalf-fst-width-10000.log

# Final run, with new unbiased Nei estimator. We use the "old" grenedalf syntax here with spence-nei as the method, which will soon change to "unbiased-nei" or something similar. But at the time of writing this script, this has not happened yet.
$GRENEDALF fst --pileup-path $DATA --window-width 10000 --omit-na-windows --pool-sizes 2500 --sample-name-prefix S --file-suffix "-width-10000" --method spence-nei --threads 4 > grenedalf-fst-width-10000.log
