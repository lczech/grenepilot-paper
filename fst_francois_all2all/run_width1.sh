#!/bin/bash

# Need modules. Load from outside the script though.
# ml cmake gnu7

# memex cluster
#GRENEDALF=/Carnegie/DPB/Homes/Users/lczech/grenedalf/bin/grenedalf
#DATA="/lustre/scratch/lczech/grenepipe-runs/ath-evo-francois-rerun/mpileup/all-merged-units.mpileup.gz"
#NAMES="/lustre/scratch/lczech/grenepipe-runs/ath-evo-francois-rerun/mpileup/all-merged-units.names.txt"
#POOLS="/home/lczech/safedata/ath_evo/grenepilot_lucas/fst_francois_all2all/poolsizes.txt"

# calc cluster
GRENEDALF="/home/lczech/grenedalf/bin/grenedalf"
DATA="/Carnegie/DPB/Data/Shared/Labs/Moi/Everyone/ath_evo/grenepilot_lucas_quick_for_paper_submission/ath-evo-francois-rerun/mpileup/all-merged-units.mpileup.gz"
NAMES="/Carnegie/DPB/Data/Shared/Labs/Moi/Everyone/ath_evo/grenepilot_lucas_quick_for_paper_submission/ath-evo-francois-rerun/mpileup/all-merged-units.names.txt"
POOLS="/home/lczech/safedata/ath_evo/grenepilot_lucas/fst_francois_all2all/poolsizes.txt"

# caltech cluster
#GRENEDALF=/central/home/lczech/grenedalf/bin/grenedalf
#DATA="/central/groups/carnegie_poc/lczech/lustre_scratch/grenepipe-runs/ath-evo-francois-rerun/mpileup/all-merged-units.mpileup.gz"
#NAMES="/central/groups/carnegie_poc/lczech/lustre_scratch/grenepipe-runs/ath-evo-francois-rerun/mpileup/all-merged-units.names.txt"
#POOLS="/central/groups/carnegie_poc/lczech/grenepilot-paper/fst_francois_all2all/poolsizes.txt"

# pool sizes (I think... this is the "number of flowers" field from the overview table...)
#S1  1_0 56  
#S4  1_1 80  
#S7  1_2 65  
#S10 1_3 19  
#S2  2_0 69  
#S5  2_1 160 
#S8  2_2 205 
#S11 2_3 50  
#S3  3_0 101 
#S6  3_1 200 
#S9  3_2 296 
#S12 3_3 97  

# Tests
#$GRENEDALF fst --pileup-file $DATA --window-width 1 --omit-na-windows --pool-sizes "56,80,65,19,69,160,205,50,101,200,296,97" --sample-name-prefix S --file-suffix "-width-1" > grenedalf-fst-width-1.log
#$GRENEDALF fst --pileup-file $DATA --window-width 10000 --omit-na-windows --pool-sizes "56,80,65,19,69,160,205,50,101,200,296,97" --sample-name-list ${NAMES} --file-suffix "-width-10k" > grenedalf-fst-width-10k.log

# Old run with grenedalf before unbiased fst was implemented. This hence uses Kofler instead.
#$GRENEDALF fst --pileup-file $DATA --window-width 10000 --omit-na-windows --pool-sizes ${POOLS} --sample-name-list ${NAMES} --file-suffix "-width-10k" > grenedalf-fst-width-10k.log

# Final run, with new unbiased Nei estimator. We use the "old" grenedalf syntax here with spence-nei as the method, which will soon change to "unbiased-nei" or something similar. But at the time of writing this script, this has not happened yet.
$GRENEDALF fst --pileup-path $DATA --window-width 1 --method spence-nei --omit-na-windows --pool-sizes ${POOLS} --sample-name-list ${NAMES} --file-suffix "-width-1" --threads 4 > grenedalf-fst-width-1.log
