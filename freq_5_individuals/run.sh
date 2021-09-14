#!/bin/bash

ml cmake gnu7

GRENEDALF=/Carnegie/DPB/Homes/Users/lczech/grenedalf/bin/grenedalf
DATA="/lustre/scratch/lczech/grenepipe-runs/ath-evo-flowerpools-rerun/mpileup/all-merged-units.mpileup.gz"

#        { "S1",  "Flowerpool1001und2" },
#        { "S2",  "Flowerpool50A" },
#        { "S3",  "Flowerpool50B" },
#        { "S4",  "Flowerpool25A" },
#        { "S5",  "Flowerpool25B" },
#        { "S6",  "Flowerpool10B2" },
#        { "S7",  "Flowerpool5B" },
#        { "S8",  "Plantpool100x" },
#        { "S9",  "Plantpool50a" },
#        { "S10", "Plantpool50b" },
#        { "S11", "Plantpool25a" },
#        { "S12", "Plantpool25b" },
#        { "S13", "Plantpool10b1" },
#        { "S14", "Plantpool10b2" },
#        { "S15", "Plantpool5x" },
#        { "S16", "PlantpoolB12345" }

$GRENEDALF frequency --pileup-file $DATA --omit-invariants --sample-name-prefix S --filter-samples-include "S7,S15" > grenedalf.log

#$GRENEDALF frequency --pileup-file $DATA --omit-invariants --sample-name-prefix S > grenedalf.log


