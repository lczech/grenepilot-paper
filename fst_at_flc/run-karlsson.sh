#!/bin/bash
#SBATCH --time=1-0
#SBATCH --mem=5G
#SBATCH --cpus-per-task=8
#SBATCH --nodes=1
#SBATCH --ntasks=1

echo "Running karlsson"

# debugging
module load cmake/3.19.4 gcc/9.2.0	
ulimit -c unlimited

GRENEDALF=/home/lczech/grenedalf/bin/grenedalf

MPILEUP=/central/groups/carnegie_poc/lczech/lustre_scratch/grenepipe-runs/ath-evo-francois-rerun/mpileup/all-merged-units.mpileup.gz
NAMES=/central/groups/carnegie_poc/lczech/lustre_scratch/grenepipe-runs/ath-evo-francois-rerun/mpileup/all-merged-units.names.txt

# Poolsizes S1..S12
# 56, 69, 101, 80, 160, 200, 65, 205, 296, 19, 50, 97

# In the order of the sampes: S1 S4 S7 S10 S2 S5 S8 S11 S3 S6 S9 S12
# 56,80,65,19,69,160,205,50,101,200,296,97
POOLSIZES="8,80,65,19,69,160,205,50,101,200,296,97"

method="karlsson"
region="5:3165000-3190000"

${GRENEDALF} fst --pileup-path $MPILEUP --sample-name-list $NAMES --filter-region ${region} --window-width 1 --pool-sizes ${POOLSIZES} --method "${method}" --omit-na-windows --file-suffix "-${method}" --threads 4 > "log-${method}.log" 
        
