Seed freqs were extracted using grenedalf by running `freq.sh`,
on the seed pileup data at /lustre/scratch/lczech/grenepipe-runs/ath-evo-seeds-pileups/
The resultingfile is called frequency.csv, and was edited to contain a combined chrom_pos column
(also done in the shell script).

The 1001g frequencies were created by Moi using plink, and are originally at ~/safedata/1001g/plink.frq
We copied this file over to here as plink.frq

The 231g frequencies were also created by Moi using plink, and are originally at ~/safedata/ath_evo/grenephase1/data-big
We copied them over here as 231g.frq

Then, using the plot.R the two scatter plots were created.
