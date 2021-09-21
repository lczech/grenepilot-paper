Seed freqs were extracted using the genesis app seed_freqs,
with the command ./seed_freqs /lustre/scratch/lczech/grenepipe-runs/ath-evo-seeds-rerun/annotated/all.vcf.gz
The resulting file is called seed_freqs.csv

We first tried to extract freqs from the 1001g from scratch, but plink did not yield a file with positions
in them, so instead we used an existing file from ~/safedata/1001g/plink.frq, in the hope that this is correct.
See below for the plink commands we tried.

1000g frews were extracted with plink, using
plink --freq --vcf 1001genomes_snp-short-indel_only_ACGTN.vcf.gz
in the ~/safedata/1001g/vcf data dir
The resulting file is called plink.frq

We also tried using the plink option --biallelic-only, in order to filter out non biallelic snps.
But that somehow does not work, and instead gives a table without position information,
but still containing non-biallelic parts... No idea what this is about.

231g copied from ~/safedata/ath_evo/grenephase1/data-big
