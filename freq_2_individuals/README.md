data created from run
~/lscratch/grenepipe-runs/2021-08-26-ath-grenenet-release-10/frequencies

on the all-merged-units.mpileup.gz file
using grenedalf frequency --omit-invariants

turn long sample names coming from grenedalf into shorter form that R does not fail with:
`sed "s/\t\([^_]*\)_[^\\.]*\\.\([^\t]*\)/\t\1\\.\2/g" frequency_long_cols.csv > frequency.csv`
