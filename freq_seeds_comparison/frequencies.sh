#!/bin/bash

# Need these modules to be loaded, in order to run things on the cluster...
ml cmake gnu7

# Paths
GRENEDALF=/Carnegie/DPB/Homes/Users/lczech/grenedalf/bin/grenedalf
DATADIR=/lustre/scratch/lczech/grenepipe-runs

# Clean up from previous runs, just in case
rm -f frequency.csv
#rm -f grenepipe.log

# Make all frequency tables, based on filtered and unfiltered calls
LIST="ath-evo-seeds-bcftools-c ath-evo-seeds-bcftools-m ath-evo-seeds-freebayes-kv ath-evo-seeds-freebayes-nkv ath-evo-seeds-haplotypecaller-kv ath-evo-seeds-haplotypecaller-nkv"
for RUN in $LIST ; do
    echo "Run $RUN"
    for TYPE in genotyped filtered ; do
        echo "Type $TYPE"
        # Get the frequencies per position from grenedalf, using only the totals,
        # and then add a column CHROM_POS that combines both, so that we can easily merge them in R.
        #$GRENEDALF frequency --vcf-file "${DATADIR}/${RUN}/${TYPE}/all.vcf.gz" --type total >> grenepipe.log
        #cat frequency.csv | sed "s/^\([^\t]*\)\t\([^\t]*\)\t/\1\t\2\t\1_\2\t/g" > frequency-${RUN}-${TYPE}.csv
        #rm frequency.csv
    done
done

# Pileup output is a bit different - we don't have a vcf, so we need different processing
PILEUPS="all-individual-units all-merged-samples all-merged-units"
for PILEUP in $PILEUPS ; do
    echo "Run $PILEUP"
    $GRENEDALF frequency --pileup-file "${DATADIR}/ath-evo-seeds-pileups/mpileup/${PILEUP}.mpileup.gz" --type total --omit-invariants >> grenepipe.log
    cat frequency.csv | sed "s/^\([^\t]*\)\t\([^\t]*\)\t/\1\t\2\t\1_\2\t/g" > frequency-ath-evo-seeds-pileups-${PILEUP}.csv
    rm frequency.csv
done
