#!/bin/bash

ml SAMtools/1.9
ml cmake gnu7

GRENEDALF=/Carnegie/DPB/Homes/Users/lczech/grenedalf/bin/grenedalf

DATADIR="/lustre/scratch/lczech/grenepipe-runs/2021-08-26-ath-grenenet-release-10/mapped"
#OUTDIR="/lustre/scratch/lczech/grenepipe-runs/2021-08-26-ath-grenenet-release-10-mapq-60"
OUTDIR="/Carnegie/DPB/Data/Shared/Labs/Moi/Everyone/ath_evo/grenepilot_lucas/freq_2_individuals/data"

# Not good: produces individual files. We want a merged file with samples in it.
#cd $DATADIR
#for f in `ls *.sorted.bam` ; do
#    echo "At $f"
#    samtools mpileup -q 60 $f | gzip > ${DIR60}/${f/.sorted.bam/.mpileup.gz}
#done

# Merge the files first, in the same order as the list in all-merged-units.names.txt
#cd $DATADIR
#samtools merge -O BAM all-merged.bam \
#    COL0L1B_CKDL210018053-1a-AK10689-AK9012_HGKVHCCX2_L1-1.sorted.bam \
#    COLRUM_POOL1_CKDL210018053-1a-AK34070-AK9050_HGKVHCCX2_L1-1.sorted.bam \
#    COLRUM_POOL2_CKDL210018053-1a-7UDI1539-AK8191_HGKVHCCX2_L1-1.sorted.bam \
#    COLRUM_POOL3_CKDL210018053-1a-7UDI1971-AK20024_HGKVHCCX2_L1-1.sorted.bam \
#    GRENE231MM1_CKDL210018053-1a-AK33987-AK9050_HGKVHCCX2_L1-1.sorted.bam \
#    GRENE231MM4_CKDL210018053-1a-AK28238-AK8191_HGKVHCCX2_L1-1.sorted.bam \
#    GRENE231MM8_CKDL210018053-1a-AK9787-AK20024_HGKVHCCX2_L1-1.sorted.bam \
#    UM20L1A_CKDL210018053-1a-AK34068-AK8944_HGKVHCCX2_L1-1.sorted.bam \
#    TWOFH1_CKDL210018053-1a-AK31398-AK8965_HGKVHCCX2_L1-1.sorted.bam \
#    TWOFH3_CKDL210018053-1a-AK7838-AK24534_HGKVHCCX2_L1-1.sorted.bam \
#    TWOFH5_CKDL210018053-1a-AK34069-AK26896_HGKVHCCX2_L1-1.sorted.bam 

#samtools mpileup -q 60 -Q 15 all-merged.bam | gzip > ${OUTDIR}/all-merged-samples-q60-Q15.mpileup.gz

#samtools mpileup -q 60 -Q 15  \
#    COL0L1B_CKDL210018053-1a-AK10689-AK9012_HGKVHCCX2_L1-1.sorted.bam \
#    COLRUM_POOL1_CKDL210018053-1a-AK34070-AK9050_HGKVHCCX2_L1-1.sorted.bam \
#    COLRUM_POOL2_CKDL210018053-1a-7UDI1539-AK8191_HGKVHCCX2_L1-1.sorted.bam \
#    COLRUM_POOL3_CKDL210018053-1a-7UDI1971-AK20024_HGKVHCCX2_L1-1.sorted.bam \
#    GRENE231MM1_CKDL210018053-1a-AK33987-AK9050_HGKVHCCX2_L1-1.sorted.bam \
#    GRENE231MM4_CKDL210018053-1a-AK28238-AK8191_HGKVHCCX2_L1-1.sorted.bam \
#    GRENE231MM8_CKDL210018053-1a-AK9787-AK20024_HGKVHCCX2_L1-1.sorted.bam \
#    UM20L1A_CKDL210018053-1a-AK34068-AK8944_HGKVHCCX2_L1-1.sorted.bam \
#    TWOFH1_CKDL210018053-1a-AK31398-AK8965_HGKVHCCX2_L1-1.sorted.bam \
#    TWOFH3_CKDL210018053-1a-AK7838-AK24534_HGKVHCCX2_L1-1.sorted.bam \
#    TWOFH5_CKDL210018053-1a-AK34069-AK26896_HGKVHCCX2_L1-1.sorted.bam | \
#    gzip > ${OUTDIR}/all-merged-units-q60-Q15.mpileup.gz


SETS=()
NAMS=()

SETS+=("-q 0")
SETS+=("-q 20")
SETS+=("-q 40")
SETS+=("-q 60")
SETS+=("-q 0 --rf 0x002 --ff 0x004 --ff 0x008")
SETS+=("-q 20 --rf 0x002 --ff 0x004 --ff 0x008")
SETS+=("-q 40 --rf 0x002 --ff 0x004 --ff 0x008")
SETS+=("-q 60 --rf 0x002 --ff 0x004 --ff 0x008")
SETS+=("-Q 30 -q 0")
SETS+=("-Q 30 -q 20")
SETS+=("-Q 30 -q 40")
SETS+=("-Q 30 -q 60")
SETS+=("-Q 30 -q 0 --rf 0x002 --ff 0x004 --ff 0x008")
SETS+=("-Q 30 -q 20 --rf 0x002 --ff 0x004 --ff 0x008")
SETS+=("-Q 30 -q 40 --rf 0x002 --ff 0x004 --ff 0x008")
SETS+=("-Q 30 -q 60 --rf 0x002 --ff 0x004 --ff 0x008")

NAMS+=("q0")
NAMS+=("q20")
NAMS+=("q40")
NAMS+=("q60")
NAMS+=("q0-f")
NAMS+=("q20-f")
NAMS+=("q40-f")
NAMS+=("q60-f")
NAMS+=("Q30-q0")
NAMS+=("Q30-q20")
NAMS+=("Q30-q40")
NAMS+=("Q30-q60")
NAMS+=("Q30-q0-f")
NAMS+=("Q30-q20-f")
NAMS+=("Q30-q40-f")
NAMS+=("Q30-q60-f")

BAMFILE="COLRUM_POOL1_CKDL210018053-1a-AK34070-AK9050_HGKVHCCX2_L1-1.sorted.bam"
REFFILE="/lustre/scratch/lczech/grenepipe-runs/ath-evo-ref/TAIR10_chr_all.fa"

run_tools() {

    {
    INDEX=$1

    SETTINGS=${SETS[${INDEX}]}
    NAME=${NAMS[${INDEX}]}

    echo "${NAME} with ${SETTINGS}"

    echo "samtools: `date`"
    cd ${DATADIR}
    samtools mpileup ${SETTINGS} --fasta-ref ${REFFILE} \
        ${BAMFILE} | \
        gzip > ${OUTDIR}/COLRUM_POOL1-${NAME}.mpileup.gz

    echo "grenedalf: `date`"
    cd ${OUTDIR}
    $GRENEDALF frequency \
        --pileup-file "COLRUM_POOL1-${NAME}.mpileup.gz" \
        --omit-invariants \
        --file-suffix "-${NAME}" \
        --sample-name-prefix "S" \
        >> grenedalf-${NAME}.log

    } &
}

export -f run_tools

# Waaaaaahh memex does not have `parallel`, such a basic tool...
# So, use a different approach instead. Ugly, but gets the job done.
#parallel run_tools ::: ${!NAMS[@]}
for i in ${!NAMS[@]} ; do
    echo $i
    run_tools $i
done
wait


echo "Done `date`"


#echo "COLRUM_POOL1-q20"
#samtools mpileup -q 20 COLRUM_POOL1_CKDL210018053-1a-AK34070-AK9050_HGKVHCCX2_L1-1.sorted.bam | \
#    gzip > ${OUTDIR}/COLRUM_POOL1-q20.mpileup.gz
#
#echo "COLRUM_POOL1-q40"
#samtools mpileup -q 40 COLRUM_POOL1_CKDL210018053-1a-AK34070-AK9050_HGKVHCCX2_L1-1.sorted.bam | \
#    gzip > ${OUTDIR}/COLRUM_POOL1-q40.mpileup.gz
#
#echo "COLRUM_POOL1-q20-fF"
#samtools mpileup -q 20 -f 0x002 -F 004 -F 0x008 COLRUM_POOL1_CKDL210018053-1a-AK34070-AK9050_HGKVHCCX2_L1-1.sorted.bam | \
#    gzip > ${OUTDIR}/COLRUM_POOL1-q20-fF.mpileup.gz
#
#echo "COLRUM_POOL1-q40-fF"
#samtools mpileup -q 40 -f 0x002 -F 004 -F 0x008 COLRUM_POOL1_CKDL210018053-1a-AK34070-AK9050_HGKVHCCX2_L1-1.sorted.bam | \
#    gzip > ${OUTDIR}/COLRUM_POOL1-q40-fF.mpileup.gz
#
#cd ${OUTDIR}
#
#echo "COLRUM_POOL1-q20"
#$GRENEDALF frequency --pileup-file "COLRUM_POOL1-q20.mpileup.gz" --omit-invariants --file-suffix "-q20" --sample-name-prefix "S" >> grenedalf.log
#
#echo "COLRUM_POOL1-q40"
#$GRENEDALF frequency --pileup-file "COLRUM_POOL1-q40.mpileup.gz" --omit-invariants --file-suffix "-q40" --sample-name-prefix "S" >> grenedalf.log
#
#echo "COLRUM_POOL1-q20-fF"
#$GRENEDALF frequency --pileup-file "COLRUM_POOL1-q20-fF.mpileup.gz" --omit-invariants --file-suffix "-q20-fF" --sample-name-prefix "S" >> grenedalf.log
#
#echo "COLRUM_POOL1-q40-fF"
#$GRENEDALF frequency --pileup-file "COLRUM_POOL1-q40-fF.mpileup.gz" --omit-invariants --file-suffix "-q40-fF" --sample-name-prefix"S" >> grenedalf.log

