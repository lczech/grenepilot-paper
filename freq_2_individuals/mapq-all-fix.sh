#!/bin/bash

# moinode
GRENEDALF=/home/lczech/grenedalf/bin/grenedalf
DATADIR="/Carnegie/DPB/Data/Shared/Labs/Moi/Everyone/ath_evo/grenepilot_lucas/freq_2_individuals/mapped"
OUTDIR="/Carnegie/DPB/Data/Shared/Labs/Moi/Everyone/ath_evo/grenepilot_lucas/freq_2_individuals/data-all"

SETS=()
NAMS=()

SETS+=("-q 0")
SETS+=("-Q 30 -q 60 --rf 0x002 --ff 0x004 --ff 0x008")

NAMS+=("q0")
NAMS+=("Q30-q60-f")

run_tools() {

    {
    BAMFILE=$1
    INDEX=$2
    BASENAME="${BAMFILE%.*}"
    BASENAME="${BASENAME%.*}"

    SETTINGS=${SETS[${INDEX}]}
    NAME=${NAMS[${INDEX}]}

    echo "${BASENAME} : ${NAME} with ${SETTINGS}"

    echo "fix samtools: `date`"
    cd ${OUTDIR}
    zcat tabtab/${BASENAME}-${NAME}.mpileup.gz | \
        sed "s/\t\t/\t*\t*/g" | \
        gzip > ${OUTDIR}/${BASENAME}-${NAME}.mpileup.gz

    echo "grenedalf: `date`"
    cd ${OUTDIR}
    $GRENEDALF frequency \
        --pileup-file "${BASENAME}-${NAME}.mpileup.gz" \
        --omit-invariants \
        --file-prefix "${BASENAME}-${NAME}-" \
        --sample-name-prefix "S" \
        >> grenedalf-${BASENAME}-${NAME}.log

    } &
}

export -f run_tools

#INFILES="TWOFH1_CKDL210018053-1a-AK31398-AK8965_HGKVHCCX2_L1-1.sorted.bam"

# Waaaaaahh memex does not have `parallel`, such a basic tool...
# So, use a different approach instead. Ugly, but gets the job done.
#parallel run_tools ::: ${!NAMS[@]}
cd ${DATADIR}
#for j in ${INFILES} ; do
for j in `ls *.sorted.bam` ; do
    for i in ${!NAMS[@]} ; do
        echo $j : $i
        run_tools $j $i
    done
    #wait
done
wait

echo "Done `date`"

