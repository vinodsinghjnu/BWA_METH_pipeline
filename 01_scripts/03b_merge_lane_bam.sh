#!/bin/bash


INPUT="$1"
DATANAME=$(basename $INPUT)
NCPUS="$2"
REAL_SAMPLE="$3"
SAMPLE_CLASS="$4"


source $HOME/miniforge3/etc/profile.d/conda.sh
conda  activate samtools_1.16.1

export LC_CTYPE=en_US.UTF-8
export LC_ALL=en_US.UTF-8


#Global variables
#MARKDUPS="/home/vinod/Software/picard-tools-1.119/MarkDuplicates.jar"
DEDUPFOLDER="06_deduplicated/${DATANAME}"
METRICSFOLDER="98_metrics"
GENOME="02_reference/hg38/hg38.fa"
#GENOME="02_reference/GrCh38_EMseq/grch38_core_plus_bs_controls.fa"

# Copy script to log folder
TIMESTAMP=$(date +%Y-%m-%d_%Hh%Mm%Ss)
SCRIPT=$0
NAME=$(basename $0)
LOG_FOLDER="10_logfiles/${DATANAME}"
#cp "$SCRIPT" "$LOG_FOLDER"/"$TIMESTAMP"_"$NAME"

start=`date +%s.%N`



# Merge BAM files ##
echo -e  "\n\n==Merger BAM files == \n"
echo -e "==BAM MERGE:==\nsamtools merge -o ${DEDUPFOLDER}/${REAL_SAMPLE}_Sperm_${SAMPLE_CLASS}_${REAL_SAMPLE}_final.bam -@ $NCPUS $(ls -1 ${DEDUPFOLDER}/*${REAL_SAMPLE}*dedup.filt.bam)"
samtools merge -o "$DEDUPFOLDER"/"$REAL_SAMPLE"_Sperm_"$SAMPLE_CLASS"_"$REAL_SAMPLE"_final.bam -@ "$NCPUS" $(ls -1 "$DEDUPFOLDER"/*$REAL_SAMPLE*dedup.filt.bam)

## Varify read counts after merge ##
echo -e  "\n\n==Varify merged files read counts == \n"
reads_before=$(ls -1 ${DEDUPFOLDER}/*${REAL_SAMPLE}*dedup.filt.bam | parallel -J 6 samtools flagstat  {} | grep -E total | cut -d ' '  -f1 | awk '{ sum += $1 } END { print sum }')
reads_after_merge=$(samtools flagstat "$DEDUPFOLDER"/"$REAL_SAMPLE"_Sperm_"$SAMPLE_CLASS"_"$REAL_SAMPLE"_final.bam | grep -E total | cut -d ' '  -f1)

if [ "$reads_before" == "$reads_after_merge" ]; then
    echo "==> successful"
else
    echo '==> Error'
fi



## Indexing #
echo -e  "\n\n==Index Merged BAM == \n"

echo -e "==Indexing:==\nsamtools index -@ $NCPUS  ${DEDUPFOLDER}/${REAL_SAMPLE}_Sperm_${SAMPLE_CLASS}_${REAL_SAMPLE}_final.bam"
samtools index -@ "$NCPUS"  "$DEDUPFOLDER"/"$REAL_SAMPLE"_Sperm_"$SAMPLE_CLASS"_"$REAL_SAMPLE"_final.bam


# validate BAM ##

echo -e  "\n\n==validate BAM after merging == \n"

ls -1 "$DEDUPFOLDER"/*_final.bam | grep -E "$REAL_SAMPLE" |
parallel -j "$NCPUS" \
    picard ValidateSamFile \
    -I {} \
    -MODE SUMMARY \
    -BISULFITE true \
    --VALIDATE_INDEX true \
    -R $GENOME \; echo


conda deactivate



end=`date +%s.%N`
runtime=$( echo "$end - $start" | bc -l )
runtime=${runtime%.*}

hours=$((runtime / 3600)); minutes=$(( (runtime % 3600) / 60 )); seconds=$(( (runtime % 3600) % 60 ))

echo "==> completed Runtime: $hours:$minutes:$seconds (hh:mm:ss)"