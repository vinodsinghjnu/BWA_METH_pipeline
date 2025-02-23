#!/bin/bash


ALIGNED_BAM="$1"
DEDUP_BAM="$2"
GENOME="$3"
NCPUS="$4"
SAMPLE="$5"


#samtools fixmate -m "$ALIGNED_BAM" "${ALIGNED_BAM%.bam}.fixmate.bam"
#samtools sort  -@ "$NCPUS" "${ALIGNED_BAM%.bam}.fixmate.bam" -o "${ALIGNED_BAM%.bam}.sorted_fixmate.bam"
#mv "${ALIGNED_BAM%.bam}.sorted_fixmate.bam" "$ALIGNED_BAM"
#samtools markdup -r -@ "$NCPUS" "$ALIGNED_BAM" "$DEDUP_BAM"
#samtools index "$DEDUP_BAM"


start=`date +%s.%N`


## Add NM tag to BAM using picard
NMtag_BAM="${ALIGNED_BAM%.bam}.picardNMTag.bam"

addMNtag_toBam_command="picard SetNmMdAndUqTags \
    -I \"${ALIGNED_BAM}\" \
    -O \"${NMtag_BAM}\" \
    -R \"$GENOME\" \
    --IS_BISULFITE_SEQUENCE true \
    --CREATE_INDEX true"

echo ">> Executing SetNmMdAndUqT command: $addMNtag_toBam_command"
eval $addMNtag_toBam_command

# validate BAM ##
VALIDATE_REPORT="${ALIGNED_BAM%.bam}.validateReport.picardNMTag.txt"
validate_command="picard ValidateSamFile \
    -I \"${NMtag_BAM}\" \
    -MODE SUMMARY \
    -BISULFITE true \
    -R \"$GENOME\" \
    -O \"$VALIDATE_REPORT\""

echo ">> Executing validate command: $validate_command"
eval $validate_command


# Mark Duplicates
METRICS_FILE="${DEDUP_BAM%.bam}.metrics.txt"
markDuplicate_command="picard MarkDuplicates \
    INPUT=\"$NMtag_BAM\" \
    OUTPUT=\"$DEDUP_BAM\" \
    METRICS_FILE=\"$METRICS_FILE\" \
    VALIDATION_STRINGENCY=SILENT \
    MAX_FILE_HANDLES=100 \
    REMOVE_DUPLICATES=true"

echo ">> Executing MarkDuplicates command: $markDuplicate_command"
eval $markDuplicate_command



# INDEXING ##
echo "> Indexing"
samtools index  "$DEDUP_BAM"


# sambamba this is faster to markduplicates and multithreaded, but Picard is wideli used
# sambamba markdup -l 1 -t 16 --sort-buffer-size 16000 --overflow-list-size 10000000 

end=`date +%s.%N`
runtime=$( echo "$end - $start" | bc -l )
runtime=${runtime%.*}

hours=$((runtime / 3600)); minutes=$(( (runtime % 3600) / 60 )); seconds=$(( (runtime % 3600) % 60 ))
echo "==> completed Runtime: $hours:$minutes:$seconds (hh:mm:ss)"