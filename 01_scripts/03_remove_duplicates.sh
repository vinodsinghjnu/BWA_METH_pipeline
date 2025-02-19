#!/bin/bash


## 1 CPU
## 100 Go

INPUT="$1"
DATANAME=$(basename $INPUT)
SAMPLE="$2"

source $HOME/miniforge3/etc/profile.d/conda.sh
conda  activate samtools_1.16.1
#conda activate ngs_packages

export LC_CTYPE=en_US.UTF-8
export LC_ALL=en_US.UTF-8



# Global variables
#MARKDUPS="/home/vinod/Software/picard-tools-1.119/MarkDuplicates.jar"
ALIGNEDFOLDER="05_aligned/${DATANAME}"
DEDUPFOLDER="06_deduplicated/${DATANAME}"
METRICSFOLDER="98_metrics"
GENOME="02_reference/hg38/hg38.fa"
#GENOME="02_reference/GrCh38_EMseq/grch38_core_plus_bs_controls.fa"

# Copy script to log folder
TIMESTAMP=$(date +%Y-%m-%d_%Hh%Mm%Ss)
SCRIPT=$0
NAME=$(basename $0)
LOG_FOLDER="10_logfiles/${DATANAME}"
#cp "$SCRIPT" "$LOG_FOLDER"/"$TIMESTAMP"_"$SCRIPT"
start=`date +%s.%N`





# validate BAM ##
ls -1 "$ALIGNEDFOLDER"/*.bam | grep -E "$SAMPLE" |
parallel -v -j 1 \
    picard ValidateSamFile \
    -I {} \
    -MODE SUMMARY \
    -BISULFITE true \
    --VALIDATE_INDEX true \
    -R $GENOME \; echo

# Remove duplicates from bam alignments with Gnu Parallel
ls -1 "$ALIGNEDFOLDER"/*.bam | grep -E "$SAMPLE" |
parallel -v -j 1 \
    picard MarkDuplicates \
    INPUT={} \
    OUTPUT="$DEDUPFOLDER"/{/.}.dedup.bam \
    METRICS_FILE="$METRICSFOLDER"/{/.}.metrics.txt \
    VALIDATION_STRINGENCY=SILENT \
    MAX_FILE_HANDLES=100 \
    REMOVE_DUPLICATES=true \; echo


# validate BAM ##
ls -1 "$DEDUPFOLDER"/*.dedup.bam | grep -E "$SAMPLE" |
parallel -v -j 1 \
    picard ValidateSamFile \
    -I {} \
    -MODE SUMMARY \
    -BISULFITE true \
    --VALIDATE_INDEX true \
    -R $GENOME \; echo

ls -1 "$DEDUPFOLDER"/*.dedup.bam | grep -E "$SAMPLE" | 
parallel -v -j 1 \
    samtools index {}

conda deactivate





# ls -1 "$ALIGNEDFOLDER"/*.bam |
# parallel -j 1  \
#     java -j 1ar "$MARKDUPS" \
#     INPUT={} \
#     OUTPUT="$DEDUPFOLDER"/{/.}.dedup.bam \
#     METRICS_FILE="$METRICSFOLDER"/{/.}.metrics.txt \
#     VALIDATION_STRINGENCY=SILENT \
#     MAX_FILE_HANDLES=100 \
#     REMOVE_DUPLICATES=true \; echo


# sambamba this is faster to markduplicates and multithreaded
# sambamba markdup -l 1 -t 16 --sort-buffer-size 16000 --overflow-list-size 10000000 

end=`date +%s.%N`
runtime=$( echo "$end - $start" | bc -l )
runtime=${runtime%.*}

hours=$((runtime / 3600)); minutes=$(( (runtime % 3600) / 60 )); seconds=$(( (runtime % 3600) % 60 ))
echo "==> completed Runtime: $hours:$minutes:$seconds (hh:mm:ss)"