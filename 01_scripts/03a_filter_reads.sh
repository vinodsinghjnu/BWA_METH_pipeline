#!/bin/bash


INPUT="$1"
DATANAME=$(basename $INPUT)
SAMPLE="$2"

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

start=`date +%s.%N`


echo ">> $SAMPLE"

# # filter out duplicates, low quality reads, not mapped in a proper pair 
ls -1  "$DEDUPFOLDER"/*.dedup.bam | grep -E "$SAMPLE" |
parallel -j 1 \
    samtools view -F 1796 -q 10  {} -o "$DEDUPFOLDER"/{/.}.filt.bam \; echo


ls -1  "$DEDUPFOLDER"/*.filt.bam | grep -E "$SAMPLE" |
parallel -j 1 \
    samtools index  {}  \; echo


# ## Validate BAM ##
echo -e  "\n\n==validate BAM after filter == \n"

ls -1 "$DEDUPFOLDER"/*.filt.bam | grep -E "$SAMPLE" |
parallel -j 1 \
    picard ValidateSamFile \
    -I {} \
    -MODE SUMMARY \
    -BISULFITE true \
    --VALIDATE_INDEX true \
    -R $GENOME \; echo

# BAM file coverage stats ##    

echo '<< BAM COVERAGE STATS <<'
bam=$(ls -1  "$DEDUPFOLDER"/*.filt.bam | grep -E "$SAMPLE")
tot_size=$(samtools view -H $bam | grep -P '^@SQ' | cut -f 3 -d ':' | awk '{sum+=$1} END {print sum}')
samtools depth $bam |  awk -v var="$tot_size" '{sum+=$3; sumsq+=$3*$3} END { print "Average = ",sum/var; print "Stdev = ",sqrt(sumsq/var - (sum/var)**2)}'
reads=$(samtools flagstat $bam | grep -E total | cut -d ' '  -f1)
echo "Reads: $reads"

conda deactivate


end=`date +%s.%N`
runtime=$( echo "$end - $start" | bc -l )
runtime=${runtime%.*}
hours=$((runtime / 3600)); minutes=$(( (runtime % 3600) / 60 )); seconds=$(( (runtime % 3600) % 60 ))
echo "==> completed Runtime: $hours:$minutes:$seconds (hh:mm:ss)"