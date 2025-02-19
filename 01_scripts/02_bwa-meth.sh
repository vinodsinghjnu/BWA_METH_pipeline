#!/bin/bash


source $HOME/miniforge3/etc/profile.d/conda.sh
conda activate ngs_packages

SCRIPT=$0
INPUT=$1
DATANAME=$(basename $INPUT)
NCPUS="$2"
SAMPLE="$3"


NAME=$(basename $0)
LLOG_FOLDER="10_logfiles/${DATANAME}"

# Define options
GENOME="02_reference/hg38/hg38.fa"  # Genomic reference .fasta
#GENOME="02_reference/GrCh38_EMseq/grch38_core_plus_bs_controls.fa"
TRIMMED_FOLDER="04_trimmed/${DATANAME}"
ALIGNED_FOLDER="05_aligned/${DATANAME}"
mkdir $ALIGNED_FOLDER

TEMP_FOLDER="99_tmp/"



#ls -1 "$INPUT"/*_R1_001.fastq.gz | perl -pe 's/_R[12]_001.fastq\.gz//g' | parallel basename {}  > $(basename $INPUT)_samples_for_alignment.txt
SAMPLE_FILE=$(basename $INPUT)_samples_for_alignment.txt

start=`date +%s.%N`

# Align reads
cat "$SAMPLE_FILE" | grep -E "$SAMPLE" |
while read file
do
    base=$(basename $file)
    echo "Aligning $base"

    # Align
    bwameth.py --threads "$NCPUS" \
        --reference "$GENOME" \
        "$TRIMMED_FOLDER"/"$base"_R1.fastq.gz \
        "$TRIMMED_FOLDER"/"$base"_R2.fastq.gz |
        samtools view -Sb -q 10 - |
        samtools sort - > "$ALIGNED_FOLDER"/"$base".bam

    samtools index "$ALIGNED_FOLDER"/"$base".bam
done

# Cleanup temp folder
rm -r "$TEMP_FOLDER"/* 2>/dev/null

conda deactivate

end=`date +%s.%N`
runtime=$( echo "$end - $start" | bc -l )
runtime=${runtime%.*}
hours=$((runtime / 3600)); minutes=$(( (runtime % 3600) / 60 )); seconds=$(( (runtime % 3600) % 60 ))
echo "==> completed Runtime: $hours:$minutes:$seconds (hh:mm:ss)"