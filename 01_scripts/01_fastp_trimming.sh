#!/bin/bash

# 4 CPU
# 10 Go


source $HOME/miniforge3/etc/profile.d/conda.sh
conda activate ngs_packages

INPUT=$1
NUMCPUS=$2
DATANAME=$3


# Copy script as it was run
TIMESTAMP=$(date +%Y-%m-%d_%Hh%Mm%Ss)
SCRIPT=$0
NAME=$(basename $0)
LOG_FOLDER="10_logfiles/${DATANAME}"
mkdir $LOG_FOLDER 
cp $SCRIPT $LOG_FOLDER/"$TIMESTAMP"_"$NAME"

# Global variables
LENGTH=60 # reads below this lengths will be removed
QUAL=25
#INPUT="03_raw_data/PRJNA357085_WGBS_rawFastq_human_spermatogonial_stem_cells_SSC"
OUTPUT="04_trimmed/${DATANAME}"
#mkdir $OUTPUT

#export LC_CTYPE=en_US.UTF-8
#export LC_ALL=en_US.UTF-8

## check if all samples have paired fastq
# ls "$INPUT"/*_R1_001.fastq.gz | perl -pe 's/R[12]\_001.fastq\.gz//g' | parallel 'echo {}; ls {}* | wc -l'

tart=`date +%s.%N`


# Trim reads with fastp and Gnu Parallel
#ls "$INPUT"/*_R1.fastq.gz | perl -pe 's/R[12]\.fastq\.gz//g' |
#ls "$INPUT"/*_R1_001.fastq.gz | perl -pe 's/R[12]\_001.fastq\.gz//g' |
ls "$INPUT"/*_R1.fastq.gz | perl -pe 's/R[12]\.fastq\.gz//g' |
parallel -v -j "$NUMCPUS" \
    fastp -i {}R1.fastq.gz -I {}R2.fastq.gz \
        -o $OUTPUT/{/}R1.fastq.gz \
        -O $OUTPUT/{/}R2.fastq.gz \
        --length_required="$LENGTH" \
        --qualified_quality_phred="$QUAL" \
        --correction \
        --trim_tail1=1 \
        --trim_tail2=1 \
        --trim_poly_x \
        --detect_adapter_for_pe \
        --overrepresentation_analysis \
        --json $OUTPUT/{/}.json \
        --html $OUTPUT/{/}.html  \
        --report_title={/}.html


# Quality check
#fastqc='/home/vinodsingh/installed_softwares/FastQC/./fastqc'
#mkdir ${OUTPUT}/fastqc_analysis

#echo "==>Quality check"

#$fastqc  -f fastq --outdir=${OUTPUT}/fastqc_analysis --memory 10000 "$INPUT"/*.fastq.gz

#multiqc fastqc_analysis/

end=`date +%s.%N`
runtime=$( echo "$end - $start" | bc -l )
runtime=${runtime%.*}

hours=$((runtime / 3600)); minutes=$(( (runtime % 3600) / 60 )); seconds=$(( (runtime % 3600) % 60 ))
echo "==> completed Runtime: $hours:$minutes:$seconds (hh:mm:ss)"

conda deactivate

### Defaults FASTP ###  # https://github.com/OpenGene/fastp#global-trimming

# trim_tail1 :  the last cycle of Illumina sequencing is uaually with low quality, and it can be 
# dropped with -t 1 or --trim_tail1=1 option.

# trim polyG (--trim_poly_g, enabled by default for NovaSeq/NextSeq data)
# --trim_poly_x (trim polyA)
# trim adapter by overlap analysis (enabled by default for PE data)