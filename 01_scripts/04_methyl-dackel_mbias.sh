#!/bin/bash


INPUT="$1"
DATANAME=$(basename $INPUT)
NCPUS="$2"
GENOME="$3"
BAM="$4"




# conda env
source $HOME/miniforge3/etc/profile.d/conda.sh
#conda deactivate
conda activate samtools_1.16.1
#conda activate ngs_packages

SCRIPT=$0
NAME=$(basename $0)
LOG_FOLDER="10_logfiles/${DATANAME}"

# Define options
#GENOME="02_reference/hg38/hg38.fa"  # Genomic reference .fasta
#GENOME="02_reference/hg38/ncbi_dataset/data/GCF_000001405.33/GCF_000001405.33_GRCh38.p7_genomic.fa"  # Genomic reference .fasta
DEDUPLICATED_FOLDER="06_deduplicated/${DATANAME}/"
METHYLDACKEL_FOLDER="07_methyl_dackel/${DATANAME}/"

start=`date +%s.%N`

# Gnu Parallel
#ls -1 "$DEDUPLICATED_FOLDER"/*"$BAM_PATTERN" | grep -E "$REAL_SAMPLE" |
#parallel -j "$NCPUS" \
#    MethylDackel mbias "$GENOME" {} {.}_mbias

echo -e  "\n==Command == \n"
MethylDackel mbias -@ $NCPUS $GENOME {} {.}_mbias
parallel -v MethylDackel mbias -@ "$NCPUS" "$GENOME" {} {.}_mbias && mv {.}_mbias.svg  "$METHYLDACKEL_FOLDER" ::: $BAM



conda deactivate


end=`date +%s.%N`
runtime=$( echo "$end - $start" | bc -l )
runtime=${runtime%.*}

hours=$((runtime / 3600)); minutes=$(( (runtime % 3600) / 60 )); seconds=$(( (runtime % 3600) % 60 ))

echo "==> completed Runtime: $hours:$minutes:$seconds (hh:mm:ss)"




