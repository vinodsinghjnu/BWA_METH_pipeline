#!/bin/bash


INPUT_R1="$1"
INPUT_R2="$2"
GENOME="$3"
ALIGNED_BAM="$4"
NCPUS="$5"
SAMPLE="$6"

start=`date +%s.%N`

# Check that all required variables are set
if [[ -z "$INPUT_R1" || -z "$INPUT_R2" || -z "$GENOME" || -z "$ALIGNED_BAM" || -z "$SAMPLE" ]]; then
    echo "Error: One or more required variables are unset!"
    exit 1
fi

echo "Processing: $INPUT_R1 and $INPUT_R2"
echo "Output BAM: $ALIGNED_BAM"


# Align

bwameth.py --threads "$((NCPUS / 2))" \
    --reference "$GENOME" \
    "$INPUT_R1" "$INPUT_R2" | \
    samtools view -@ "$((NCPUS / 4))" -Sb -q 10 -F 4 -F 256 -F 2048 - | \
    samtools sort -@ "$((NCPUS / 2))" -o "$ALIGNED_BAM" -T "${SAMPLE}_temp" && \
    samtools index "$ALIGNED_BAM"

## Add MD and NM tags
#samtools calmd -@ "$NCPUS" -b "$ALIGNED_BAM" "$GENOME" > "${ALIGNED_BAM%.bam}.calmd.bam" 2> "${ALIGNED_BAM%.bam}.calmd.log"

## Re-index the updated BAM
#samtools index "${ALIGNED_BAM%.bam}.calmd.bam"

## Remove original BAM and its index
#rm -v "$ALIGNED_BAM" "$ALIGNED_BAM.bai"

## Rename the new BAM to match the original filename (optional)
#mv -v "${ALIGNED_BAM%.bam}.calmd.bam" "$ALIGNED_BAM"
#mv -v "${ALIGNED_BAM%.bam}.calmd.bam.bai" "$ALIGNED_BAM.bai"

##

# samtools view -@ "$((NCPUS / 4))" -Sb -q 10 -f 1 - | \
# samtools view -@ "$((NCPUS / 4))" -Sb -q 1 -F 4 -F 256 -F 2048 -
# -q 1: Keeps nearly all mapped reads (removes unaligned ones).
# -f 1: Keeps properly paired reads
# -F 4: Excludes unmapped reads.
# -F 256: Removes secondary alignments.
# -F 2048: Removes supplementary alignments.

end=`date +%s.%N`
runtime=$( echo "$end - $start" | bc -l )
runtime=${runtime%.*}
hours=$((runtime / 3600)); minutes=$(( (runtime % 3600) / 60 )); seconds=$(( (runtime % 3600) % 60 ))
echo "==> completed Runtime: $hours:$minutes:$seconds (hh:mm:ss)"

