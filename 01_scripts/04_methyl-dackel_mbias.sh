#!/bin/bash


DEDUP_BAM="$1"
MBIAS_OUT="$2"
GENOME="$3"
NCPUS="$4"
SAMPLE="$5"


start=`date +%s.%N`

# Gnu Parallel

echo -e  "\n==Command == \n"
MethylDackel mbias -@ "$NCPUS" "$GENOME" "$DEDUP_BAM" "$MBIAS_OUT"

# Get Alignment trimming parameters in "${MBIAS_OUT%.svg}_cpgOTOB.txt"
MethylDackel mbias -@ "$NCPUS" "$GENOME" "$DEDUP_BAM" "${MBIAS_OUT%.svg}_cpg" >& "${MBIAS_OUT%.svg}_cpgOTOB.txt"
#MethylDackel mbias -@ "$NCPUS" --noCpG --CHH --CHG "$GENOME" "$DEDUP_BAM" "${MBIAS_OUT%.svg}_chn">& "${MBIAS_OUT%.svg}_chnOTOB.txt"


if [[ ! -s "${MBIAS_OUT%.svg}_cpgOTOB.txt" ]]; then
    echo "Error: Input file ${MBIAS_OUT%.svg}_cpgOTOB.txt is missing or empty."
    exit 1
fi

## Create perCPG Methyl matrix in methylKit format
# MethylDackel extract deafult : MAPQ >= 10 and Phred >= 5, can changed with parameters "-q 10 -p 5"
MethylDackel extract $(cat "${MBIAS_OUT%.svg}_cpgOTOB.txt" | cut -d ':' -f 2) --minOppositeDepth 4 --maxVariantFrac 0.5 --methylKit "$GENOME" "$DEDUP_BAM" -o ${MBIAS_OUT%.svg} \
&& gzip -f ${MBIAS_OUT%.svg}_CpG.methylKit 

## Create perCPG Methyl matrix in bedGraph format, +ve and -ve CpGs at a position are merged together
MethylDackel extract $(cat "${MBIAS_OUT%.svg}_cpgOTOB.txt" | cut -d ':' -f 2) --minOppositeDepth 4 --maxVariantFrac 0.5 --mergeContext "$GENOME" "$DEDUP_BAM" -o ${MBIAS_OUT%.svg} \
&& gzip -f ${MBIAS_OUT%.svg}_CpG.bedGraph 






end=`date +%s.%N`
runtime=$( echo "$end - $start" | bc -l )
runtime=${runtime%.*}

hours=$((runtime / 3600)); minutes=$(( (runtime % 3600) / 60 )); seconds=$(( (runtime % 3600) % 60 ))

echo "==> completed Runtime: $hours:$minutes:$seconds (hh:mm:ss)"




