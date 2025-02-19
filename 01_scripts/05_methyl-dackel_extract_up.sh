#!/bin/bash


INPUT="$1"
DATANAME=$(basename $INPUT)
NCPUS="$2"
BAM_PATTERN="$3"
REAL_SAMPLE="$4"
#BAM_PATTERN="_final.bam"
#BAM_PATTERN=".dedup.filt.bam"



# Conda env
source $HOME/miniforge3/etc/profile.d/conda.sh
conda activate samtools_1.16.1

# keep some info
TIMESTAMP=$(date +%Y-%m-%d_%Hh%Mm%Ss)
SCRIPT=$0
NAME=$(basename $0)
LOG_FOLDER="10_logfiles/${DATANAME}"
#cp "$SCRIPT" "$LOG_FOLDER"/"$TIMESTAMP"_"$NAME"

# Define options
GENOME="02_reference/hg38/hg38.fa"  # Genomic reference .fasta
#GENOME="02_reference/hg38/ncbi_dataset/data/GCF_000001405.33/GCF_000001405.33_GRCh38.p7_genomic.fa"  # Genomic reference .fasta
DEDUPLICATED_FOLDER="06_deduplicated/${DATANAME}"
METHYLDACKEL_FOLDER="07_methyl_dackel/${DATANAME}"
TEMP_FOLDER="99_tmp/"

start=`date +%s.%N`


#l=$(ls -1 "$DEDUPLICATED_FOLDER"/*"$BAM_PATTERN" | grep -E "$REAL_SAMPLE" |  parallel "MethylDackel mbias -@ $NCPUS -r chr1 $GENOME {} $(basename {} _final.bam)_chr1_cpg")

ls -1 "$DEDUPLICATED_FOLDER"/*"$BAM_PATTERN" | grep -E "$REAL_SAMPLE" |  xargs -I {} sh -c "MethylDackel mbias -@ $NCPUS $GENOME {} $(basename {} ${BAM_PATTERN})_cpg>&OTOB$REAL_SAMPLE.txt"

ls -1 "$DEDUPLICATED_FOLDER"/*"$BAM_PATTERN" | grep -E "$REAL_SAMPLE" |  xargs -I {} sh -c "MethylDackel mbias -@ $NCPUS --noCpG --CHH --CHG $GENOME {} $(basename {} ${BAM_PATTERN})_chn>&OTOBchn$REAL_SAMPLE.txt"

#l=$(cat OTOB$REAL_SAMPLE.txt | cut -d ':' -f 2)


# Gnu Parallel
ls -1 "$DEDUPLICATED_FOLDER"/*"$BAM_PATTERN" | grep -E "$REAL_SAMPLE" |
parallel -v -j "$NCPUS" \
    MethylDackel extract $(cat OTOB$REAL_SAMPLE.txt | cut -d ':' -f 2) --minOppositeDepth 4 --maxVariantFrac 0.5 "$GENOME" {} \; \
    gzip {.}_CpG.bedGraph \; mv {.}_CpG.bedGraph.gz "$METHYLDACKEL_FOLDER" \; \
    MethylDackel extract $(cat OTOB$REAL_SAMPLE.txt | cut -d ':' -f 2) --minOppositeDepth 4 --maxVariantFrac 0.5 --methylKit "$GENOME" {} \; \
    gzip {.}_CpG.methylKit \; mv {.}_CpG.methylKit.gz "$METHYLDACKEL_FOLDER" \; \
    MethylDackel extract $(cat OTOB$REAL_SAMPLE.txt | cut -d ':' -f 2) --minOppositeDepth 4 --maxVariantFrac 0.5 --mergeContext "$GENOME" {} \; \
    mv {.}_CpG.bedGraph {.}_CpG_merged.bedGraph \; gzip {.}_CpG_merged.bedGraph \; mv {.}_CpG_merged.bedGraph.gz "$METHYLDACKEL_FOLDER" \; \
    MethylDackel extract $(cat OTOBchn$REAL_SAMPLE.txt | cut -d ':' -f 2) --noCpG --CHH --CHG --minOppositeDepth 4 --maxVariantFrac 0.5 "$GENOME" {} \; \
    gzip "{.}_CH*.bedGraph" \; mv "{.}_CH*.bedGraph.gz" "$METHYLDACKEL_FOLDER" \; \
    MethylDackel extract $(cat OTOBchn$REAL_SAMPLE.txt | cut -d ':' -f 2) --noCpG --CHH --CHG --minOppositeDepth 4 --maxVariantFrac 0.5 --methylKit "$GENOME" {} \; \
    gzip "{.}_CH*.methylKit" \; mv "{.}_CH*.methylKit.gz" "$METHYLDACKEL_FOLDER" \; \
    

# Move files to METHYLDACKEL_FOLDER
#ls -1 "$DEDUPLICATED_FOLDER"/ | grep -E "$REAL_SAMPLE" | grep -v \.bam | parallel -v -j "$NCPUS" mv "$DEDUPLICATED_FOLDER"/{} "$METHYLDACKEL_FOLDER" # mv 

ls OTOB*.txt | grep -E "$REAL_SAMPLE" | parallel -j "$NCPUS" mv {} "$METHYLDACKEL_FOLDER"
ls "$DEDUPLICATED_FOLDER"/*.svg | grep -E "$REAL_SAMPLE" | parallel -j "$NCPUS" mv {} "$METHYLDACKEL_FOLDER"


conda deactivate



end=`date +%s.%N`
runtime=$( echo "$end - $start" | bc -l )
runtime=${runtime%.*}

hours=$((runtime / 3600)); minutes=$(( (runtime % 3600) / 60 )); seconds=$(( (runtime % 3600) % 60 ))
echo "==> completed Runtime: $hours:$minutes:$seconds (hh:mm:ss)"
