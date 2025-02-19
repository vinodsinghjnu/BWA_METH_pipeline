#!/bin/bash


INPUT="$1"
DATANAME=$(basename $INPUT)
NCPUS="$2"
REAL_SAMPLE="$3"
BAM_PATTERN="_final.bam"


source $HOME/miniforge3/etc/profile.d/conda.sh
conda  activate samtools_1.16.1

export LC_CTYPE=en_US.UTF-8
export LC_ALL=en_US.UTF-8


#Global variables
#MARKDUPS="/home/vinod/Software/picard-tools-1.119/MarkDuplicates.jar"
DEDUPFOLDER="06_deduplicated/${DATANAME}"
METRICSFOLDER="98_metrics"
GENOME="02_reference/hg38/hg38.fa"
METHYLDACKEL_FOLDER="07_methyl_dackel/${DATANAME}"


#bamFile="06_deduplicated/Muenster/179956_Sperm_UndiffSpermatogonia_179956_final.bam" 
bamFile=$(ls -1 ${DEDUPFOLDER}/*${BAM_PATTERN} | grep -E $REAL_SAMPLE)
sample=$(basename "$bamFile" _final.bam)


echo -e "chr\tcontext\tstrand\tRead\tPosition\tnMethylated\tnUnmethylated\tnMethylated(+dups)\tnUnmethylated(+dups)" > "${METHYLDACKEL_FOLDER}"/"${sample}_combined_mbias.tsv"
chrs=($(samtools view -H "$bamFile" | grep @SQ | cut -f 2 | sed 's/SN://'| grep -v _random | grep -v "chrUn"  | sed 's/|/\\|/'))

    for chr in ${chrs[*]}; do
        for context in CHH CHG CpG; do
            arg=''
            if [ $context = 'CHH' ]; then
                arg='--CHH --noCpG'
            elif [ $context = 'CHG' ]; then
                arg='--CHG --noCpG'
            fi
            # need two calls to add columns containing the counts without filtering duplicate reads (for rrEM-seq where start/end is constrained)
            # not sure why we need both --keepDupes and -F, probably a bug in mbias
            join -t $'\t' -j1 -o 1.2,1.3,1.4,1.5,1.6,2.5,2.6 -a 1 -e 0 \
            <( \
                MethylDackel mbias --noSVG $arg -@ "$NCPUS" -r "$chr" "$GENOME" "$bamFile" | \
                tail -n +2 | awk '{print $1"-"$2"-"$3"\t"$0}' | sort -k 1b,1
            ) \
            <( \
                MethylDackel mbias --noSVG -@ "$NCPUS" --keepDupes -F 2816 $arg -@ "$NCPUS" -r "$chr" "$GENOME" "$bamFile" | \
                tail -n +2 | awk '{print $1"-"$2"-"$3"\t"$0}' | sort -k 1b,1
            ) \
            | sed "s/^/${chr}\t${context}\t/" \
            >> "${METHYLDACKEL_FOLDER}"/"${sample}_combined_mbias.tsv"
        done
    done

    # makes the svg files for trimming checks
    MethylDackel mbias -@ "$NCPUS" --noCpG --CHH --CHG -r ${chrs[0]} "$GENOME" "$bamFile" "${sample}_chn"
    for f in *chn*.svg; do sed -i "s/Strand<\\/text>/Strand $f ${chrs[0]} CHN <\\/text>/" "$f"; done;

    MethylDackel mbias -@ "$NCPUS" -r ${chrs[0]} "$GENOME" "$bamFile" "${sample}_cpg"
    for f in *cpg*.svg; do sed -i "s/Strand<\\/text>/Strand $f ${chrs[0]} CpG<\\/text>/" "$f"; done;

ls -1 | grep -E  "*_OB.svg|*_OT.svg" | parallel -j "$NCPUS" mv "$DEDUPLICATED_FOLDER"/{} "$METHYLDACKEL_FOLDER" # mv 