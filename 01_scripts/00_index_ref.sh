#!/bin/bash

# 1 CPU
# 12 Go

# Copy script to log folder
TIMESTAMP=$(date +%Y-%m-%d_%Hh%Mm%Ss)
SCRIPT=$0
NAME=$(basename $0)
LOG_FOLDER="10_logfiles"

# Global variables
REF="02_reference/hg38/hg38.fa"
#INDREF="02_reference"
#REF="02_reference/GrCh38_EMseq/grch38_core_plus_bs_controls.fa"
INDREF="02_reference/GrCh38_EMseq"


cd bwa-meth-master
export PATH=${PATH}:$PWD
cd ..

cd toolshed-0.4.0
export PATH=${PATH}:$PWD
cd ..



echo "$SCRIPT"
cp "$SCRIPT" "$LOG_FOLDER"/"$TIMESTAMP"_"$NAME"

# Index the reference with BSseeker2
bwameth.py index "$REF"
samtools faidx "$REF"
