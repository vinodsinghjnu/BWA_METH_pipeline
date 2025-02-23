#!/bin/bash


wd='/home/vinodsingh/mutationalscanning/Workspaces/vinod/codes/bwa-meth_pipeline-master'
cd $wd

source $HOME/miniforge3/etc/profile.d/conda.sh
conda activate ngs_packages
#conda activate ngs_packages

# SLURM cluster account name #
HPC_accountName='mutationalscanning'

## required softwares #
fastqc='/home/vinodsingh/installed_softwares/FastQC/./fastqc'

# dataset name"
DATANAME="PRJEB34432_NormalSperm"


# Make important directories
mkdir -p 03_raw_data/${DATANAME}
mkdir -p 04_trimmed/${DATANAME}
mkdir -p 05_aligned/${DATANAME}
mkdir -p 06_deduplicated/${DATANAME}
mkdir -p 07_methyl_dackel/${DATANAME}
mkdir -p 10_logfiles/${DATANAME}

mkdir -p 02_reference/hg38 && wget -c https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz -O 02_reference/hg38/hg38.fa.gz && gunzip 02_reference/hg38/hg38.fa.gz
#wget -c https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz
#mkdir 02_reference/hg38/ ; gunzip hg38.fa.gz;  mv hg38.fa.gz 02_reference/hg38/
REF="${BASEDIR}/02_reference/hg38/hg38.fa"
BASEDIR='/home/vinodsingh/mutationalscanning/Workspaces/vinod/BWA_METH_pipeline/'
CODESDIR="${BASEDIR}/01_scripts/"
FASTQDIR="${BASEDIR}"/03_raw_data/"$DATANAME"/

wd=$(BASEDIR)

mkdir 04_trimmed/${DATANAME}; mkdir 05_aligned/${DATANAME}; mkdir 06_deduplicated/${DATANAME}; mkdir 07_methyl_dackel/${DATANAME}; mkdir 10_logfiles/${DATANAME}

#DATADIR='/home/vinodsingh/mutationalscanning/Workspaces/vinod/codes/bwa-meth_pipeline-master/03_raw_data/MuensterMethyFastq/'

# List Samples for alignment #
# sample name format: ERR3523436_1.fastq.gz, ERR3523436_2.fastq.gz
fastqPat=".fastq.gz"; grepPat="_R1${fastqPat}|_1${fastqPat}"
  ls "$FASTQDIR"/*$fastqPat | grep -E $grepPat | parallel basename {} | \
  awk -v pat="_$grepPat" -v repl="" '{ sub(pat,repl,$0); print $0 }' \
  > $(basename $FASTQDIR)_samples_for_alignment.txt  # for muenster published data

cat $(basename $FASTQDIR)_samples_for_alignment.txt


# # Get read length 
# readLengthFile="$FASTQDIR"/read_length.txt
# cd $FASTQDIR
# for f in *.fastq.gz; do
#   rl=$(zcat $f | awk '{if(NR%4==2) print length($1)}' | sort -n | uniq -c)
#   echo "$f $rl"
#   #echo $(zcat $f|wc -l)/4|bc
# done >$readLengthFile
# cd $wd


### Create QC report ###

mkdir ${FASTQDIR}/fastqc_analysis
srun --account $HPC_accountName --time=02:00:00 -c $n_cpus --mem 40G \
  conda run -n ngs_packages \
  fastqc  -f fastq --outdir=${FASTQDIR}/fastqc_analysis -t $n_cpus --memory 10000 "$FASTQDIR"/*.fastq.gz

# MULTIQC report
conda run -n multiqc_env multiqc ${FASTQDIR}/fastqc_analysis -o ${FASTQDIR}/fastqc_analysis --force

### Reads Triming using FASTp ###

LENGTH=60 # reads below this lengths will be removed
QUAL=25 # reads below this lengths will be removed
n_cpus=4

srun --account $HPC_accountName --time=03:00:00 -c $n_cpus --mem 40G \
  conda run -n ngs_packages \
  bash -c "echo 'Using environment: ' && conda info --envs && ls '$FASTQDIR'/*_1.fastq.gz | perl -pe 's/_[12]\.fastq\.gz//g' | \
  parallel -v -j '$n_cpus' \
      fastp -i {}_1.fastq.gz -I {}_2.fastq.gz \
          -o ${BASEDIR}/04_trimmed/${DATANAME}/{/}_1.fastq.gz \
          -O ${BASEDIR}/04_trimmed/${DATANAME}/{/}_2.fastq.gz \
          --length_required='$LENGTH' \
          --qualified_quality_phred='$QUAL' \
          --correction \
          --trim_tail1=1 \
          --trim_tail2=1 \
          --trim_poly_x \
          --detect_adapter_for_pe \
          --overrepresentation_analysis \
          --json ${BASEDIR}/04_trimmed/${DATANAME}/{/}.json \
          --html ${BASEDIR}/04_trimmed/{/}.html  \
          --report_title=${BASEDIR}/04_trimmed/${DATANAME}/{/}.html"

### Create QC report of trimmed reads ###

mkdir ${BASEDIR}/04_trimmed/${DATANAME}/fastqc_analysis
srun --account $HPC_accountName --time=02:00:00 -c $n_cpus --mem 40G \
  conda run -n ngs_packages \
  fastqc  -f fastq --outdir=${BASEDIR}/04_trimmed/${DATANAME}/fastqc_analysis -t $n_cpus --memory 10000 ${BASEDIR}/04_trimmed/${DATANAME}/*.fastq.gz

# MULTIQC report
conda run -n multiqc_env multiqc ${BASEDIR}/04_trimmed/${DATANAME}/fastqc_analysis -o ${BASEDIR}/04_trimmed/${DATANAME}/fastqc_analysis --force



###########################
#### MWA-METH alignment ### 
###########################

### index the reference with BSseeker2 ###

srun --account $HPC_accountName --time=03:00:00 -c $n_cpus --mem 40G \
  conda run -n ngs_packages \
  bash -c "echo 'Using environment: ' && echo \$CONDA_DEFAULT_ENV && bwameth.py index $REF && samtools faidx $REF"

#REF=02_reference/GrCh38_EMseq/grch38_core_plus_bs_controls.fa

## Alignment BWA-METH ##
n_cpus=24
parallel -v sbatch \
  --account $HPC_accountName \
  -J 02_bwa-meth.{} \
  --time 20:00:00 \
  -c $n_cpus \
  --mem 80G \
  ${CODESDIR}/runOn_masterCluster.sh ngs_packages \
  01_scripts/02_bwa-meth.sh \
  ${BASEDIR}/04_trimmed/${DATANAME}/{}_1.fastq.gz \
  ${BASEDIR}/04_trimmed/${DATANAME}/{}_2.fastq.gz \
  $REF \
  ${BASEDIR}/05_aligned/${DATANAME}/{}.bam \
  $n_cpus \
  {} ::::  $(basename $FASTQDIR)_samples_for_alignment.txt


# varify run is complete or not ##
grep -E error $(ls -1 xx-02_bwa-meth*.out) | wc -l
grep -E "==> completed" $(ls -1  xx-02_bwa-meth*.out) | wc -l

mv $(ls -1 xx-02_bwa-meth*.*) 10_logfiles/${DATANAME}/


## Remove duplicates ##

n_cpus=3
parallel sbatch \
  --account $HPC_accountName \
  -J 03_remove_duplicates_{} \
  --time 10:00:00 \
  --cpus-per-task $n_cpus \
  --mem=60GB \
  ${CODESDIR}/runOn_masterCluster.sh samtools_1.16.1 \
  ${CODESDIR}/03_remove_duplicates.sh  \
  ${BASEDIR}/05_aligned/${DATANAME}/{}.bam \
  ${BASEDIR}/06_deduplicated/${DATANAME}/{}.dedup.bam \
  $REF \
  $n_cpus \
  {} :::: $(basename $FASTQDIR)_samples_for_alignment.txt


grep -E "==> completed" $(ls  xx-03_remove_duplicates_*.out)
grep -E error $(ls -1 xx-03_remove_duplicates_*.out) | wc -l
grep -E "ERROR:" $(ls  xx-03_remove_duplicates_*.out)
mv $(ls -1 xx-03_remove_duplicates_*.out) 10_logfiles/${DATANAME}/



## RUN methyl-dackel_mbias plot ##

/Users/au734493/mount/genomedk_mutScanning/BWA_METH_pipeline/01_scripts/04_methyl-dackel_mbias.sh

n_cpus=10
parallel sbatch \
  --account $HPC_accountName \
  -J 04_methyl-dackel_mbias.{} \
  --time 10:00:00 \
  --cpus-per-task $n_cpus \
  --mem=60GB \
  ${CODESDIR}/runOn_masterCluster.sh samtools_1.16.1 \
  ${CODESDIR}/04_methyl-dackel_mbias.sh  \
  ${BASEDIR}/06_deduplicated/${DATANAME}/{}.dedup.bam \
  ${BASEDIR}/07_methyl_dackel/${DATANAME}/{}.dedup.mbias.svg \
  $REF \
  $n_cpus \
  {} :::: $(basename $FASTQDIR)_samples_for_alignment.txt


grep -E "==> completed" $(ls  xx-04_methyl-dackel_mbias.*.out)
grep -E error $(ls -1 xx-04_methyl-dackel_mbias.*.out) | wc -l
grep -E "ERROR:" $(ls  xx-04_methyl-dackel_mbias.*.out)
mv $(ls -1 xx-04_methyl-dackel_mbias.*.out) 10_logfiles/${DATANAME}/


